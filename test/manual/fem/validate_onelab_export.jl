# Explicit physical acceptance experiments; excluded from automatic discovery.
# Runs real FEM solves and detached exports. Usage:
# julia --project=test \
#   test/manual/fem/validate_onelab_export.jl /new/evidence/directory
# The maintained CI parity tests live in test/extensions/fem_export.jl.
using LineCableModels,Gmsh,LinearAlgebra,SHA,Printf
FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
getdp=FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
root=abspath(only(ARGS))
ispath(root) && error("Choose a new validation directory")
mkpath(root)
settings=(pml_layers=8,solver_threads=1,frequency_workers=1,keep_run_directory=true,trace=true,gmsh_verbosity=0,getdp_verbosity=0)
function matrixfile(path)
 rows=[split(line,'\t') for line in readlines(path)[3:end]]
 out=zeros(ComplexF64,maximum(parse(Int,r[1]) for r in rows),maximum(parse(Int,r[2]) for r in rows))
 for r in rows
  out[parse(Int,r[1]),parse(Int,r[2])]=complex(parse(Float64,r[5]),parse(Float64,r[6]))
 end
 out
end
function compare(a,b,label)
 @assert size(a)==size(b)
 for part in (real,imag)
  x,y=part.(a),part.(b);scale=maximum(abs,y);err=maximum(abs,x-y)
  @assert all(abs.(x-y).<=2e-9.*abs.(y).+100eps(Float64)*scale) "$label $part error=$err scale=$scale"
  println(label," ",part," maxabs=",err," scale=",scale)
 end
 flush(stdout)
end
function wiresystem(name,rho;insulated=false,connections=[1,2],poses=[(0.,.1),(.2,-.1)])
 copper=Material(kind=:conductor,rho=rho)
 body=insulated ? terminal(:core,core(copper;r=.005),insulation(Material(kind=:insulator,rho=Inf,eps_r=2.3);t=.002)) : terminal(:core,core(copper;r=.005))
 wire=build(CableDesign,name,body)
 build(LineCableSystem,fill(wire,length(poses)),poses;connections=[Dict(:core=>i) for i in connections])
end
problem(system)=LineParametersProblem(system;frequencies=[50.],earth_props=homogeneous(rho=100.,eps_r=10.))
function native_run(bundle;mesh=nothing,basis=0,meshonly=false)
 entry=joinpath(bundle,"study.pro")
 if mesh===nothing || meshonly
  gmsh=Gmsh.gmsh_jll.gmsh()
  run(Cmd(`$gmsh $(joinpath(bundle,"study.geo")) -setnumber BuildMesh 1 -0 -v 2`;dir="/tmp"))
  mesh=joinpath(bundle,"study.msh")
 end
 meshonly && return mesh
 run(Cmd(`$getdp $entry -msh $mesh -solve LineCableModelsFEM -setnumber BasisTerminal $basis -v 2`;dir="/tmp"))
 out=joinpath(bundle,"results","f0001-quasi-fw-b"*lpad(basis,4,'0'))
 @assert isfile(joinpath(out,"completed.txt"))
 out
end

function exportref(label,p,form)
 ref=compute(p,form;options=settings)
 folder=joinpath(root,label)
 export_data(:onelab,p,form;file_name=joinpath(folder,"study.pro"),mesh_options=(pml_layers=8,))
 mesh=joinpath(details(ref).data.fem.run.run_directory,"mesh","model.msh")
 out=native_run(folder;mesh)
 for (q,expected) in (("Z",Z(ref)[:,:,1]),("Y",Y(ref)[:,:,1]),("P-primitive",details(ref).data.fem.primitive.P_primitive[:,:,1]))
  compare(matrixfile(joinpath(out,"matrices",q*".tsv")),expected,label*"/"*q)
 end
 folder,ref,mesh,out
end
# Grounded terminal plus nontrivial phase ordering, then bundled insulated cables.
for (label,insulated,connections,bundle) in (("grounded",false,[2,0,1],false),("composite-bundle",true,[1,0,1],true))
 p=problem(wiresystem(label,1.72e-8;insulated,connections,poses=[(0.,.1),(.2,-.1),(.4,.1)]))
 form=Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,reduce_bundle=bundle,kron_reduction=true,ideal_transposition=false))
 exportref(label,p,form)
end
p=problem(wiresystem("edited",1.72e-8))
form=Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
folder,ref,mesh,base=exportref("edits",p,form)
cp(base,joinpath(root,"unedited-reference"));base=joinpath(root,"unedited-reference")
datafile=joinpath(folder,"study_data.pro");original=read(datafile,String)
# Edit the native named coefficient; the same mesh isolates the constitutive change.
edited=replace(original,r"Material_1_conductor_Sigma\(\) = \{[^}]*\};"=>"Material_1_conductor_Sigma() = {"*@sprintf("%.17g",inv(2*1.72e-8))*"};")
@assert edited!=original
write(datafile,edited);hash=bytes2hex(open(sha256,datafile))
changed=compute(problem(wiresystem("edited",2*1.72e-8)),form;options=(settings...,mesh_path=mesh))
out=native_run(folder;mesh)
@assert bytes2hex(open(sha256,datafile))==hash
compare(matrixfile(joinpath(out,"matrices","Z.tsv")),Z(changed)[:,:,1],"material-edit/Z")
compare(matrixfile(joinpath(out,"matrices","Y.tsv")),Y(changed)[:,:,1],"material-edit/Y")
# Double the actual current amplitude and retain visible coefficient normalization.
write(datafile,replace(original,"UnitSource = 1.;"=>"UnitSource = 2.;"))
hash=bytes2hex(open(sha256,datafile));out=native_run(folder;mesh,basis=1)
@assert bytes2hex(open(sha256,datafile))==hash
compare(matrixfile(joinpath(out,"matrices","Z-primitive.tsv")),details(ref).data.fem.primitive.Z_primitive[:,1:1,1],"amplitude-edit/normalized-Z")
for q in ("az","b","e")
 a=import_data(:pos,joinpath(out,"maps",q*"_f0001_b0001.pos"))
 b=import_data(:pos,joinpath(base,"maps",q*"_f0001_b0001.pos"))
 @assert length(a.blocks)==length(b.blocks)
 for (i,(x,y)) in enumerate(zip(a.blocks,b.blocks))
  compare(x.values,2 .* y.values,"amplitude-edit/"*q*"/$i")
 end
end
# Remesh via an editable native setting, then compare both pipelines on that mesh.
write(datafile,replace(original,"MeshScale = {1.,"=>"MeshScale = {1.5,"))
newmesh=native_run(folder;meshonly=true)
@assert bytes2hex(open(sha256,newmesh))!=bytes2hex(open(sha256,mesh))
remeshed=compute(p,form;options=(settings...,mesh_path=newmesh))
out=native_run(folder;mesh=newmesh)
compare(matrixfile(joinpath(out,"matrices","Z.tsv")),Z(remeshed)[:,:,1],"remesh/Z")
compare(matrixfile(joinpath(out,"matrices","Y.tsv")),Y(remeshed)[:,:,1],"remesh/Y")
@assert !isfile(joinpath(out,"paths.pro"))
# Total-length output changes final Z/Y units, preserving primitive quantities.
perlength=read(datafile,String)
cp(out,joinpath(root,"per-metre-reference"));out=joinpath(root,"per-metre-reference")
try
 write(datafile,replace(perlength,"LineLength = 1;"=>"LineLength = 7;",
     "OutputTotal = {0,"=>"OutputTotal = {1,"))
 total=native_run(folder;mesh=newmesh)
 compare(matrixfile(joinpath(total,"matrices","Z.tsv")),7 .* Z(remeshed)[:,:,1],"total/Z")
 compare(matrixfile(joinpath(total,"matrices","Y.tsv")),7 .* Y(remeshed)[:,:,1],"total/Y")
 for (quantity,unit) in (("Z","ohm"),("Y","S"),("P","ohm m"))
  @assert startswith(first(readlines(joinpath(total,"matrices",quantity*".tsv"))),"# "*unit*";")
 end
 for quantity in ("Z-primitive","P-primitive","P")
  compare(matrixfile(joinpath(total,"matrices",quantity*".tsv")),
      matrixfile(joinpath(out,"matrices",quantity*".tsv")),"total/"*quantity)
 end
finally
 write(datafile,perlength)
end
# Native constraint and connection edits remain authoritative too.
source=joinpath(folder,"formulations","quasi-full.pro");equations=read(source,String)
try
 write(source,replace(equations,"Value \$FEM_I~{t};"=>"Value -\$FEM_I~{t};"))
 local out=native_run(folder;mesh=newmesh)
 compare(matrixfile(joinpath(out,"matrices","Z.tsv")),-Z(remeshed)[:,:,1],"constraint-edit/Z")
 compare(matrixfile(joinpath(out,"matrices","Y.tsv")),-Y(remeshed)[:,:,1],"constraint-edit/Y")
 @assert read(source,String)!=equations
finally
 write(source,equations)
end
try
 write(datafile,replace(perlength,"Connection_2 = 2;"=>"Connection_2 = 0;",
     "KronReduction = 0;"=>"KronReduction = 1;"))
 local out=native_run(folder;mesh=newmesh)
 primitive=details(remeshed).data.fem.primitive
 options=Formulation(:LineCableModelsFEM;options=(reduce_bundle=false,kron_reduction=true,ideal_transposition=false)).options
 reduced=LineCableModels.Engine.reduce_primitive_matrices(primitive.Z_primitive,primitive.P_primitive,[1,0],options)
 compare(matrixfile(joinpath(out,"matrices","Z.tsv")),reduced.Z[:,:,1],"connection-edit/Z")
 compare(matrixfile(joinpath(out,"matrices","Y.tsv")),inv(reduced.P[:,:,1]),"connection-edit/Y")
finally
 write(datafile,perlength)
end
println("ALL EXTRA VALIDATIONS PASSED: ",root)
