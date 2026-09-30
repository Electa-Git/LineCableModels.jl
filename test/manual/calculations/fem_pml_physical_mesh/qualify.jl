# Retained scientific qualification; public compute only, no prototype override.
using LineCableModels, Gmsh, JSON3, TOML, LinearAlgebra, Printf, Dates
const ROOT=joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-physical-mesh")
const REPO=pkgdir(LineCableModels)
include(joinpath(REPO,"test/support/scenarios.jl"))
const REDUCTIONS=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
const FEM=Formulation(:LineCableModelsFEM;options=(;REDUCTIONS...,physics=:quasi_fw))
const ANALYTICAL=Formulation(earth_impedance=formula(:unified;options=(Γ=0.,)),
    earth_admittance=formula(:unified;options=(Γ=0.,));options=REDUCTIONS)
const MESH=(domain_skin_depths=24.,pml_resolution=(interpolation_cells=72,coefficient_change=.12),volume_quadrature=12,
    mesh_size_factor=3.,exterior_mesh_size_factor=8.,conductor_geometry_tolerance=1e-3,
    conductor_skin_depth_elements=6.,conductor_mesh_growth=sqrt(1.25),
    conductor_skin_depths=5.,conductor_thickness_elements=4)
const OPTIONS=(;MESH...,mesh_policy=:remesh,resume_run_directory=:latest,
    keep_run_directory=true,trace=true,timing=true,output_basis=:pul,
    frequency_workers=4,solver_threads=1,plot_field_maps=false,gmsh_verbosity=4,
    getdp_verbosity=4,verbosity=(default=1,),getdp_executable=joinpath(ROOT,"getdp-live"))

say(xs...)=(println(Dates.now()," ",xs...);flush(stdout))
function snapshot(name)
    path=joinpath(REPO,".linecablemodels/fem/runs",name,"input/problem.json")
    LineCableModels.ImportExport.deserialize_value(JSON3.read(read(path,String),Dict{String,Any}))
end
function csv_matrices(path,Z,Y,frequencies)
    open(path,"w") do io
        println(io,"quantity,frequency_hz,receiver,source,real,imaginary")
        for (label,A) in (("Z",Z),("Y",Y)),k in eachindex(frequencies),j in axes(A,2),i in axes(A,1)
            @printf(io,"%s,%.17g,%d,%d,%.17g,%.17g\n",label,frequencies[k],i,j,real(A[i,j,k]),imag(A[i,j,k]))
        end
    end
end
function saved_reference(path,frequencies,n)
    values=Dict(q=>fill(complex(NaN,NaN),n,n,length(frequencies)) for q in ("Z","Y"))
    for line in readlines(path)[2:end]
        row=split(line,',');f=parse(Float64,row[2]);k=findfirst(==(f),frequencies)
        k===nothing && continue
        values[row[1]][parse(Int,row[3]),parse(Int,row[4]),k]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
    end
    all(isfinite,values["Z"])&&all(isfinite,values["Y"]) || error("Incomplete saved reference: $path")
    return values["Z"],values["Y"]
end
function study(label,problem;reference_csv=nothing)
    directory=joinpath(ROOT,"qualification",label)
    mkpath(directory)
    marker=joinpath(directory,"complete.toml")
    if isfile(marker)
        status=TOML.parsefile(marker)
        status["G_sign_mismatches"]==0 || error("Recorded sign regression in $label; inspect before continuing")
        say("REUSE QUALIFICATION ",label)
        return
    end
    n=length(problem.system.terminal_order)
    frequencies=problem.frequencies
    zref,yref=if reference_csv===nothing
        reference=compute(problem,ANALYTICAL)
        Z(reference),Y(reference)
    else
        saved_reference(reference_csv,frequencies,n)
    end
    csv_matrices(joinpath(directory,"reference.csv"),zref,yref,frequencies)
    write(joinpath(directory,"options.txt"),repr(OPTIONS)*"\n")
    say("QUALIFY BEGIN ",label," frequencies=",frequencies," terminals=",n,"; native stdout is streamed to this same log")
    measured=@timed compute(problem,FEM;options=OPTIONS)
    result=measured.value
    runinfo=details(result).data.fem.run
    z,y=Z(result),Y(result)
    csv_matrices(joinpath(directory,"matrices.csv"),z,y,frequencies)
    for q in ("Z","P")
        cp(joinpath(runinfo.run_directory,"raw",q*".tsv"),joinpath(directory,q*"-primitive.tsv");force=true)
    end
    mismatches=0
    open(joinpath(directory,"components.csv"),"w") do io
        println(io,"quantity,frequency_hz,i,j,value,reference,signed_error,absolute_error,relative_error,sign_match")
        for (q,a,b) in (("R",real.(z),real.(zref)),("X",imag.(z),imag.(zref)),
                ("G",real.(y),real.(yref)),("B",imag.(y),imag.(yref)))
            for k in eachindex(frequencies),j in 1:n,i in 1:n
                v,r=a[i,j,k],b[i,j,k]
                match=sign(v)==sign(r)
                # A zero reference does not assert the sign of an unresolved zero;
                # preserve its raw absolute error. No data are clipped.
                if q=="G" && !iszero(r) && !match
                    mismatches+=1
                    say("G SIGN MISMATCH ",label," f=",frequencies[k]," [",i,",",j,"] FEM=",v," reference=",r)
                end
                println(io,join((q,frequencies[k],i,j,v,r,v-r,abs(v-r),iszero(r) ? NaN : abs((v-r)/r),match),','))
            end
        end
    end
    status=Dict("label"=>label,"wall_seconds"=>measured.time,"compile_seconds"=>measured.compile_time,
        "recompile_seconds"=>measured.recompile_time,"gc_seconds"=>measured.gctime,
        "run_directory"=>runinfo.run_directory,"reused"=>runinfo.reused,
        "frequencies"=>length(frequencies),"columns"=>n*length(frequencies),
        "G_sign_mismatches"=>mismatches,"reference"=>reference_csv===nothing ? "matching unified Gamma=0" : reference_csv)
    open(io->TOML.print(io,status),marker*".tmp","w");mv(marker*".tmp",marker;force=true)
    say("QUALIFY DONE ",label," wall_seconds=",measured.time," compile_seconds=",measured.compile_time," G_sign_mismatches=",mismatches)
    mismatches==0 || error("Qualification stopped at a demonstrated conductance sign regression: $label. Completed data are preserved.")
end
say("QUALIFICATION: public physical PML prescription ",MESH.pml_resolution)
# Start with the actual failing radius, then the remaining ordinary study.
for (label,name) in (("radius-0.085","run-WFpqS5"),("radius-0.001","run-KgMjbL"),
        ("radius-0.01","run-YWmojw"),("rho-0.1","run-ZYYfqP"),
        ("rho-1","run-SCTzPh"),("rho-100","run-rDJaLr"),("rho-1000","run-WAGP7m"))
    study(label,snapshot(name))
end
say("ORDINARY STUDY COMPLETE: 70 frequency solves / 140 columns; proceeding to cross-layout preservation")
for layout in (:all_earth,:air_1)
    problem=CurrentScenarios.three_bare_wires_problem(;heights=getproperty(CurrentScenarios.three_bare_wires_layouts,layout),
        frequencies=[.1,21.544346900318832,1e6],name="three_bare_wires_$layout")
    study("three-$layout",problem)
end
say("THREE-WIRE PRESERVATION COMPLETE: six frequency solves / eighteen columns")
# Use the exact saved physical designs from completed conductor qualification.
for (name,reference_folder) in (("screen","interface-delivery"),("tube","cable-fixtures"),("sector","interface-delivery"))
    reference=joinpath(REPO,".linecablemodels/fem/conductor-mesh-qualification",reference_folder,name,"normal/quasi_fw")
    oldrun=basename(TOML.parsefile(joinpath(reference,"cost.toml"))["run_directory"])
    original=snapshot(oldrun)
    problem=LineParametersProblem(original.system;frequencies=[1e6],temperature=original.temperature,earth_props=original.earth_props)
    study("cable-$name",problem;reference_csv=joinpath(reference,"matrices.csv"))
end
say("MANAGED QUALIFICATION COMPLETE: 79 frequency solves / 167 columns; detached preservation and final assessment remain")
