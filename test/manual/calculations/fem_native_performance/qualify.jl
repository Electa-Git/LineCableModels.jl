# Manual qualification, never imported by production. Serial and resumable.
module Reference
include("../fem_local_interface_mesh/qualify.jl")
end
using LineCableModels, Gmsh, SHA, TOML, Dates, Printf, Statistics
const B = Reference
const FEM = B.FEM
const gmsh = Gmsh.gmsh
const ROOT = joinpath(pkgdir(LineCableModels), ".linecablemodels/fem/native-performance-20260930")
const PILOTS = ((:air,.1,0.), (:mixed,.1,.99), (:mixed,1e6,.99))
const VARIANTS = ("baseline", "harmonic", "physical3", "quad9", "quad16", "binary", "amd")
label(c) = "$(c[1])-f$(c[2])-gamma$(c[3])"
say(xs...) = B.say(xs...)
record(args...) = B.record(args...)
digest(path) = B.digest(path)

function freeze!()
    marker = joinpath(ROOT,"snapshot.toml")
    isfile(marker) && return
    mkpath(ROOT)
    paths = ["src/engine/options.jl", "test/manual/calculations/run_two_bare_wires_fem.jl"]
    base = joinpath(pkgdir(LineCableModels),"ext/LineCableModelsGmshExt")
    for (dir,_,files) in walkdir(base), file in files
        push!(paths,relpath(joinpath(dir,file),pkgdir(LineCableModels)))
    end
    for file in paths
        dst=joinpath(ROOT,"snapshot",file); mkpath(dirname(dst))
        cp(joinpath(pkgdir(LineCableModels),file),dst)
    end
    record(marker,Dict("files"=>Dict(f=>digest(joinpath(ROOT,"snapshot",f)) for f in paths)))
    say("FROZEN ",length(paths)," current source files; original sources and runner untouched")
end

function harmonic!(path)
    s = read(path,String)
    firstterm = findfirst("      Galerkin {\n        DtDof [sigmaZ[]",s)
    ending = findnext("        GlobalTerm { [Dof{I}, {U}];",s,last(firstterm))
    replacement = """
      // Harmonic admittivity combines the identical conductivity/permittivity supports.
      Galerkin { [Complex[0,omega[]]*seZ[]*Dof{a}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [seZ[]*Dof{ur}, {a}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [Complex[0,omega[]]*seZ[]*Dof{a}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1; }
      Galerkin { [seZ[]*Dof{ur}, {ur}];
        In Domain_Mag; Jacobian Vol; Integration I1; }

"""
    write(path,s[1:first(firstterm)-1]*replacement*s[first(ending):end])
end

function integration!(path, variant)
    old = read(path,String)
    old_i1 = first(split(old,"  { Name I2;"))
    pml = startswith(variant,"quad")
    n = pml ? parse(Int,variant[5:end]) : 3
    # Criterion is a native, zero-based rule selector; it changes no weak terms.
    fun = "Function { NativeRule[All] = 0; NativeRule[Region[{AIR_PML, EARTH_PML}]] = 1; }\n"
    physical = pml ? "VolumeQuadrature" : "3"
    second = pml ? "{ Type GaussLegendre; Case { { GeoElement Quadrangle; NumberOfPoints $n; } } }" :
        "{ Type Gauss; Case { { GeoElement Triangle; NumberOfPoints VolumeQuadrature; } } }"
    head = """
Integration {
  { Name I1; Criterion NativeRule[];
    Case {
      { Type Gauss; Case {
        { GeoElement Point; NumberOfPoints 1; }
        { GeoElement Line; NumberOfPoints 4; }
        { GeoElement Triangle; NumberOfPoints $physical; }
        { GeoElement Quadrangle; NumberOfPoints 4; }
        { GeoElement Triangle2; NumberOfPoints 7; }
      } }
      $second
    }
  }
"""
    write(path,fun*replace(old,old_i1=>head;count=1))
end

function entity_hash(entities)
    io=IOBuffer()
    for (d,s) in sort(entities)
        print(io,(d,s),gmsh.model.mesh.get_elements(d,s))
        nt,xyz,_=gmsh.model.mesh.get_nodes(d,s,true,false)
        order=sortperm(nt); print(io,nt[order],reshape(xyz,3,:)[:,order])
    end
    bytes2hex(sha256(take!(io)))
end

function recombine_pml!(pml)
    converted=[]
    for s in pml
        nt,xyz,_=gmsh.model.mesh.get_nodes(2,s,true,false)
        owned,owned_xyz,_=gmsh.model.mesh.get_nodes(2,s,false,false)
        types,tags,nodes=gmsh.model.mesh.get_elements(2,s)
        types==[2] || error("Reference PML must be linear triangles")
        coords=Dict(t=>Tuple(x) for (t,x) in zip(nt,eachcol(reshape(xyz,3,:))))
        # Qualification copies: remove the known rectangular-cell diagonals.
        # Gmsh's general Blossom recombiner is not the transfinite recombiner:
        # it can pair across grid lines on strongly stretched cells. Production
        # qualification must subsequently check native Transfinite+Recombine.
        unmatched=Dict{Tuple{UInt64,UInt64},Vector{UInt64}}()
        qn=UInt64[]
        for tri in eachcol(reshape(only(nodes),3,:))
            # Imported transfinite grids have small CAD interpolation drift.
            # This only identifies grid diagonals; no node is moved or rounded.
            edges=((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
            tol=1e-4minimum(hypot(coords[a][1]-coords[b][1],coords[a][2]-coords[b][2]) for (a,b) in edges)
            diagonals=[(min(a,b),max(a,b)) for (a,b) in
                edges if
                abs(coords[a][1]-coords[b][1])>tol && abs(coords[a][2]-coords[b][2])>tol]
            length(diagonals)==1 || error("Nonrectangular source triangle: surface=$s points=$([coords[t] for t in tri]) tolerance=$tol")
            diagonal=only(diagonals)
            if !haskey(unmatched,diagonal)
                unmatched[diagonal]=collect(tri); continue
            end
            quad=unique([pop!(unmatched,diagonal);tri])
            length(quad)==4 || error("Bad rectangular triangle pair")
            cx=sum(coords[t][1] for t in quad)/4; cy=sum(coords[t][2] for t in quad)/4
            sort!(quad;by=t->atan(coords[t][2]-cy,coords[t][1]-cx))
            points=[coords[t] for t in quad]
            all(k->abs(points[k][1]-points[mod1(k+1,4)][1])<=tol ||
                abs(points[k][2]-points[mod1(k+1,4)][2])<=tol,1:4) ||
                error("Recombination changed rectangular grid cells")
            append!(qn,quad)
        end
        isempty(unmatched) || error("Unpaired PML triangles")
        push!(converted,(s,owned,owned_xyz,qn))
    end
    gmsh.model.mesh.clear([(2,s) for s in pml])
    for (s,nt,xyz,qn) in converted
        gmsh.model.mesh.add_nodes(2,s,nt,xyz)
        gmsh.model.mesh.add_elements_by_type(s,3,Int[],qn)
    end
end

function transform_mesh!(dir, variant)
    session=FEM._start_gmsh(3)
    try
        gmsh.parser.clear(); gmsh.onelab.clear()
        gmsh.parser.set_string("OnelabAction",["audit"])
        gmsh.open(joinpath(dir,"study.msh"))
        pml=Int.(gmsh.model.get_entities_for_physical_group(2,1005))
        physical=[(d,s) for (d,s) in gmsh.model.get_entities(2) if s ∉ pml]
        paths=[(1,s) for tag in (7001,7002) for s in gmsh.model.get_entities_for_physical_group(1,tag)]
        before=entity_hash([physical;paths])
        measured=@timed begin
            if startswith(variant,"quad")
                recombine_pml!(pml)
            end
        end
        before==entity_hash([physical;paths]) || error("Physical mesh/path changed")
        counts=Dict{String,Int}()
        for s in pml
            types,tags,_=gmsh.model.mesh.get_elements(2,s)
            for (t,ids) in zip(types,tags); counts[string(t)]=get(counts,string(t),0)+length(ids); end
        end
        if startswith(variant,"quad")
            Set(keys(counts))==Set(["3"]) || error("PML not entirely quadrangular: $counts")
        end
        gmsh.option.set_number("Mesh.Binary",variant=="binary" ? 1 : 0)
        gmsh.option.set_number("Mesh.MshFileVersion",4.1)
        gmsh.option.set_number("Mesh.SaveAll",1)
        write_time=@elapsed gmsh.write(joinpath(dir,"study.msh"))
        record(joinpath(dir,"mesh.toml"),Dict("mesh_seconds"=>measured.time,
            "compile_seconds"=>measured.compile_time,"write_seconds"=>write_time,
            "pml_elements_by_type"=>counts,"physical_path_hash"=>before,
            "sha256"=>digest(joinpath(dir,"study.msh")),"bytes"=>filesize(joinpath(dir,"study.msh"))))
    finally
        gmsh.option.set_number("Mesh.MeshOnlyEmpty",0); gmsh.option.set_number("Mesh.Binary",0)
        FEM._finish_gmsh(session)
    end
end

function prepare!(case,variant)
    dir=joinpath(ROOT,label(case),variant)
    isfile(joinpath(dir,"sources.toml")) && return dir
    src=joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-corner-mesh",label(case),"baseline")
    if !isdir(src)
        src=joinpath(B.ROOT,label(case),case[3]==0 ? "localized" : "localized-terminal-shift")
    end
    files=readlines(joinpath(src,".onelab-export-files"))
    for f in [files;".onelab-export-files";"study.msh"]
        dst=joinpath(dir,f); mkpath(dirname(dst)); cp(joinpath(src,f),dst;force=true)
    end
    for f in readdir(joinpath(dir,"formulations"))
        frozen=joinpath(ROOT,"snapshot/ext/LineCableModelsGmshExt/getdp",f)
        isfile(frozen) && cp(frozen,joinpath(dir,"formulations",f);force=true)
    end
    variant in ("harmonic","combined") && harmonic!(joinpath(dir,"formulations/quasi-full.pro"))
    (variant in ("physical3","combined") || startswith(variant,"quad")) && integration!(joinpath(dir,"formulations/integration.pro"),variant)
    if startswith(variant,"quad") || variant=="binary"
        transform_mesh!(dir,variant)
    else
        record(joinpath(dir,"mesh.toml"),Dict("mesh_seconds"=>0.,"compile_seconds"=>0.,
            "sha256"=>digest(joinpath(dir,"study.msh")),"bytes"=>filesize(joinpath(dir,"study.msh"))))
    end
    record(joinpath(dir,"sources.toml"),Dict(f=>digest(joinpath(dir,f)) for f in files))
    return dir
end

function pilots!()
    freeze!()
    for case in PILOTS
        problem,form=B.fixture(case...)
        for variant in VARIANTS
            dir=prepare!(case,variant)
            say("BEGIN ",label(case)," / ",variant)
            B.solve!(dir;ordering=variant=="amd" ? 0 : nothing)
            if variant!="baseline" && !isfile(joinpath(dirname(dir),"comparison-$variant.toml"))
                B.compare!(dirname(dir),problem,form;candidate=variant,tolerance=.02)
            end
        end
    end
    say("COMPLETE initial native pilots; no production changes promoted")
end

if abspath(PROGRAM_FILE)==@__FILE__
    pilots!()
end
