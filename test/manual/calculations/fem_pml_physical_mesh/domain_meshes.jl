# Mesh-only preparation of the pending two-skin-depth probes. No GetDP call.
# Preserve the full frequency vector so conductor grading remains identical.
using LineCableModels, Gmsh, JSON3, TOML, SHA, Dates

const ROOT = joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-physical-mesh")
const FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
const problem = LineCableModels.ImportExport.deserialize_value(JSON3.read(read(
    joinpath(dirname(ROOT),"runs/run-WFpqS5/input/problem.json"),String),Dict{String,Any}))
const form = Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,Γ=0,
    reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
const controls = (domain_skin_depths=2.,
    pml_resolution=(interpolation_cells=72,coefficient_change=.12),
    mesh_size_factor=3.,exterior_mesh_size_factor=8.,volume_quadrature=12,
    conductor_geometry_tolerance=1e-3,conductor_skin_depth_elements=6.,
    conductor_mesh_growth=sqrt(1.25),conductor_skin_depths=5.,conductor_thickness_elements=4)
say(xs...) = (println(Dates.now()," ",xs...);flush(stdout))

function counts(path)
    session = FEM._start_gmsh(0)
    try
        gmsh.open(path)
        nodes = length(first(gmsh.model.mesh.get_nodes()))
        elements = sum(length,gmsh.model.mesh.get_elements(2)[2];init=0)
        groups = Dict{String,Int}()
        for (dimension,tag) in gmsh.model.get_physical_groups(2)
            name = gmsh.model.get_physical_name(dimension,tag)
            startswith(name,"LCM/domain/") || continue
            groups[name] = sum(gmsh.model.get_entities_for_physical_group(2,tag);init=0) do entity
                sum(length,gmsh.model.mesh.get_elements(2,entity)[2];init=0)
            end
        end
        return (;nodes,elements,groups)
    finally
        FEM._finish_gmsh(session)
    end
end

mkpath(joinpath(ROOT,"domain-probes"))
open(joinpath(ROOT,"domain-probes/mesh-costs.csv"),"w") do csv
    println(csv,"case,frequency_hz,domain_halfwidth_m,pml_thickness_m,nodes,triangles,pml_triangles,mesh_seconds")
    for (extent,factor) in (("l24",12.),("l2",1.)), (name,index) in (("low",1),("mid",4))
        label = "d2-$extent-$name"
        directory = joinpath(ROOT,"domain-probes",label)
        entry = joinpath(directory,"study.pro")
        mesh = joinpath(directory,"study.msh")
        marker = joinpath(directory,"mesh.toml")
        options = (;controls...,pml_thickness_factor=factor)
        model = FEM._resolved_fem_model(problem,form,
            computation_options(LineCableModelsFEM,ComputationOptions(;options...)))
        plan = model.mesh_plans[index]
        if !isfile(marker)
            say("DOMAIN MESH BEGIN ",label," f=",plan.frequency,"; no field solve")
            isfile(entry) || export_data(:onelab,problem,form;file_name=entry,mesh_options=options)
            mesher = Gmsh.gmsh_jll.gmsh()
            command = `$mesher $(joinpath(directory,"study.geo")) -setnumber FrequencyIndex $index -setnumber BuildMesh 1 -0 -v 4`
            streaming = Cmd(Cmd(["stdbuf","-oL","-eL",command.exec...]);env=command.env)
            seconds = @elapsed run(pipeline(streaming;stdout,stderr=stdout))
            metrics = counts(mesh)
            record = Dict("frequency_index"=>index,"frequency_hz"=>plan.frequency,
                "domain_halfwidth_m"=>plan.domain_halfwidth,"pml_thickness_m"=>collect(plan.pml_thickness),
                "pml_strength"=>collect(plan.pml_strength),"normal_intervals"=>collect(plan.pml_layers),
                "domain_mesh_size_m"=>plan.domain_mesh_size,
                "exterior_mesh_sizes_m"=>collect(plan.exterior_mesh_sizes),
                "mesh_seconds"=>seconds,"nodes"=>metrics.nodes,"triangles"=>metrics.elements,
                "groups"=>metrics.groups,"mesh_sha256"=>bytes2hex(open(sha256,mesh)),
                "field_solve_performed"=>false)
            open(io -> TOML.print(io,record),marker*".tmp","w")
            mv(marker*".tmp",marker;force=true)
        end
        result = TOML.parsefile(marker)
        result["mesh_sha256"] == bytes2hex(open(sha256,mesh)) || error("Prepared mesh changed: $label")
        pml = result["groups"]["LCM/domain/pml"]
        println(csv,join((label,result["frequency_hz"],result["domain_halfwidth_m"],
            first(result["pml_thickness_m"]),result["nodes"],result["triangles"],pml,result["mesh_seconds"]),','))
        flush(csv)
        say("DOMAIN MESH DONE ",label," nodes=",result["nodes"]," triangles=",result["triangles"],
            " PML triangles=",pml," intervals=",result["normal_intervals"])
    end
    for (label,name) in (("d24-l24-low","diagnose-strict-low"),("d24-l24-mid","diagnose-strict-mid"))
        marker = TOML.parsefile(joinpath(ROOT,name,"complete.toml"))
        metrics = counts(joinpath(ROOT,name,"detached/study.msh"))
        delta = sqrt(.1/(pi*marker["frequency"]*4pi*1e-7))
        println(csv,join((label,marker["frequency"],24delta,24delta,metrics.nodes,
            metrics.elements,metrics.groups["LCM/domain/pml"],marker["mesh_seconds"]),','))
    end
end
say("DOMAIN MESH PREPARATION COMPLETE; four meshes prepared; zero additional field solves")
