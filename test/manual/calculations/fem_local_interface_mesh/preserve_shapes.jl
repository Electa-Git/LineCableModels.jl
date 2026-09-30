# Mesh-only checks of the already qualified screen, tubular and sector fixtures.
isdefined(@__MODULE__,:ROOT) || include("qualify.jl")
include("../fem_conductor_mesh/fixtures.jl")

function shape_mesh(dir, problem, candidate; media=:wave)
    session = FEM._start_gmsh(2)
    try
        gmsh.parser.clear(); gmsh.onelab.clear()
        gmsh.parser.set_string("OnelabAction",["audit"])
        gmsh.open(joinpath(dir,"study.geo"))
        candidate && localize!(problem;media=media==:decay ? :both : media,decay_footprint=media==:decay)
        gmsh.model.mesh.generate(2)
        n = length(problem.system.terminal_order)
        group(d,t) = [(d,c) for c in gmsh.model.get_entities_for_physical_group(d,t)]
        # Conductor sizing is prescribed; its unstructured interior need not
        # reproduce identical node coordinates on an independent Gmsh run.
        return Dict("pml"=>coordinate_hash(group(2,1005)),
            "contours"=>coordinate_hash(vcat([group(1,4000+i) for i in 1:n]...)),
            "paths"=>coordinate_hash(vcat([group(1,7000+i) for i in 1:n]...)),
            "metal"=>coordinate_hash(group(2,6001)),
            "metal_triangles"=>sum(length(block) for (_,s) in group(2,6001)
                for block in gmsh.model.mesh.get_elements(2,s)[2]),
            "conductor_controls"=>Dict(name=>gmsh.parser.get_number(name)
                for name in gmsh.parser.get_names() if startswith(name,"ConductorRegion")),
            "nodes"=>length(first(gmsh.model.mesh.get_nodes())))
    finally
        FEM._finish_gmsh(session)
    end
end

function preserve_shapes(; media=:wave)
    for (name,design,height) in (("screen",ConductorMeshFixtures.screened_cable(),-1.),
            ("tube",ConductorMeshFixtures.tubular_cable(),-1.),
            ("sector",ConductorMeshFixtures.sector_cable(),1.))
        dir = mkpath(joinpath(ROOT,"shapes",name))
        suffix = media==:wave ? "" : "-$media"
        marker = joinpath(dir,"preservation$suffix.toml")
        isfile(marker) && (say("REUSE SHAPE ",name); continue)
        system = build(LineCableSystem,design,Pose2(0.,height);
            connections=Dict(t=>i for (i,t) in enumerate(design.terminal_order)))
        problem = LineParametersProblem(system;frequencies=[1e6],temperature=20.,
            earth_props=homogeneous(rho=100.,eps_r=1.,mu_r=1.))
        form = Formulation(:LineCableModelsFEM;options=REDUCTIONS)
        source = media==:production ? joinpath(dir,"production") : dir
        prepare_bundle(source,problem,form)
        saved = joinpath(dir,"mesh-comparison.toml")
        a = isfile(saved) ? TOML.parsefile(saved)["baseline"] : shape_mesh(dir,problem,false)
        b = shape_mesh(source,problem,media!=:production;media)
        record(joinpath(dir,"mesh-comparison$suffix.toml"),Dict("baseline"=>a,"localized"=>b,
            "metal_mesh_identical"=>a["metal"]==b["metal"]))
        isempty(a["conductor_controls"]) && error("Missing conductor prescription")
        for key in ("pml","contours","paths","conductor_controls")
            a[key] == b[key] || error("Shape preservation failure $name: $key")
        end
        record(marker,Dict("baseline"=>a,"localized"=>b,"passed"=>true))
        say("SHAPE PRESERVATION PASS ",name," nodes ",a["nodes"]," -> ",b["nodes"],
            " metal triangles ",a["metal_triangles"]," -> ",b["metal_triangles"],
            " identical metal coordinates=",a["metal"]==b["metal"])
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    preserve_shapes(;media=isempty(ARGS) ? :wave : Symbol(only(ARGS)))
end
