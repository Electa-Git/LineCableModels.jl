# Bounded native mesh/domain study using the existing model and plan types.
# It changes test-local plans only; no production meshing policy is replaced.
module FEMMeshControl
using Gmsh,Printf,TOML
using LineCableModels: computation_options,LineCableModelsFEM

rebuild(value;kwargs...)=typeof(value)((get(kwargs,name,getfield(value,name))
    for name in fieldnames(typeof(value)))...)

function study(FEM,problem,formulation,controls;directory,observable)
    model=FEM._resolved_fem_model(FEM._preflight_fem_problem(problem),formulation)
    plan=only(model.mesh_plans)
    # Start one mesh scale coarser, then complete all three refinements. The
    # former 500000-triangle allocation could not supply four mesh levels and
    # three domains. Only this test-owned construction allocation changes.
    plan=rebuild(plan;domain_mesh_size=2plan.domain_mesh_size,
        infinite_mesh_size=2plan.infinite_mesh_size,
        interface_mesh_size=2plan.interface_mesh_size,
        cable_interface_mesh_sizes=2plan.cable_interface_mesh_sizes,
        wave_mesh_sizes=map(x->2x,plan.wave_mesh_sizes))
    model=rebuild(model;region_plans=[rebuild(region;mesh_size=2region.mesh_size)
            for region in model.region_plans],
        fine_mesh_size=2model.fine_mesh_size,coarse_mesh_size=2model.coarse_mesh_size,
        cable_outer_mesh_sizes=2model.cable_outer_mesh_sizes,mesh_plans=[plan])
    execution=computation_options(LineCableModelsFEM,controls)
    executable=FEM._getdp_selection(execution).path
    records=Dict{String,Any}[];mesh_values=Any[];domain_values=Any[]
    affordable=3
    for factor in (1,2,4)
        changed_plan=rebuild(plan;domain_radius=factor*plan.domain_radius,
            shell_outer_radius=factor*plan.shell_outer_radius,
            domain_mesh_size=factor*plan.domain_mesh_size,
            infinite_mesh_size=factor*plan.infinite_mesh_size)
        changed_model=rebuild(model;domain_radius=changed_plan.domain_radius,
            shell_outer_radius=changed_plan.shell_outer_radius,mesh_plans=[changed_plan])
        run=FEM._create_run(directory)
        session=FEM._start_gmsh(0)
        try
            geometry=FEM._build_geometry!(changed_model,"refinement-$factor",changed_plan)
            base_mesh=only(FEM._select_meshes!(run,changed_model,geometry,execution,directory))
            FEM._prepare_run_inputs!(run,changed_model)
            # Retain the CAD model while refining so curved boundaries remain
            # tied to their present geometric definition.
            for level in 0:affordable
                kinds,tags,_=gmsh.model.mesh.get_elements(2)
                all(==(2),kinds) || error("inconclusive mesh study: non-first-order triangles")
                triangles=sum(length,tags)
                triangles<=2000000 || error("inconclusive mesh study: 2000000-triangle bound")
                mesh=joinpath(run.path,"mesh","level$level.msh")
                gmsh.write(mesh)
                FEM._inspect_loaded_mesh(changed_model,mesh)
                if factor==1 || level==affordable
                    work=joinpath(run.path,"study-level$level");mkpath(work)
                    if formulation.options.physics===Symbol("quasi-fw")
                        path=FEM._voltage_path_file(run,changed_plan.frequency_index)
                        FEM._write_voltage_paths(path,mesh,changed_plan,changed_model;shell_segments=512)
                    end
                    source=joinpath(run.path,"input/getdp/model.pro")
                    command=FEM._getdp_command(executable,source,mesh,run,formulation,
                        execution,changed_plan,[1,2],work;reuse_factorization=false)
                    open(joinpath(work,"getdp.log"),"w") do io
                        Base.run(pipeline(command;stdout=io,stderr=io))
                    end
                    matrices=map(("Z","P")) do quantity
                        matrix=zeros(ComplexF64,2,2)
                        for basis in 1:2
                            raw=joinpath(work,"raw/jobs",@sprintf("getdp-f0001-b%04d-%s.tsv",basis,quantity))
                            FEM._valid_job_raw(raw,2,1,plan.frequency,basis) || error("invalid study output")
                            for line in eachline(raw)
                                row=split(line)
                                matrix[parse(Int,row[3]),basis]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
                            end
                        end
                        matrix
                    end
                    value=observable(matrices...)
                    factor==1 && push!(mesh_values,value)
                    level==affordable && push!(domain_values,value)
                    push!(records,Dict("factor"=>factor,"level"=>level,"triangles"=>triangles,
                        "mesh"=>mesh,"run"=>work,"real"=>real.(vec(value)),"imag"=>imag.(vec(value))))
                    open(joinpath(directory,"mesh-domain.toml"),"w") do io
                        TOML.print(io,Dict("cases"=>records,"scope"=>"native discretization and domain convergence"))
                    end
                end
                level==affordable && break
                4triangles<=2000000 || error("inconclusive mesh study: cannot complete the allocated four levels")
                gmsh.model.mesh.refine()
            end
        finally
            FEM._finish_gmsh(session)
        end
    end
    return (;mesh_values,domain_values,records)
end
end
