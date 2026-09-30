# Reproduce every mesh of the interactive runner without starting its plots or
# solver grid. This specifically retains frequency order and conductor CAD reuse.
module RunnerSetup
    using LineCableModels
    runner = joinpath(pkgdir(LineCableModels),"test/manual/calculations/run_two_bare_wires_fem.jl")
    setup = first(split(read(runner,String),"# One detached case:"))
    include_string(@__MODULE__,replace(setup,"using Revise\n"=>"","using GLMakie\n"=>""),runner)
end

using LineCableModels, Gmsh

function verify_runner_meshes()
    fem = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    output = mkpath(joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-corner-mesh"))
    cases = [(radius=RunnerSetup.fixed_external_radius,rho=rho)
        for rho in RunnerSetup.soil_resistivity_grid]
    append!(cases,[(radius=r,rho=RunnerSetup.fixed_soil_resistivity)
        for r in RunnerSetup.external_radius_grid])
    execution = computation_options(LineCableModelsFEM,ComputationOptions(;RunnerSetup.fem_options...))
    open(joinpath(output,"runner-mesh-check.csv"),"w") do table
        println(table,"radius_m,rho_ohm_m,frequency_hz,nodes,triangles,seconds")
        for case in cases
            problem = RunnerSetup.build_two_bare_wires_problem(case.radius,case.rho,RunnerSetup.frequency_grid)
            model = fem._resolved_fem_model(problem,RunnerSetup.fem_formulation,execution)
            session = fem._start_gmsh(2)
            try
                geometry = fem._build_geometry!(model,"runner-mesh-check")
                for plan in model.mesh_plans
                    elapsed = @elapsed begin
                        fem._update_exterior_mesh!(model,geometry,plan)
                        fem._configure_mesh!(model,geometry,plan)
                        gmsh.model.mesh.generate(2)
                        fem._inspect_loaded_mesh(model,"runner-$(case)-$(plan.frequency_index)")
                    end
                    nodes = length(first(gmsh.model.mesh.get_nodes()))
                    triangles = sum(length,gmsh.model.mesh.get_elements(2)[2];init=0)
                    println(table,join((case.radius,case.rho,plan.frequency,nodes,triangles,elapsed),','))
                    flush(table)
                    println("RUNNER MESH DONE radius=",case.radius," rho=",case.rho,
                        " f=",plan.frequency," nodes=",nodes," triangles=",triangles)
                    flush(stdout)
                end
            finally
                fem._finish_gmsh(session)
            end
        end
    end
    println("COMPLETE all manual-runner meshes; no solver grid or plotting was started")
    flush(stdout)
end

if abspath(PROGRAM_FILE)==@__FILE__
    output = mkpath(joinpath(pkgdir(LineCableModels),".linecablemodels/fem/pml-corner-mesh"))
    open(joinpath(output,"live.log"),"a") do log
        redirect_stdio(stdout=log,stderr=log) do
            verify_runner_meshes()
        end
    end
end
