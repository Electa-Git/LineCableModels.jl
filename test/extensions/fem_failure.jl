@testitem "Gmsh FEM / public solver failure retains log and run evidence" tags=[:extension] begin
    using Gmsh
    using JSON3

    if Sys.isunix()
        extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
        copper = Material(kind = :conductor, rho = 1.72e-8)
        design = build(CableDesign, "fem-public-failure",
            terminal(:core, core(copper; r = 0.005)))
        system = build(LineCableSystem, design, (0.0, -0.1);
            connections = Dict(:core=>1), system_id = "fem-public-failure")
        problem = LineParametersProblem(system; frequencies = [50.0],
            earth_props = homogeneous(rho = 100.0, eps_r = 10.0))
        @test !Bool(Gmsh.gmsh.is_initialized())
        mktempdir() do directory
            executable = joinpath(directory, "failing-getdp")
            write(executable, raw"""
                #!/bin/sh
                if [ "$1" = "-info" ]; then
                    echo 'GetDP Version 3.5.0 failure fixture'
                    exit 0
                fi
                echo 'deliberate solver failure' >&2
                exit 17
                """)
            chmod(executable, 0o700)
            log_file = joinpath(directory, "logs", "computation.log")
            formulation = Formulation(:LineCableModelsFEM;
                options = (ideal_transposition = false,),
                fem_options = (getdp_executable = executable, gmsh_verbosity = 0,
                    getdp_verbosity = 0, keep_run_directory = true))
            failure = try
                compute(problem, formulation;
                    options = (log_file = log_file, verbosity = (default = 0,)))
                nothing
            catch exception
                exception
            end
            @test failure isa LineCableModelsFEMError
            if failure isa LineCableModelsFEMError && failure.run_directory !== nothing
                run_path = failure.run_directory
                try
                    @test failure.category === :getdp
                    @test failure.field === :client
                    @test occursin(run_path, sprint(showerror, failure))
                    state = JSON3.read(read(joinpath(run_path, "run.json"), String))
                    @test state.state == "failed"
                    @test state.getdp_invocations == 1
                    @test occursin("GetDP", state.message)
                    @test isfile(joinpath(run_path, "input", "problem.json"))
                    @test isfile(joinpath(run_path, "input", "getdp", "quasi-tem.pro"))
                    @test isfile(joinpath(run_path, "logs", "getdp.log"))
                    @test !isfile(joinpath(run_path, "raw", "checksums.json"))
                    messages = read(log_file, String)
                    @test occursin("FEM geometry ready", messages)
                    @test occursin("FEM mesh ready", messages)
                    @test occursin("Starting isolated GetDP", messages)
                    @test !occursin("FEM scan completed successfully", messages)
                finally
                    # Remove only this fixture's precisely identified failed run.
                    @assert realpath(dirname(run_path)) ==
                            realpath(joinpath(extension._runtime_root(), "runs"))
                    rm(run_path; recursive = true)
                end
            end
            @test !Bool(Gmsh.gmsh.is_initialized())
        end
    end
end
