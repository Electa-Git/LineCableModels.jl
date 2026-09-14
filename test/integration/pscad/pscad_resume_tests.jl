@testitem "PSCAD / verified completed-run reuse" tags=[:integration] begin
    using LineCableModels
    using TOML, SHA, Logging
    const P=LineCableModels.PSCAD
    const launches=Ref(0)
    const mismatched_readback=Ref(false)
    const incomplete_matrix=Ref(false)
    const station_identity=Dict("schema"=>"1", "version"=>"5.1.0",
        "pscad_sha256"=>repeat("1", 64), "line_constants_sha256"=>repeat("2", 64),
        "master_library_sha256"=>repeat("3", 64))
    # Exercise the real export/stage/parse/checkpoint path with an isolated
    # transport fixture. These matrices are protocol fixtures, not references.
    function P.remote_command(::Val{:completed_run_fixture}, config::P.RemoteConfig, command::AbstractString)
        if occursin("-SharedCase", command)
            launches[]+=1
            root=replace(only(match(r"-SharedCase '([^']+)'", command).captures), '\\'=>'/')
            script="""
              using TOML
              root = $(repr(root))
              input = TOML.parsefile(joinpath(root, "computation.toml"))
              output = joinpath(root, "outputs")
              for (name, value) in (("zm", 2.0), ("zp", 30.0), ("ym", 0.001), ("yp", 75.0))
                  open(joinpath(output, "result_" * name * ".out"), "w") do io
                      println(io, "LOG10(FN) FN element")
                      for frequency in input["frequencies"]
                          println(io, log10(frequency), " ", frequency, " ", value)
                      end
                  end
              end
              open(joinpath(output, "native-settings.toml"), "w") do io
                  observed = Dict(component => Dict(field => control["readback"]
                      for (field, control) in controls) for (component, controls) in input["native_settings"])
                  $(mismatched_readback[]) && (observed["ground"]["EarthForm"] = "UNEXPECTED")
                  TOML.print(io, observed)
              end
              $(incomplete_matrix[]) && write(joinpath(output, "result_zm.out"), "LOG10(FN) FN element\n")
              write(joinpath(output, "timing.txt"), "3.25")
              write(joinpath(output, "pscad-console.txt"), "protocol fixture")
              open(joinpath(output, "solver.toml"), "w") do io
                  TOML.print(io, input["solver"])
              end
              """
        else
            script="print("*repr(sprint(TOML.print, station_identity))*")"
        end
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end

    directory=mktempdir()
    try
        copper=Material(:conductor, 1.72e-8, 1.0)
        design=build(CableDesign,
            "resume-fixture",
            terminal(:core,
                solid(copper, Disk(0.004)), insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
        system=build(LineCableSystem, [design], [Pose2(0, -1)];
            connections = [Dict(:core=>1)], system_id = basename(directory), line_length = 42.0)
        frequency=collect(10.0 .^ range(-1, 6; length = 101))
        problem=LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = frequency)
        remote=P.RemoteConfig(
            "fixture", directory, "unused", "julia", "python";
            local_root=directory, transport = :completed_run_fixture)
        options=(; remote, resume_run_directory = :latest)
        for value in (:invalid, "", 1)
            @test_throws ArgumentError computation_options(P.PSCADFormulation,
                (; remote, resume_run_directory = value))
        end
        @test_throws ArgumentError computation_options(P.PSCADFormulation, (;
            remote, solver_identity = 1))
        quiet_log=Test.TestLogger()
        first_result=with_logger(quiet_log) do
            compute(problem, Formulation(:pscad); options)
        end
        @test isempty(quiet_log.logs)
        @test launches[] == 1
        @test !details(first_result).execution.reused
        @test details(first_result).execution.elapsed_scope == P.PSCAD_TIMING_SCOPE
        source=details(first_result).execution.source_run
        files=[joinpath(root, file) for (root, _, names) in walkdir(source)
               for file in names]
        original=Dict(path=>(bytes2hex(open(sha256, path)), mtime(path)) for path in files)
        @test isfile(joinpath(source, "complete.toml"))
        verbose_log=Test.TestLogger()
        reused=with_logger(verbose_log) do
            compute(problem, Formulation(:pscad; earth_impedance = :Pollaczek1926);
                options=(; options..., verbosity=(default=0, PSCAD=1)))
        end
        @test any(record->record.message == "Exporting PSCAD computation project", verbose_log.logs)
        @test any(record->record.message == "PSCAD reuses a verified completed run", verbose_log.logs)
        @test launches[] == 1
        @test Z(reused) == Z(first_result)
        @test Y(reused) == Y(first_result)
        @test details(reused).execution.reused
        @test details(reused).execution.elapsed_seconds == 0
        @test details(reused).execution.source_elapsed_seconds == 3.25
        @test details(reused).execution.source_elapsed_scope == P.PSCAD_TIMING_SCOPE
        @test occursin("no solver execution", details(reused).execution.elapsed_scope)
        @test details(reused).formulations.requested.earth_impedance.identifier === :Pollaczek1926
        sample_log=Test.TestLogger()
        total=with_logger(sample_log) do
            LineCableModels.with_performance_sample() do
                compute(problem, Formulation(:pscad); options = (;
                    options..., output_basis = :total, verbosity=(default=0, PSCAD=2)))
            end
        end
        @test isempty(sample_log.logs)
        @test launches[] == 1
        @test Z(total) == 42 .* Z(first_result)
        @test Y(total) == 42 .* Y(first_result)
        @test original ==
              Dict(path=>(bytes2hex(open(sha256, path)), mtime(path)) for path in files)
        callbacks=Int[]
        callback_log=Test.TestLogger()
        batch=with_logger(callback_log) do
            compute(problem,
                [Formulation(:pscad), Formulation(:pscad; earth_impedance = :Pollaczek1926)];
                options = (;
                    options..., on_result = (problem, index, result)->begin
                        push!(callbacks, index)
                        @info "optional callback diagnostic"
                        @warn "visible callback warning"
                    end))
        end
        @test count(record->record.level == Logging.Warn, callback_log.logs) == 2
        @test all(record->record.level >= Logging.Warn, callback_log.logs)
        @test launches[] == 1
        @test callbacks == [1, 2]
        @test details(batch[1]).formulations.requested.earth_impedance.identifier === :default
        @test details(batch[2]).formulations.requested.earth_impedance.identifier === :Pollaczek1926
        @test all(record -> record.formula === :Pollaczek1926,
            details(batch[2]).native_setting.interactions.earth_impedance)
        @test Z(batch[1]) == Z(batch[2])
        @test Z(batch[1]) !== Z(batch[2])
        @test occursin("no solver execution", details(batch[2]).execution.elapsed_scope)
        homogeneous_choices = (air = :default, earth = :default, mixed = :default)
        represented = compute(problem,
            [Formulation(:pscad), Formulation(:pscad; earth_impedance = homogeneous_choices)]; options)
        @test launches[] == 1
        @test details(represented[1]).formulations.requested.earth_impedance.identifier === :default
        @test map(record -> record.identifier,details(represented[2]).formulations.requested.earth_impedance) == homogeneous_choices
        @test Z(represented[1]) == Z(represented[2])
        @test Z(represented[1]) !== Z(represented[2])
        changed_earth=LineParametersProblem(
            system; earth_props = homogeneous(rho = 200.0), frequencies = frequency)
        @test_throws ArgumentError compute(changed_earth, Formulation(:pscad);
            options = (; remote, resume_run_directory = source))
        changed_frequencies=LineParametersProblem(
            system; earth_props = homogeneous(rho = 100.0),
            frequencies = collect(10.0 .^ range(-1, 5; length = 101)))
        @test_throws ArgumentError compute(changed_frequencies, Formulation(:pscad);
            options = (; remote, resume_run_directory = source))

        @test_throws ArgumentError compute(
            problem, Formulation(:pscad; earth_impedance = :Saad1996);
            options = (; remote, resume_run_directory = source))
        @test_throws ArgumentError compute(problem, Formulation(:pscad);
            options = (; remote, solver_identity = Dict("version"=>"changed")))
        @test launches[] == 1
        compute(problem, Formulation(:pscad; earth_impedance = :Saad1996); options)
        @test launches[] == 2
        station_identity["line_constants_sha256"]=repeat("4", 64)
        compute(problem, Formulation(:pscad); options)
        @test launches[] == 3
        station_identity["line_constants_sha256"]=repeat("2", 64)
        # A partial run must never be accepted as a completed numerical result.
        partial=mktempdir(directory)
        cp(joinpath(source, "computation.toml"), joinpath(partial, "computation.toml"))
        @test_throws ArgumentError compute(problem, Formulation(:pscad);
            options = (; remote, resume_run_directory = partial))
        open(joinpath(source, "outputs", "result_zm.out"), "a") do io
            println(io, "damaged")
        end
        @test_throws ArgumentError compute(problem, Formulation(:pscad); options)
        @test launches[] == 3
        @test endswith(read(joinpath(source, "outputs", "result_zm.out"), String), "damaged\n")
        mismatched_readback[] = true
        @test_throws ArgumentError compute(problem, Formulation(:pscad); options = (; remote))
        @test launches[] == 4
        mismatched_readback[] = false
        incomplete_matrix[] = true
        caught = try
            compute(problem, Formulation(:pscad); options = (; remote))
            nothing
        catch error
            error
        end
        @test caught isa ErrorException
        @test occursin("PSCAD result validation failed", sprint(showerror, caught))
        @test occursin("Full PSCAD diagnostics:", sprint(showerror, caught))
        incomplete_matrix[] = false
        # Distinct physical settings in one batch require distinct native runs.
        before = launches[]
        heterogeneous = compute(problem,
            [Formulation(:pscad), Formulation(:pscad; earth_impedance = :Saad1996)];
            options = (; remote))
        @test launches[] == before + 2
        @test !details(heterogeneous[1]).execution.reused
        @test !details(heterogeneous[2]).execution.reused
        @test details(heterogeneous[1]).native_setting.ground !=
            details(heterogeneous[2]).native_setting.ground
    finally
        rm(directory; recursive = true)
    end
end
