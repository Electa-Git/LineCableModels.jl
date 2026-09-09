@testitem "PSCAD / transport failure retains diagnostics and interruption cancels" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    using Logging
    function P.remote_command(::Val{:local_failure_probe}, config::P.RemoteConfig, command::AbstractString)
        script = occursin("-SharedCase", command) ?
            "println(stderr, \"synthetic transport failure\"); exit(9)" :
            "write(" * repr(joinpath(config.local_root, "cancelled")) * ", \"cancelled\")"
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
    function P.remote_command(::Val{:local_interrupt_probe}, config::P.RemoteConfig, command::AbstractString)
        script = "write(" * repr(joinpath(config.local_root, "started")) *
            ", \"started\"); println(\"running\"); println(stderr, \"diagnostic\"); flush(stdout); flush(stderr); sleep(30)"
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
    mktempdir() do root
        config = P.RemoteConfig("local-fixture", root, "scratch", "julia", "python";
            local_root = root, transport = :local_failure_probe)
        directory = mkpath(joinpath(root, "nested", "case"))
        project = joinpath(directory, "generated.pscx")
        write(project, "protocol fixture")
        output = joinpath(directory, "outputs")
        caught = try
            P.run_remote_pscad(config, project, output, Formulation(:pscad),
                10.0.^range(-1, 5; length = 101); output_stem = "fixture")
            nothing
        catch error
            error
        end
        @test caught isa ErrorException
        @test occursin("Full PSCAD diagnostics:", sprint(showerror, caught))
        @test occursin("Transport stderr:", sprint(showerror, caught))
        @test read(joinpath(output, "transport-stderr.txt"), String) == "synthetic transport failure\n"
        @test isfile(joinpath(root, "cancelled"))
        @test occursin("before the remote failure", read(joinpath(output, "pscad-console.txt"), String))
        @test_throws ArgumentError P.run_remote_pscad(config, project, output,
            Formulation(:pscad), 10.0.^range(-1, 5; length = 101); output_stem = "fixture")
        missing = joinpath(root, "missing.log")
        @test P._diagnostic_tail(missing) == "PSCAD produced no diagnostic log."
        write(missing, "\n   \n")
        @test P._diagnostic_tail(missing) == "PSCAD diagnostic log is empty."
        write(missing, "one\n\ntwo\nthree\n")
        @test P._diagnostic_tail(missing; count = 2) == "two\nthree"
    end
    for callback_fails in (false, true)
        mktempdir() do root
            config = P.RemoteConfig("local-fixture", root, "scratch", "julia", "python";
                local_root = root, transport = :local_interrupt_probe)
            cancelled = Ref(false)
            callback = () -> begin
                cancelled[] = true
                callback_fails && error("synthetic cancellation diagnostic")
            end
            log = Test.TestLogger()
            task = with_logger(log) do
                @async try
                    P._run_remote(config, "ignored";
                        stdout_path = joinpath(root, "out.log"),
                        stderr_path = joinpath(root, "err.log"),
                        stream = true, on_interrupt = callback)
                catch error
                    error
                end
            end
            ready = timedwait(() -> isfile(joinpath(root, "started")) || istaskdone(task), 15)
            @test ready === :ok
            @test !istaskdone(task)
            schedule(task, InterruptException(); error = true)
            @test fetch(task) isa InterruptException
            @test cancelled[]
            @test callback_fails == any(record -> occursin("cancellation could not be confirmed",
                string(record.message)), log.logs)
        end
    end
end

@testitem "PSCAD / unsupported indexed equations fail without fallback" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    for equation in (P.earth_impedance, P.earth_potential_coefficient),
            (kind, s, t) in ((:self, 1, 2), (:mutual, 2, 3), (:self, 3, 3))
        caught = try
            equation(Val(:default), Val(kind), Val(s), Val(t), Val(:pscad))
            nothing
        catch error
            error
        end
        @test caught isa ArgumentError
        @test occursin("source in layer $s and target in layer $t", sprint(showerror, caught))
    end
    @test_throws ArgumentError P.internal_impedance(Val(:default), Val(:invalid), Val(:pscad))
    @test_throws ArgumentError P.insulation_impedance(Val(:not_registered), Val(:pscad))
    for (s, t) in ((1, 2), (2, 1))
        @test P.earth_impedance(Val(:Ametani2009), Val(:mutual), Val(s), Val(t), Val(:pscad)) ==
            (EarthForm3 = (value = 0, readback = "AMETANIL"),)
        @test P.earth_impedance(Val(:Lucca1994), Val(:mutual), Val(s), Val(t), Val(:pscad)) ==
            (EarthForm3 = (value = 2, readback = "LUCCA"),)
    end
    const pipe = LineCableModels.Engine.PipeImpedance.Formula(:default)
    @test_throws ArgumentError Formulation(Val(:pscad), Val(:default), pipe, Val(:pipe))
end

@testitem "Mixed analytical selections / reciprocity and restricted applicability" tags=[:integration] begin
    const E = LineCableModels.Engine
    rho = [Inf, 100.0]
    epsilon = 8.8541878128e-12 .* [1, 10]
    mu = fill(4pi * 1e-7, 2)
    selected = E.EarthImpedance.Formula(:Ametani2009)
    @test occursin("mixed", description(selected))
    forward = E.EarthPair(1, 2, (2.0, -1.0), 0.75, (1, 2))
    reverse = E.EarthPair(2, 1, (-1.0, 2.0), 0.75, (2, 1))
    for frequency in (0.1, 50.0, 1e5)
        jω = complex(0.0, 2pi * frequency)
        z = selected(rho, epsilon, mu, jω, forward)()
        reversed = selected(rho, epsilon, mu, jω, reverse)()
        @test isfinite(z) && !iszero(z)
        @test z ≈ reversed rtol = 5eps(Float64)
    end
    for pair in (E.EarthPair(1, 1, (2.0, 2.0), 0.0, (1, 1); radius = 0.01),
            E.EarthPair(1, 2, (-1.0, -2.0), 0.75, (2, 2)))
        @test_throws ArgumentError validate(selected, pair)
    end
    copper = Material(:conductor, 1.72e-8, 1.0)
    design = build(CableDesign, "mixed-order", terminal(:core,
        solid(copper, Disk(0.004)), insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
    choices = (air = :Carson1926, earth = :Pollaczek1926, mixed = :Ametani2009)
    configuration = Formulation(earth_impedance = choices)
    function problem(poses)
        system = build(LineCableSystem, [design, design], poses;
            connections = [Dict(:core => 1), Dict(:core => 2)])
        LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    end
    # No mixed default Y is registered: a Z choice cannot authorize a Y fallback.
    @test_throws ArgumentError compute(problem([Pose2(0, 2), Pose2(0.75, -1)]), configuration)
    @test_throws ArgumentError compute(problem([Pose2(0.75, -1), Pose2(0, 2)]), configuration)
end
