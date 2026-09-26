@testitem "PSCAD / deterministic numeric validation precedes staging and transport" tags=[:integration] begin
    using Measurements, Calculus, Serialization
    const P = LineCableModels.PSCAD
    const contacts = Ref(0)
    function P.remote_command(::Val{:preflight_probe}, ::P.RemoteConfig, ::AbstractString)
        contacts[] += 1
        error("preflight contacted transport")
    end
    function model(; radius=0.004, rho=1.72e-8, soil=100.0, temperature=20.0,
            frequency=10.0 .^ range(-1, 5; length=101))
        # Materialize the system at the problem's numeric type; this test targets
        # PSCAD numeric preparation, not the shared constructor's conversion methods.
        T = promote_type(typeof(radius), typeof(rho), typeof(soil), typeof(temperature), eltype(frequency))
        design = build(CableDesign, "numeric-preflight", terminal(:core,
            solid(Material(:conductor, convert(T,rho), 1.0), Disk(convert(T,radius))),
            insulation(Material(:insulator, 1e14, 2.3); t=0.002)))
        system = build(LineCableSystem, [design], [Pose2(zero(T), -one(T))]; connections=[Dict(:core=>1)])
        LineParametersProblem(system; earth_props=homogeneous(rho=soil), temperature, frequencies=frequency)
    end
    mktempdir() do directory
        root = joinpath(directory, "untouched")
        remote = P.RemoteConfig("fixture", "shared", "scratch", "julia", "python";
            local_root=root, transport=:preflight_probe)
        selected = Formulation(:pscad)
        for sigma in (0.0, 0.001)
            uncertain(x) = measurement(x, abs(x) * sigma)
            for inputs in ((radius=uncertain(0.004),), (rho=uncertain(1.72e-8),),
                    (soil=uncertain(100.0),), (temperature=uncertain(20.0),),
                    (frequency=uncertain.(10.0 .^ range(-1, 5; length=101)),))
                problem = model(; inputs...)
                @test_throws r"PSCAD does not support Measurement" validate(problem, selected)
                @test_throws r"PSCAD does not support Measurement" compute(problem, selected; options=(;remote))
                @test !ispath(root)
                @test_throws r"PSCAD does not support Measurement" export_data(:pscad,
                    problem.system, problem.earth_props; temperature=problem.temperature,
                    file_name=joinpath(root, "invalid.pscx"))
            end
            @test_throws r"PSCAD does not support Measurement" Formulation(:pscad;
                options=(base_frequency=uncertain(50.0),))
            plain = model()
            @test_throws r"PSCAD does not support Measurement" export_data(:pscad,
                plain.system, plain.earth_props; base_freq=uncertain(50.0),
                file_name=joinpath(root, "invalid.pscx"))
        end
        @test_throws ArgumentError compute(model(frequency=[50.0]), selected; options=(;remote))
        @test_throws ArgumentError compute(model(), [selected, Formulation(:pscad; earth_impedance=:gary1976)]; options=(;remote))
        plain = model()
        snapshot(value) = (io=IOBuffer(); serialize(io,value); take!(io))
        original = snapshot(plain)
        @test validate(plain,selected) === plain
        @test snapshot(plain) == original
        metal = Material(:conductor,1.72e-8,1.0)
        dielectric = Material(:insulator,1e14,2.3)
        design = build(CableDesign,"five-components",(terminal(Symbol("component$i"),
            i == 1 ? core(metal;r=0.004) : sheath(metal;t=0.001),
            insulation(dielectric;t=0.002)) for i in 1:5)...)
        system = build(LineCableSystem,[design],[Pose2(0,-1)];
            connections=[Dict(name=>i for (i,name) in enumerate(design.terminal_order))])
        excessive = LineParametersProblem(system;earth_props=plain.earth_props,
            frequencies=plain.frequencies)
        @test_throws ArgumentError compute(excessive,selected;options=(;remote))
        @test !ispath(root)
        @test contacts[] == 0
    end
end

@testitem "PSCAD / caller deadline terminates a silent transport" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    function P.remote_command(::Val{:deadline_probe}, ::P.RemoteConfig, ::AbstractString)
        `$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e 'sleep(30)'`
    end
    mktempdir() do directory
        remote = P.RemoteConfig("fixture", "shared", "scratch", "julia", "python";
            local_root=directory, transport=:deadline_probe)
        started = time()
        @test_throws r"transport exceeded its timeout" P._run_remote(remote, "unused"; timeout_seconds=0.5)
        @test time() - started < 10
    end
end
