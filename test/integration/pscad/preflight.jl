@testitem "PSCAD / deterministic numeric validation precedes staging and connecting" tags=[:integration, :measurements] setup=[PSCADStandIn] begin
    using Measurements, Calculus
    const P = LineCableModels.PSCAD
    const S = PSCADStandIn
    function model(; radius=0.004, rho=1.72e-8, soil=100.0, temperature=20.0,
            frequency=10.0 .^ range(-1, 5; length=101))
        # Materialize the system at the problem's numeric type. This test targets
        # PSCAD numeric input construction, not the shared constructor's conversion methods.
        T = promote_type(typeof(radius), typeof(rho), typeof(soil), typeof(temperature), eltype(frequency))
        design = build(CableDesign, "numeric-preflight", terminal(:core,
            solid(Material(:conductor, convert(T,rho), 1.0), Disk(convert(T,radius))),
            insulation(Material(:insulator, 1e14, 2.3); t=0.002)))
        system = build(LineCableSystem, [design], [Pose2(zero(T), -one(T))]; connections=[Dict(:core=>1)])
        LineParametersProblem(system; earth_props=homogeneous(rho=soil), temperature, frequencies=frequency)
    end
    mktempdir() do directory
        root = joinpath(directory, "untouched")
        application = S.install(root)
        remote = S.remote(root)
        selected = Formulation(:pscad)
        for sigma in (0.0, 0.001)
            uncertain(x) = measurement(x, abs(x) * sigma)
            for inputs in ((radius=uncertain(0.004),), (rho=uncertain(1.72e-8),),
                    (soil=uncertain(100.0),), (temperature=uncertain(20.0),),
                    (frequency=uncertain.(10.0 .^ range(-1, 5; length=101)),))
                problem = model(; inputs...)
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
        @test_throws ArgumentError compute(model(),
            Grid((selected, Formulation(:pscad; earth_impedance=:gary1976))); options=(;remote))
        plain = model()
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
        @test S.connects(application) == 0
    end
end

@testitem "PSCAD / a start that never listens fails within the timeout" tags=[:integration, :pscad] setup=[PSCADStandIn] begin
    const P = LineCableModels.PSCAD
    const S = PSCADStandIn
    mktempdir() do directory
        application = S.install(directory; listening = false)
        touch(joinpath(S.STATION[], "never-listens"))
        remote = S.remote(directory; timeout = 5)
        started = time()
        caught = try
            P.identify(remote)
            nothing
        catch error
            error
        end
        @test caught isa ErrorException
        message = sprint(showerror, caught)
        @test occursin("PSCAD did not start on station with $(S.PSCAD_EXECUTABLE) in " *
                       "$(S.SHARED_ROOT) and serve automation on port $(S.PORT) within 5 " *
                       "seconds", message)
        @test time() - started < 20
        @test S.connects(application) == 0
        @test all(S.ended, [S.processes("forward"); S.processes("start")])
    end
end
