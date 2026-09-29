@testitem "Engine / analytical earth impedances / source equations and scalar types" tags=[:unit] begin
    const E = LineCableModels.Engine
    # Independently evaluated with mpmath at 80 decimal digits from Saad (5)/(6),
    # Wedepohl-Wilcox (7)/(8), and PSCAD (8-28)/(8-34). All lengths are in m.
    fixtures = (
        (:saad1996, :self, (-1, -1), 0, (2, 2), 1,
            "0.000049436050542009708503711029366076531443284669219897", "0.00071885090587064612822779262048186031909301811589599"),
        (:saad1996, :mutual, (-1, -2), 0.75, (2, 2), 1,
            "0.000049479207377905126939332583978927107226053797161543", "0.00041543494490419351553063029561713015092742757804909"),
        (:saad1996, :mutual, (-1, -2), 0, (2, 2), 1,
            "0.000049479424322246073989373455794913818704375130284438", "0.00042945544041135657926423194846619726345300197040710"),
        (:wedepohl1973, :self, (-1, -1), 0, (2, 2), 2.5,
            "0.00012383531262692624137104164156086318287741590034191", "0.0017249149591946593457986273924888321274339160401549"),
        (:wedepohl1973, :mutual, (-1, -2), 0.75, (2, 2), 2.5,
            "0.00012406794143358087068884689359206882972041322921760", "0.00096625258207122531418241849265399990355136359011624"),
        (:gary1976, :self, (10, 10), 0, (1, 1), 1,
            "0.000048477544121804667108869118283828124463065274088964", "0.00072467343587872156465893473991564770639800889262178"),
        (:gary1976, :mutual, (10, 12), 0.75, (1, 1), 1,
            "0.000048391816886193307359343300382265063290248575131592", "0.00038772498399398139118419994348600879217263733480040"),
        (:lucca1994, :mutual, (10, -2), 0.75, (1, 2), 1,
            "0.000048694756966967542207217998198776082068243035825963", "0.00027361164043745607069020182417817933205189613540840"),
        (:lucca1994, :mutual, (-2, 10), 0.75, (2, 1), 1,
            "0.000048694756966967542207217998198776082068243035825963", "0.00027361164043745607069020182417817933205189613540840"))
    function evaluate(functor)
        value = @inferred functor()
        equation = functor.binding.equation
        @test (@inferred equation(functor, functor.binding.pair, nothing)) == value
        return value
    end
    for T in (Float32, Float64, BigFloat),
        (id, kind, heights, x, layers, μr, re, im) in fixtures

        μ0 = 4T(pi) * T(10)^(-7)
        ε0 = parse(T, "8.8541878128e-12")
        s = complex(zero(T), 100T(pi))
        pair = E.EarthPair(1, kind === :self ? 1 : 2, T.(heights), T(x), layers;
            radius = kind === :self ? T(1)/100 : nothing)
        functor = E.EarthImpedance.Formula(id)(T[Inf, 100], [ε0, 10ε0],
            [μ0, T(μr)*μ0], s, pair)
        value = evaluate(functor)
        @test value isa Complex{T}
        @test isfinite(value)
        T === Float32 || (@test value ≈ complex(parse(T, re), parse(T, im)) rtol=5e-14)
    end
    # Changing radius affects only the explicit self distance for this expression.
    μ0, ε0, s = 4pi*1e-7, 8.8541878128e-12, 100pi*im
    gary = E.EarthImpedance.Formula(:gary1976)
    self(r) = gary([Inf, 100.0], [ε0, 10ε0], [μ0, μ0], s,
        E.EarthPair(1, 1, (10.0, 10.0), 0.0, (1, 1); radius = r))()
    @test self(0.02)-self(0.01) ≈ -s*μ0/(2pi)*log(2) rtol=5e-14
end

@testitem "Engine / analytical earth impedances / uncertainty dependencies" tags=[:unit] begin
    using Measurements
    const E = LineCableModels.Engine
    μ0, ε0, s = 4pi*1e-7, 8.8541878128e-12, 100pi*im
    function evaluate(id, rho, heights, layers)
        pair = E.EarthPair(1, 1, heights, 0.0, layers; radius = 0.01)
        id === :lucca1994 && (pair = E.EarthPair(1, 2, heights, 0.75, layers))
        selected = E.EarthImpedance.Formula(id)
        T = typeof(rho)
        functor = selected(T[Inf, rho], T[ε0, 10ε0], T[μ0, μ0], Complex{T}(s), pair)
        return @inferred functor()
    end
    rho = measurement(100.0, 2.0)
    for (id, heights, layers) in ((:saad1996, (-1.0, -1.0), (2, 2)),
        (:wedepohl1973, (-1.0, -1.0), (2, 2)),
        (:gary1976, (10.0, 10.0), (1, 1)), (:lucca1994, (10.0, -2.0), (1, 2)))
        z = evaluate(id, rho, heights, layers)
        shared = evaluate(id, rho, heights, layers)
        independent = evaluate(id, measurement(100.0, 2.0), heights, layers)
        @test z isa Complex{Measurement{Float64}}
        @test Measurements.value(real(z)) ≈ real(evaluate(id, 100.0, heights, layers))
        @test Measurements.value(imag(z)) ≈ imag(evaluate(id, 100.0, heights, layers))
        @test iszero(uncertainty(real(z-shared)))
        @test iszero(uncertainty(imag(z-shared)))
        @test uncertainty(real(z-independent)) > 0
        @test uncertainty(imag(z-independent)) > 0
    end
    # Gary self: p ∝ sqrt(rho), dZ/drho = (s μ0 / 2π) p / (2rho(h+p)).
    z = evaluate(:gary1976, rho, (10.0, 10.0), (1, 1))
    p = sqrt(100/(s*μ0))
    sensitivity = s*μ0/(2pi) * p/(200*(10+p))
    @test Measurements.derivative(real(z), rho) ≈ real(sensitivity) rtol=2e-13
    @test Measurements.derivative(imag(z), rho) ≈ imag(sensitivity) rtol=2e-13
end

@testitem "Engine / analytical earth impedances / public calculations and material transport" tags=[:unit] setup=[FormulaFixtures] begin
    using LinearAlgebra
    const E=LineCableModels.Engine
    metal=Material(:conductor, 1.72e-8, 1.0)
    dielectric=Material(:insulator, Inf, 2.3)
    design=build(CableDesign, "analytical-earth",
        terminal(:core,
            core(metal; r = 0.004), insulation(dielectric; t = 0.002)))
    options=(reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    for heights in ((10.0, 12.0), (-1.0, -2.0), (10.0, -2.0))
        system=build(LineCableSystem, fill(design, 2),
            [Pose2(0, heights[1]), Pose2(0.75, heights[2])];
            connections = [Dict(:core=>1), Dict(:core=>2)])
        problem=LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
            frequencies = [50.0, 500.0])
        for buried in (:saad1996, :wedepohl1973)
            selected=Formulation(
                earth_impedance = (air = :gary1976, earth = buried, mixed = :lucca1994),
                earth_admittance = :ideal; options)
            result=compute(problem, selected; options = (trace = true,))
            trace=details(result).data.trace
            @test size(Z(result)) == size(Y(result)) == (2, 2, 2)
            @test result.f == problem.frequencies
            @test all(isfinite, Z(result)) && all(isfinite, Y(result))
            @test occursin("Gary", description(selected, Z))
            for k in 1:2
                @test Z(result)[:, :, k] ≈ trace.Z[:, :, k]
                @test Y(result)[:, :, k] ≈
                      2pi*problem.frequencies[k]*im*inv(trace.P[:, :, k])
                @test trace.P[:, :, k] ≈ trace.Pin[:, :, k]+trace.Pg[:, :, k]
            end
            # The existing frequency loop supplies the law once per soil/frequency.
            fd=FormulaFixtures.DispersiveEarth()
            dispersive=Formulation(earth_impedance = selected.methods.earth_impedance,
                earth_admittance = :ideal, earth_properties = fd; options)
            actual=compute(problem, dispersive; options = (trace = true,))
            @test length(fd.seen) == 2
            for (k, frequency) in enumerate(problem.frequencies)
                reference_problem=LineParametersProblem(system;
                    earth_props = homogeneous(rho = 100/(1+frequency/100)), frequencies = [frequency])
                reference=compute(reference_problem, selected)
                @test Z(actual)[:, :, k] ≈ Z(reference)[:, :, 1]
                @test Y(actual)[:, :, k] ≈ Y(reference)[:, :, 1]
            end
        end
        # One Unified quantity still prepares its complete system when the other
        # quantity uses a direct equation. Compare each quantity to its full run.
        full=compute(problem, Formulation(; options))
        unified_Z=compute(problem, Formulation(earth_impedance = :unified,
            earth_admittance = :ideal; options))
        unified_P=compute(problem,
            Formulation(
                earth_impedance = (air = :gary1976, earth = :saad1996, mixed = :lucca1994),
                earth_admittance = :unified; options))
        @test Z(unified_Z) == Z(full)
        @test Y(unified_P) == Y(full)
        if all(<(0), heights)
            partial=Formulation(earth_impedance = (earth = :saad1996,),
                earth_admittance = (earth = :ideal,); options)
            reference=compute(problem, partial)
            space=Formulation(earth_impedance = Grid((:saad1996, :wedepohl1973)),
                earth_admittance = :ideal; options)
            results=compute(problem, space)
            @test length(results) == 2
            @test Z(results[1]) == Z(reference)
            @test Y(results[1]) == Y(reference)
            @test Z(results[2]) == Z(compute(problem,
                Formulation(earth_impedance = :wedepohl1973, earth_admittance = :ideal; options)))
        end
    end
end
