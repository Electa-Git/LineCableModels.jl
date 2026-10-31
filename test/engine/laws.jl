@testset "Engine laws" begin
    Engine = LineCableModels.Engine
    select = LineCableModels.select
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)

    @testset "the linear resistivity law" begin
        law = Engine.Resistivity.Formula(Val(:linear))
        @test law(copper, copper.T0) == copper.rho
        @test law(copper, copper.T0 + 70) ≈
              copper.rho * (1 + 70 * copper.alpha) rtol = 1.0e-15
        @test_throws DomainError law(copper, copper.T0 + 150)
        @test_throws DomainError law(copper, copper.T0 - 150)
    end

    @testset "the constant properties of the earth" begin
        layer = last(homogeneous(rho = 100, eps_r = 10, mu_r = 2).layers)
        law = Engine.EarthProperties.Formula(Val(:constant))
        @test law(layer, 50.0) == layer
    end

    @testset "the deepest layer of a stack" begin
        model = Engine.EarthProperties.Formula(Val(:bottommost))
        upper = EarthLayer(100.0, 10.0, 1.0, 5.0)
        lower = EarthLayer(1000.0, 5.0, 1.0)
        @test model((upper, lower), 50.0) === lower
    end

    @testset "the lossless admittivity of insulators and semicons" begin
        f = 1.0e3
        polymer = Material(:insulator, 1.0e14, 2.3)
        semicon = Material(:semicon, 1000.0, 1000.0)
        insulation = Engine.InsulationAdmittance.Formula(Val(:lossless))
        @test insulation(polymer, f) ==
              im * 2π * f * 8.8541878128e-12 * polymer.eps_r
        semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossless))
        @test semiconducting(semicon, f) ==
              im * 2π * f * 8.8541878128e-12 * semicon.eps_r
    end

    @testset "the lossy admittivity of insulators and semicons" begin
        insulating = Engine.InsulationAdmittance.Formula(Val(:lossy))
        semiconducting = Engine.SemiconAdmittance.Formula(Val(:lossy))
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        semicon = Material(:semicon, 1000.0, 1000.0; tan_delta = 0.05)
        for f in (50.0, 1.0e3, 1.0e6),
            (law, material) in ((insulating, polymer), (semiconducting, semicon))

            ω = 2π * f
            ε = 8.8541878128e-12 * material.eps_r
            @test @inferred(law(material, f)) ≈
                  complex(1 / material.rho + ω * ε * material.tan_delta, ω * ε) rtol = 1.0e-15
        end
        polymer32 = Material(:insulator, 1.0f14, 2.3f0; tan_delta = 1.0f-3)
        semicon32 = Material(:semicon, 1000.0f0, 1000.0f0; tan_delta = 0.05f0)
        @test real(@inferred(insulating(polymer32, 50.0f0))) > 0
        @test real(@inferred(semiconducting(semicon32, 50.0f0))) > 0
    end

    @testset "the lossy admittivity without loss is the lossless one" begin
        insulator = Material(:insulator, Inf, 2.3)
        semicon = Material(:semicon, Inf, 1000.0)
        for f in (50.0, 1.0e3, 1.0e6)
            @test Engine.InsulationAdmittance.Formula(Val(:lossy))(insulator, f) ===
                  Engine.InsulationAdmittance.Formula(Val(:lossless))(insulator, f)
            @test Engine.SemiconAdmittance.Formula(Val(:lossy))(semicon, f) ===
                  Engine.SemiconAdmittance.Formula(Val(:lossless))(semicon, f)
        end
    end

    @testset "a lossy insulation gives the shunt admittance its conductance" begin
        polymer = Material(:insulator, 1.0e14, 2.3; tan_delta = 1.0e-3)
        design = build(CableDesign, "insulated_wire",
            Terminal(:core, Region(:core_metal, Circle(0.0425), copper)),
            Region(:core_insulation, Offset(1.0e-3), polymer))
        system = build(CableSystem, [design, design],
            [Placement(0.0, -1.0), Placement(1.0, -1.0)];
            connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires",
            line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 0.1);
            frequencies = 10.0 .^ range(-1, 7; length = 101))
        lossless = compute(problem, CoaxialFormulation())
        lossy = compute(problem, CoaxialFormulation(insulation_admittance = :lossy))
        @test all(>(0), real.(lossy.Y[1, 1, :]))
        @test all(real.(lossy.Y[1, 1, :]) .> real.(lossless.Y[1, 1, :]))
    end

    @testset "a slot takes an identifier, an identifier with options or a formula" begin
        Formula = Engine.Resistivity.Formula
        law = Formula(Val(:linear))
        @test select(CoaxialFormulation, Val(:resistivity), :linear) === law
        @test select(CoaxialFormulation, Val(:resistivity), :linear => (;)) === law
        @test select(CoaxialFormulation, Val(:resistivity), law) === law
        @test_throws "the formulas are (:linear,)" select(CoaxialFormulation, Val(:resistivity), :quadratic)
        @test_throws "the formulas are (:martins2026, :ideal)" select(
            CoaxialFormulation, Val(:earth_impedance), :carson1926)
        @test_throws ArgumentError select(CoaxialFormulation, Val(:resistivity), :linear => (rate = 1.0,))
    end
end
