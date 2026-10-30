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

    @testset "a slot takes an identifier, an identifier with options or a formula" begin
        Formula = Engine.Resistivity.Formula
        law = Formula(Val(:linear))
        @test select(Formula, :linear) === law
        @test select(Formula, :linear => (;)) === law
        @test select(Formula, law) === law
        @test_throws "the formulas are (:linear,)" select(Formula, :quadratic)
        @test_throws "the formulas are (:unified,)" select(
            Engine.EarthReturn.Formula, :carson1926)
        @test_throws ArgumentError select(Formula, :linear => (rate = 1.0,))
    end
end
