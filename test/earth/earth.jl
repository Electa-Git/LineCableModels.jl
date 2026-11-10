@testset "Earth" begin
    @testset "a homogeneous earth is air above one earth layer" begin
        earth = homogeneous(rho = 100, eps_r = 10)
        @test earth isa EarthModel{Float64, 2}
        air, soil = earth.layers
        @test (air.rho, air.eps_r, air.mu_r, air.thickness) == (Inf, 1.0, 1.0, Inf)
        @test (soil.rho, soil.eps_r, soil.mu_r, soil.thickness) == (100.0, 10.0, 1.0, Inf)
    end

    @testset "a layered earth is a stack from the air to a semi-infinite layer" begin
        air = EarthLayer(Inf, 1, 1)
        upper = EarthLayer(100, 10, 1, 5.0)
        bottom = EarthLayer(500, 10, 1)
        @test EarthModel((air, upper, bottom)).layers === (air, upper, bottom)
        @test_throws "must be air and at least one earth layer" EarthModel((air,))
        @test_throws "layers[1] must be semi-infinite air" EarthModel((upper, bottom))
        @test_throws "layers[2] must be semi-infinite" EarthModel((air, upper))
        @test_throws "layers[3] must be semi-infinite" EarthModel((air, upper, upper))
        @test_throws "finite thickness" EarthModel((air, bottom, bottom))
        @test_throws "EarthLayer.thickness must be positive" EarthLayer(100, 10, 1, 0.0)
    end
    @testset "the layers of a block build an earth model" begin
        model = @layers begin
            EarthLayer(Inf, 1.0, 1.0)
            EarthLayer(100.0, 10.0, 1.0, 5.0)
            EarthLayer(1000.0, 10.0, 1.0)
        end
        @test model == EarthModel((EarthLayer(Inf, 1.0, 1.0),
            EarthLayer(100.0, 10.0, 1.0, 5.0), EarthLayer(1000.0, 10.0, 1.0)))
    end
    @testset "a declared earth model takes a relative permittivity of at least 1" begin
        air = EarthLayer(Inf, 1.0, 1.0)
        low = EarthLayer(52.4, -300334.5, 1.0)
        @test low.eps_r == -300334.5
        @test_throws "EarthModel.layers[2].eps_r must be at least 1" EarthModel((air, low))
        @test_throws "EarthModel.layers[2].eps_r must be at least 1" EarthModel((
            air, EarthLayer(100.0, 0.5, 1.0)))
        @test_throws "EarthModel.layers[1].eps_r must be at least 1" EarthModel((
            EarthLayer(Inf, 0.5, 1.0), EarthLayer(100.0, 10.0, 1.0)))
        @test_throws DomainError EarthLayer(100.0, Inf, 1.0)
        @test_throws DomainError EarthLayer(100.0, NaN, 1.0)
    end
    @testset "a computed medium takes any resistivity but NaN" begin
        @test EquivalentEarthLayer(-393, 333 // 2, 1) ===
              EquivalentEarthLayer{Float64}(-393.0, 166.5, 1.0)
        @test EquivalentEarthLayer(Inf, 1.0, 1.0).rho == Inf
        @test EquivalentEarthLayer(52.4, -300334.5, 1.0).eps_r == -300334.5
        @test_throws "EquivalentEarthLayer.rho must not be NaN" EquivalentEarthLayer(
            NaN, 1.0, 1.0)
        @test_throws "EquivalentEarthLayer.eps_r must be finite" EquivalentEarthLayer(
            100.0, Inf, 1.0)
        @test_throws "EquivalentEarthLayer.eps_r must be finite" EquivalentEarthLayer(
            100.0, NaN, 1.0)
        @test_throws "EquivalentEarthLayer.mu_r must be positive and finite" EquivalentEarthLayer(
            100.0, 1.0, 0.0)
        @test_throws "EquivalentEarthLayer.mu_r must be positive and finite" EquivalentEarthLayer(
            100.0, 1.0, Inf)
    end
end
