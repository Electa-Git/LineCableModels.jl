@testset "Earth" begin
    @testset "a homogeneous earth is air above one earth layer" begin
        earth = homogeneous(rho = 100, eps_r = 10)
        @test earth isa EarthModel{Float64, 2}
        air, soil = earth.layers
        @test (air.rho, air.eps_r, air.mu_r, air.thickness) == (Inf, 1.0, 1.0, Inf)
        @test (soil.rho, soil.eps_r, soil.mu_r, soil.thickness) == (100.0, 10.0, 1.0, Inf)
        @test homogeneous(rho = 100.0f0) isa EarthModel{Float32, 2}
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
    @testset "finite relative permittivity of an equivalent earth" begin
        @test EarthLayer(52.4, -300334.5, 1.0).eps_r == -300334.5
        @test_throws DomainError EarthLayer(100.0, Inf, 1.0)
        @test_throws DomainError EarthLayer(100.0, NaN, 1.0)
    end
end
