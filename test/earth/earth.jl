@testset "Earth" begin
    @testset "a homogeneous earth is air above one earth layer" begin
        earth = homogeneous(rho = 100, eps_r = 10)
        @test earth isa EarthModel{Float64, 2}
        air, soil = earth.layers
        @test (air.rho, air.eps_r, air.mu_r, air.thickness) == (Inf, 1.0, 1.0, Inf)
        @test (soil.rho, soil.eps_r, soil.mu_r, soil.thickness) == (100.0, 10.0, 1.0, Inf)
        @test homogeneous(rho = 100.0f0) isa EarthModel{Float32, 2}
    end

    @testset "a layered earth is not admitted" begin
        air = EarthLayer(Inf, 1, 1)
        @test_throws ArgumentError EarthModel((air, EarthLayer(100, 10, 1, 5.0)))
        @test_throws ArgumentError EarthModel((air, EarthLayer(100, 10, 1, 5.0),
            EarthLayer(500, 10, 1)))
    end
end
