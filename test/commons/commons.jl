@testset "Commons" begin
    @testset "validate returns what it checks" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        @test LineCableModels.validate(copper) === copper
        earth = homogeneous(rho = 100, eps_r = 10)
        @test LineCableModels.validate(earth) === earth
    end
end
