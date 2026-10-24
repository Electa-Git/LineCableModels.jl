@testset "Materials" begin
    @testset "the class is the type parameter" begin
        @test Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393) isa
              Material{Float64, :conductor}
        @test Material(:semicon, 1000.0, 1000.0) isa Material{Float64, :semicon}
        @test material(kind = :insulator, rho = Inf, eps_r = 2.3) ===
              Material(:insulator, Inf, 2.3)
    end

    @testset "Float32 and mixed inputs" begin
        @test Material(:conductor, 1.7241f-8, 1, 1, 20, 0.00393f0) isa
              Material{Float32, :conductor}
        @test Material(:conductor, 1.7241f-8, 1, 1, 20, 0.00393) isa
              Material{Float64, :conductor}
    end
end
