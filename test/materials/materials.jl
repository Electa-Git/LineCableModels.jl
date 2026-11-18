@testset "Materials" begin
    @testset "the class is the type parameter" begin
        @test Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393) isa
              Material{Float64, :conductor}
        @test Material(:semicon, 1000.0, 1000.0) isa Material{Float64, :semicon}
        @test material(kind = :insulator, rho = Inf, eps_r = 2.3) ===
              Material(:insulator, Inf, 2.3)
        @test_throws DomainError Material(:conductr, 1.7241e-8)
    end

    @testset "dielectrics side by side combine by their fractions of the area" begin
        polyester = material(kind = :insulator, rho = 1.0e15, eps_r = 3.3, tan_delta = 2.0e-3)
        compound = material(kind = :insulator, rho = 1.0e13, eps_r = 2.2, tan_delta = 5.0e-4)
        semicon = Material(:semicon, 1000.0, 1000.0)
        composite = Material([polyester, compound], [0.25, 0.75])
        @test composite isa Material{Float64, :insulator}
        @test 1 / composite.rho ≈ 0.25 / 1.0e15 + 0.75 / 1.0e13
        @test composite.eps_r ≈ 0.25 * 3.3 + 0.75 * 2.2
        @test composite.eps_r * composite.tan_delta ≈
              0.25 * 3.3 * 2.0e-3 + 0.75 * 2.2 * 5.0e-4
        @test Material([polyester, semicon], [0.5, 0.5]) isa Material{Float64, :semicon}
        @test Material([compound], [1.0]) === compound
    end
end
