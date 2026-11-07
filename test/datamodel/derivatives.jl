@testset "DataModel derivatives" begin
    @testset "a centered shape carries the derivatives of its radii" begin
        @test ForwardDiff.derivative(r -> support(Circle(r)), 0.01) == 1
        @test ForwardDiff.derivative(r -> support(Annulus(5e-3, r)), 0.01) == 1
        @test ForwardDiff.derivative(r -> distance(Circle(r)), 0.01) == 0
        @test ForwardDiff.derivative(r -> distance(Annulus(r, 2e-2)), 0.01) == 1
    end

    @testset "a shape off the origin carries the derivative of its radius" begin
        at = Placement(0.05, 0.0)
        @test ForwardDiff.derivative(r -> support(Circle(r, at)), 0.01) == 1
        @test ForwardDiff.derivative(r -> distance(Circle(r, at)), 0.01) == -1
    end

    @testset "the distance of a placement is the hypot of its position" begin
        for (x, y) in ((0.0, 0.0), (-0.0, 0.0), (0.0, -0.0), (3e-3, -4e-3), (1e-300, 0.0))
            @test distance(Placement(x, y)) === hypot(x, y)
        end
    end

    # ForwardDiff 1 reads the partials of a dual in iszero.
    if pkgversion(ForwardDiff) >= v"1"
        @testset "a center moving through the origin has no derivative of its distance" begin
            @test isnan(ForwardDiff.derivative(s -> distance(Placement(s, 0.0)), 0.0))
        end
    end
end
