using LineCableModels.BaseGeometry: Contour
using LineCableModels.DataModel: Sector

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

    @testset "a contour about the origin carries the derivatives of its radii" begin
        @test ForwardDiff.derivative(r -> support(Contour(Circle(r))), 0.01) == 1
        @test ForwardDiff.derivative(r -> distance(Contour(Annulus(1e-3, r))), 2e-3) == 0
        @test ForwardDiff.derivative(
            r -> distance(Contour(Sector(0, r, π / 3, 0, 0, Placement(0, 0)))), 10e-3) == 0
    end

    @testset "a sector with its apex on the origin has a finite derivative of its distance" begin
        @test isfinite(ForwardDiff.derivative(
            r -> distance(Sector(0, r, 2π / 3, 0, 0, Placement(0, 0))), 10e-3))
    end

    @testset "a center on the origin has a zero derivative of its distance" begin
        @test ForwardDiff.derivative(s -> distance(Placement(s, 0.0)), 0.0) == 0
    end
end
