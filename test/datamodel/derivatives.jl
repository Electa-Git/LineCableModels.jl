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

    # Julia 1.12.7 with ForwardDiff 1.4.7 on x86_64 Linux gave relative differences
    # below 1e-15 at 50 Hz.
    @testset "a tube at zero inner radius has the solid value and zero derivative" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1.0e14, 2.3)
        solid,
        tube = map((Circle(10e-3), Annulus(ForwardDiff.Dual(0.0, 1.0), 10e-3))) do shape
            design = build(CableDesign, "solid_tube",
                Terminal(:core, Region(:core_metal, shape, copper)),
                Region(:ins, Offset(8e-3), polymer))
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [(core = 1,)], id = "solid_tube", line_length = 1000.0)
            return compute(
                LineParametersProblem(system, homogeneous(rho = 100.0);
                    frequencies = [50.0]), CoaxialFormulation())
        end
        for (Q, Q0) in ((tube.Z, solid.Z), (tube.Y, solid.Y)), part in (real, imag)

            @test ForwardDiff.value(part(Q[1, 1, 1])) ≈ part(Q0[1, 1, 1]) rtol=1e-15
            @test ForwardDiff.partials(part(Q[1, 1, 1]))[1] == 0
        end
    end

    # On Julia 1.12.7, x86_64 Linux, indices 2, 3, 37, 204 and 399 gave negative
    # residues in the 400-radius scan. ForwardDiff 1.4.7 at 50 Hz and h = 1e-5r
    # gave a maximum relative error of 8.90e-10 against central differences.
    @testset "a solid core off the origin carries its radius derivative" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1.0e14, 2.3)
        for r in range(1e-3, 60e-3; length = 400)[[1, 2, 3, 37, 204, 399, 400]]
            h = 1e-5 * r
            zr, zx,
            zp,
            zm = map(((ForwardDiff.Dual(r, 1.0),
                    ForwardDiff.Dual(20e-3 + r, 1.0)),
                (r, ForwardDiff.Dual(20e-3 + r, 1.0)),
                (r + h, 20e-3 + r + h), (r - h, 20e-3 + r - h))) do (radius, x)
                design = build(CableDesign, "solid_core",
                    Terminal(:core,
                        Region(:core_metal, Circle(radius, Placement(x, 0.0)), copper)),
                    Region(:ins, Offset(8e-3), polymer))
                system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                    connections = [(core = 1,)], id = "solid_core", line_length = 1000.0)
                return compute(
                    LineParametersProblem(system, homogeneous(rho = 100.0);
                        frequencies = [50.0]),
                    CoaxialFormulation()).Z[1, 1, 1]
            end
            difference = (zp - zm) / (2h)
            for part in (real, imag)
                @test abs(ForwardDiff.partials(part(zr))[1] / part(difference) - 1) <= 1e-9
                @test ForwardDiff.partials(part(zx))[1] == 0
            end
        end
    end
end
