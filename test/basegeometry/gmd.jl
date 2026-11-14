using LineCableModels.BaseGeometry: Contour, Segment, Arc, EllipticArc, gmd

@testset "Contour geometric mean distance" begin
    # Float64 bounds measured with gmd at τ = 1e-6 on Julia 1.12.7, x86_64 Linux.
    bar = Contour(Rectangle(10e-3, 5e-3))
    strip = Contour(Rectangle(10e-3, 0.5e-3, Placement(0, 2.75e-3)))
    # Maxwell, Treatise II, §692 (6), and Rosa (1907), Eq. (8), respectively.
    @test abs(log(gmd(bar, bar)) - (-5.697596499345662)) <= 1e-9
    @test abs(log(gmd(bar, strip)) - (-5.481795822821743)) <= 1e-9

    # The armored bipole's sectors 2 and 7 share an edge across the angular origin.
    # Adaptive order 15 at τ = 1e-12 gives the logarithm below.
    # On Julia 1.12.7, x86_64 Linux, gmd at τ = 1e-6 differs by 2.2916e-11.
    A = Contour(LineCableModels.Sector(0.00183, 0.0048417248992482015,
        1.0471975511965974, 0, 0, Placement(0, 0, 0)))
    B = Contour(LineCableModels.Sector(0.00183, 0.0048417248992482015,
        1.0471975511965974, 0, 0, Placement(0, 0, 5.235987755982988)))
    @test abs(log(gmd(A, B)) - (-5.690895480626689)) <= 1e-9

    @testset "curve derivatives" begin
        for curve in (Segment((1e-3, -2e-3), (4e-3, 5e-3)),
            Arc(Placement(2e-3, -1e-3, 0.4), 3e-3, 1.7),
            EllipticArc(Placement(2e-3, -1e-3, 0.4), 3e-3, 2e-3, 0.2, 1.7, 0.3e-3))
            t, h = 0.37, 1e-4
            coarse = (curve(t + h) .- curve(t - h)) ./ (2h)
            fine = (curve(t + h / 2) .- curve(t - h / 2)) ./ h
            exact = LineCableModels.BaseGeometry.derivative(curve, t)
            # The central difference's O(h²) estimate and Float64 rounding of its points.
            bound = abs.(fine .- coarse) ./ 3 .+ 8eps(Float64) * maximum(abs, curve(t)) / h
            @test all(abs.(exact .- fine) .<= bound)
        end
    end
end
