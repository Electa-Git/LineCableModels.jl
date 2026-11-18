@testset "Plated bar permeability derivative" begin
    w, h = 10e-3, 1e-6
    values = map((ForwardDiff.Dual(w, 1.0), w + h, w - h)) do width
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, 1e14, 2.5)
        design = build(CableDesign, "plated",
            Terminal(:core,
                Stack(Region(:core_bar, Rectangle(width, 5e-3), copper),
                    Region(:core_plating, Offset(0.5e-3), copper))),
            Region(:core_insulation, Offset(1e-3), polymer))
        only(coaxial_input(design).material).mu_r
    end
    fine = (values[2] - values[3]) / (2h)
    partial = only(ForwardDiff.partials(values[1]))
    # The partial has a 1e-6 relative bound, at the scale of τ, which bounds ln GMD.
    # On Julia 1.12.7, x86_64 Linux, gmd at τ = 1e-6 gave a gap of 2.2533e-7
    # for a partial of 89.6856562353.
    @test abs(partial - fine) <= 1e-6 * abs(fine)
    @test ForwardDiff.value(values[1]) ≈ 2.050 rtol = 1e-3
end
