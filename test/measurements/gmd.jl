@testset "Plated bar permeability spread" begin
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    polymer = Material(:insulator, 1e14, 2.5)
    width = Measurements.measurement(10e-3, 1e-5)
    design = build(CableDesign, "plated",
        Terminal(:core,
            Stack(Region(:core_bar, Rectangle(width, 5e-3), copper),
                Region(:core_plating, Offset(0.5e-3), copper))),
        Region(:core_insulation, Offset(1e-3), polymer))
    mu = only(coaxial_input(design).material).mu_r
    @test isfinite(Measurements.uncertainty(mu)) && Measurements.uncertainty(mu) > 0
    # Julia 1.12.7 on x86_64 Linux, with gmd at τ = 1e-6, gave this nominal bound.
    @test Measurements.value(mu) ≈ 2.050 rtol = 1e-3
end
