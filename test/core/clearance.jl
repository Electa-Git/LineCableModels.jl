@testitem "Core / clearance / sampling does not require Measurements" tags=[:core_only] begin
    using LineCableModels
    using Random
    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
    copper = Material(kind = :conductor, rho = 1.7e-8)
    design = build(CableDesign, "core-clearance", Group(:core,
        Region(:metal, Disk(0.01), copper)))
    space = Gridspace{LineCableSystem}(spacing -> build(LineCableSystem,
        [design, design], [Pose2(-spacing / 2, -1), Pose2(spacing / 2, -1)]),
        (Grid(0.021, AbsoluteError(0.003)),))
    sampled = @test_logs (:warn, r"Sampled cable placements adjusted") rand(
        MersenneTwister(4), space; distribution = (_rng, mean, sigma) -> mean - 3sigma)
    @test sampled.positions[2].x - sampled.positions[1].x - 0.02 >= 1e-6
    @test sampled.clearances[1, 2] == 1e-6
    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
end
