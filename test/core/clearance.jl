@testitem "Core / clearance / sampling does not require Measurements" tags=[:core_only, :extension, :datamodel, :slow] begin
    using LineCableModels
    using Random
    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
    copper = Material(kind = :conductor, rho = 1.7e-8)
    design = build(CableDesign, "core-clearance", Group(:core,
        Region(:metal, Disk(0.01), copper)))
    space = Gridspace{LineCableSystem}(spacing -> build(LineCableSystem,
        [design, design], [Pose2(-spacing / 2, -1), Pose2(spacing / 2, -1)]),
        (Grid(0.021, AbsoluteError(0.003)),))
    # A draw whose cables overlap throws. Only a Monte Carlo study draws again.
    @test_throws OverlapError rand(MersenneTwister(4), space;
        distribution = (_rng, mean, sigma) -> mean - 3sigma)
    # A draw whose cables touch is raised to 1 μm.
    touching = @test_logs rand(MersenneTwister(4), space;
        distribution = (_rng, mean, sigma) -> mean - sigma / 3)
    @test touching.positions[2].x - touching.positions[1].x - 0.02 ≈ 1e-6 atol = 1e-14
    @test Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing
end
