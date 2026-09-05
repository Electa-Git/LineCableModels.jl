@testitem "Extensions / Distributions boundary / unloaded sampler fails explicitly" tags=[
    :extension,
    :core_only
] begin
    using Random
    import LineCableModels

    @test Base.get_extension(LineCableModels, :LineCableModelsDistributionsExt) === nothing
    PB=LineCableModels.ParametricBuilder
    uncertain=only(PB.Grid(10.0, PB.AbsoluteError(1.0)))
    @test_throws ArgumentError rand(
        MersenneTwister(1),
        uncertain;
        distribution = (family = :unsupported,)
    )
end

@testitem "Extensions / Distributions adapter / normalized sampling and density laws" tags=[
    :extension
] begin
    using Distributions
    using Random
    using Statistics
    import LineCableModels

    extension_module=Base.get_extension(
        LineCableModels,
        :LineCableModelsDistributionsExt
    )
    @test extension_module !== nothing

    PB=LineCableModels.ParametricBuilder
    @test Base.ispublic(LineCableModels, :sample_uncertainty)
    @test parentmodule(LineCableModels.sample_uncertainty) === LineCableModels
    @test any(
        method -> method.module === extension_module,
        methods(LineCableModels.sample_uncertainty)
    )
    uncertain=only(PB.Grid(10.0, PB.AbsoluteError(2.0)))
    standard_draw=rand(MersenneTwister(41), uncertain; distribution = Normal(0, 1))
    shifted_draw=rand(MersenneTwister(41), uncertain; distribution = Normal(10, 3))
    @test standard_draw == shifted_draw
    @test isfinite(standard_draw)

    density=LineCableModels.UQ.HistogramDensity(
        [1.0, 3.0, 5.0],
        [0.25, 0.25]
    )
    @test pdf(density, 1.0) == 0.25
    @test pdf(density, 5.0) == 0.25
    @test pdf(density, 6.0) == 0.0
    @test logpdf(density, 6.0) == -Inf
    @test cdf(density, 2.0) == 0.25
    @test cdf(density, 4.0) == 0.75
    @test quantile(density, 0.0) == 1.0
    @test quantile(density, 1.0) == 5.0
    @test insupport(density, 3.0)
    @test !insupport(density, 6.0)
    @test moment(density, 0) == 1.0
    @test mode(density) == 2.0
    @test modes(density) == [2.0, 4.0]
    @test isapprox(mean(density), 3.0)
    @test var(density) > 0
    @test first(density.edges) <= rand(MersenneTwister(9), density) <=
          last(density.edges)
    @test_throws ArgumentError extension_module._raw_moment(density, -1)
end
