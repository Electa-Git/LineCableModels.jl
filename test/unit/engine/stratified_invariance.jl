@testitem "Engine / stratified earth / artificial interfaces and reciprocity" tags=[:unit] begin
    const EN = LineCableModels.Engine
    const EI = EN.EarthImpedance
    epsilon0 = 8.8541878128e-12
    mu0 = 4pi * 1e-7
    selected = EI.Formula(:Tsiamitros2008)

    # An interface between identical media has no physical effect. Exercise
    # same-layer, cross-layer, overhead and mixed pairs, including a source in
    # the infinite bottom layer and interfaces between both source and target.
    placements = (
        ((2.0, 3.0), (1, 1)),
        ((-0.5, -3.5), (2, 5)),
        ((-1.5, -2.5), (3, 4)),
        ((-1.5, -1.8), (3, 3)),
        ((2.0, -3.5), (1, 5)),
    )
    for frequency in (50.0, 1000.0)
        s = complex(0.0, 2pi * frequency)
        homogeneous = @inferred selected([Inf, 100.0], epsilon0 .* [1.0, 12.0],
            mu0 .* [1.0, 1.0], s, nothing, nothing, [Inf, Inf])
        subdivided = @inferred selected([Inf, 100.0, 100.0, 100.0, 100.0],
            epsilon0 .* [1.0, 12.0, 12.0, 12.0, 12.0], fill(mu0, 5),
            s, nothing, nothing, [Inf, 1.0, 1.0, 1.0, Inf])
        for (heights, layers) in placements
            whole_pair = EN.EarthPair(1, 2, heights, 0.75,
                map(layer -> layer == 1 ? 1 : 2, layers))
            split_pair = EN.EarthPair(1, 2, heights, 0.75, layers)
            reciprocal = EN.EarthPair(2, 1, reverse(heights), 0.75, reverse(layers))
            expected = @inferred homogeneous(Val(:mutual), whole_pair)
            actual = @inferred subdivided(Val(:mutual), split_pair)
            @test isfinite(actual)
            @test actual ≈ expected rtol=5e-7
            @test subdivided(Val(:mutual), reciprocal) ≈ actual rtol=1e-12
        end
        for (depth, layer) in ((0.5, 2), (1.5, 3), (3.5, 5))
            whole_self = EN.EarthPair(1, 1, (-depth, -depth), 0.02, (2, 2))
            split_self = EN.EarthPair(1, 1, (-depth, -depth), 0.02, (layer, layer))
            @test subdivided(Val(:self), split_self) ≈
                homogeneous(Val(:self), whole_self) rtol=5e-7
        end
    end

    # Subdivide an internal layer while retaining genuine material contrasts
    # above and below it. This exercises reflected fields, not only zero
    # reflection coefficients in the homogeneous limit.
    s = complex(0.0, 2pi * 50.0)
    unsplit = selected([Inf, 10.0, 250.0, 50.0],
        epsilon0 .* [1.0, 5.0, 15.0, 8.0], mu0 .* [1.0, 1.1, 1.3, 1.05],
        s, nothing, nothing, [Inf, 1.0, 2.0, Inf])
    split = selected([Inf, 10.0, 250.0, 250.0, 50.0],
        epsilon0 .* [1.0, 5.0, 15.0, 15.0, 8.0],
        mu0 .* [1.0, 1.1, 1.3, 1.3, 1.05],
        s, nothing, nothing, [Inf, 1.0, 1.0, 1.0, Inf])
    for (heights, layers) in placements
        merged_layers = map(layer -> layer == 4 ? 3 : layer == 5 ? 4 : layer, layers)
        original = EN.EarthPair(1, 2, heights, 0.75, merged_layers)
        subdivided = EN.EarthPair(1, 2, heights, 0.75, layers)
        reciprocal = EN.EarthPair(2, 1, reverse(heights), 0.75, reverse(layers))
        expected = unsplit(Val(:mutual), original)
        actual = split(Val(:mutual), subdivided)
        @test isfinite(actual)
        @test actual ≈ expected rtol=5e-7
        @test split(Val(:mutual), reciprocal) ≈ actual rtol=1e-12
    end
end
