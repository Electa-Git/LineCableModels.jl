@testitem "ParametricBuilder / wire patterns / stranded geometry invariants" tags = [:unit] begin
    const Wires = LineCableModels.ParametricBuilder.WirePatterns

    for awg in (-3, 0, 12, 40)
        diameter = Wires.awg_to_d_mm(awg)
        area = Wires.awg_to_area_mm2(awg)
        @test Wires.d_mm_to_awg(diameter) ≈ awg atol = 64eps(Float64)
        @test Wires.area_mm2_to_awg(area) ≈ awg atol = 64eps(Float64)
    end
    @test Wires.awg_label(-3) == "0000 (4/0)"
    @test Wires.awg_label(12) == "12"

    target_mm2 = 1_000.0
    choices = make_stranded(target_mm2)
    for choice in values(choices)
        @test choice.wires == 1 + 3 * choice.layers * (choice.layers - 1)
        @test choice.total_area_m2 > 0
        @test choice.wire_diameter_m > 0
    end
    @test abs(1e6 * choices.best_match.total_area_m2 - target_mm2) <=
          abs(1e6 * choices.min_diam.total_area_m2 - target_mm2)
    @test choices.min_layers.layers <= choices.min_diam.layers
    @test_throws AssertionError make_stranded(0.0)
    @test_throws AssertionError make_stranded(10.0; nmin = 10, nmax = 9)
end

@testitem "ParametricBuilder / wire patterns / screened constraints" tags = [:unit] begin
    required_area_mm2 = 35.0
    lay_diameter_mm = 60.0
    minimum_coverage = 85.0
    choices = make_screened(
        required_area_mm2,
        lay_diameter_mm;
        coverage_min_pct = minimum_coverage
    )

    for choice in values(choices)
        @test 1e6 * choice.total_area_m2 >= required_area_mm2
        @test minimum_coverage <= choice.coverage_pct <= 100.0
        @test choice.radius_m == (choice.lay_diameter_m + choice.wire_diameter_m) / 2
        separation = 2 * choice.radius_m * sinpi(1 / choice.wires)
        @test separation >= choice.wire_diameter_m
    end
    @test choices.min_wires.wires <= choices.min_diam.wires
    @test choices.min_diam.wire_diameter_m <= choices.min_wires.wire_diameter_m

    custom = make_screened(
        required_area_mm2,
        lay_diameter_mm;
        custom_diameters_mm = [1.2],
        nmin = 40,
        nmax = 40,
        max_overshoot_pct = Inf
    )
    @test any(choice -> startswith(choice.awg, "custom"), values(custom))
    @test any(choice -> choice.wire_diameter_m == 0.0012, values(custom))

    @test_throws AssertionError make_screened(0.0, lay_diameter_mm)
    @test_throws AssertionError make_screened(required_area_mm2, 0.0)
    @test_throws AssertionError make_screened(
        required_area_mm2,
        lay_diameter_mm;
        coverage_min_pct = 101.0
    )
end
