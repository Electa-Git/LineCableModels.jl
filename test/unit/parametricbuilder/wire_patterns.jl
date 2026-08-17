@testitem "ParametricBuilder / wire patterns / stranded estimates" tags = [:unit] begin
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
    estimate = make_stranded(target_mm2)
    @test estimate isa WireEstimate{Float64}
    @test estimate.feasible
    @test estimate.status === :feasible
    @test isempty(estimate.reasons)
    for choice in estimate
        @test choice.wires == 1 + 3 * choice.layers * (choice.layers - 1)
        @test choice.total_area_m2 > 0
        @test choice.wire_diameter_m > 0
    end
    match = estimate[Val(:match)]
    layers = estimate[Val(:layers)]
    diameter = estimate[Val(:diameter)]
    @test abs(1e6 * match.total_area_m2 - target_mm2) <=
          abs(1e6 * diameter.total_area_m2 - target_mm2)
    @test layers.layers <= diameter.layers
    @test estimate[Val(:wires)].wires == minimum(p.wires for p in estimate)

    estimate32 = make_stranded(Float32(95))
    @test estimate32 isa WireEstimate{Float32}
    @test typeof(first(estimate32).total_area_m2) === Float32

    infeasible = make_stranded(1.0e12; nmin = 40, nmax = 40)
    @test !infeasible.feasible
    @test infeasible.status === :infeasible
    @test !isempty(infeasible.reasons)
    @test infeasible[Val(:match)] isa Wires.HexaPattern

    @test_throws DomainError make_stranded(0.0)
    @test_throws ArgumentError make_stranded(10.0; nmin = 10, nmax = 9)
    @test_throws ArgumentError estimate[Val(:unknown)]
end

@testitem "ParametricBuilder / wire patterns / screened estimates" tags = [:unit] begin
    const Wires = LineCableModels.ParametricBuilder.WirePatterns

    required_area_mm2 = 35.0
    lay_diameter_mm = 60.0
    minimum_coverage = 85.0
    estimate = make_screened(
        required_area_mm2,
        lay_diameter_mm;
        coverage_min_pct = minimum_coverage
    )

    @test estimate isa WireEstimate{Float64}
    @test estimate.feasible
    for choice in estimate
        @test 1e6 * choice.total_area_m2 >= required_area_mm2
        @test minimum_coverage <= choice.coverage_pct <= 100.0
        @test choice.radius_m ==
              (choice.lay_diameter_m + choice.wire_diameter_m) / 2
        separation = 2 * choice.radius_m * sinpi(1 / choice.wires)
        @test separation >= choice.wire_diameter_m
    end
    @test estimate[Val(:wires)].wires <= estimate[Val(:diameter)].wires
    @test estimate[Val(:diameter)].wire_diameter_m <=
          estimate[Val(:wires)].wire_diameter_m

    custom = make_screened(
        required_area_mm2,
        lay_diameter_mm;
        custom_diameters_mm = [1.2],
        nmin = 40,
        nmax = 40,
        max_overshoot_pct = Inf
    )
    @test any(choice -> startswith(choice.awg, "custom"), custom)
    @test any(choice -> choice.wire_diameter_m == 0.0012, custom)

    screen = first(estimate)
    @test Wires.maxfill(
        Wires.ScreenPattern,
        screen.radius_m,
        screen.wire_diameter_m / 2
    ) >= screen.wires

    infeasible = make_screened(
        required_area_mm2,
        lay_diameter_mm;
        coverage_min_pct = 99,
        coverage_max_pct = 99
    )
    @test !infeasible.feasible
    @test !isempty(infeasible.reasons)
    @test infeasible[Val(:match)] isa Wires.ScreenPattern

    @test_throws DomainError make_screened(0.0, lay_diameter_mm)
    @test_throws DomainError make_screened(required_area_mm2, 0.0)
    @test_throws DomainError make_screened(
        required_area_mm2,
        lay_diameter_mm;
        coverage_min_pct = 101.0
    )
end
