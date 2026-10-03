@testitem "ParametricBuilder / wire patterns / stranded estimates" tags = [:unit, :parametric] begin
    const Wires = LineCableModels.ParametricBuilder.WirePatterns

    for awg in (-3, 0, 12, 40)
        diameter = Wires.awg_to_d_mm(awg)
        area = Wires.awg_to_area_mm2(awg)
        @test Wires.d_mm_to_awg(diameter) ≈ awg atol = 64eps(Float64)
        @test Wires.area_mm2_to_awg(area) ≈ awg atol = 64eps(Float64)
    end
    @test Wires.awg_label(-3) == "0000 (4/0)"
    @test Wires.awg_label(12) == "12"

    target_area = 1_000.0
    estimate = estimate_stranding(target_area)
    @test estimate isa WireEstimate{Float64}
    @test estimate.feasible
    @test estimate.status === :feasible
    @test isempty(estimate.reasons)
    for choice in estimate
        @test choice.wires == 1 + 3 * choice.layers * (choice.layers - 1)
        @test choice.total_area > 0
        @test choice.wire_diameter > 0
    end
    match = estimate[:closest_area]
    layers = estimate[:fewest_layers]
    diameter = estimate[:smallest_diameter]
    @test abs(1e6 * match.total_area - target_area) <=
          abs(1e6 * diameter.total_area - target_area)
    @test layers.layers <= diameter.layers
    @test estimate[:fewest_wires].wires == minimum(p.wires for p in estimate)
    @test estimate[Val(:closest_area)] === estimate[:closest_area]

    estimate32 = estimate_stranding(Float32(95))
    @test estimate32 isa WireEstimate{Float32}
    @test typeof(first(estimate32).total_area) === Float32

    infeasible = estimate_stranding(1.0e12; awg_min = 40, awg_max = 40)
    @test !infeasible.feasible
    @test infeasible.status === :infeasible
    @test !isempty(infeasible.reasons)
    @test infeasible[:closest_area] isa Wires.HexaPattern

    @test_throws DomainError estimate_stranding(0.0)
    @test_throws ArgumentError estimate_stranding(10.0; awg_min = 10, awg_max = 9)
    @test_throws ArgumentError estimate[:unknown]
end

@testitem "ParametricBuilder / wire patterns / screened estimates" tags = [:unit, :parametric] begin
    const Wires = LineCableModels.ParametricBuilder.WirePatterns

    target_area = 35.0
    lay_diameter = 60.0
    minimum_coverage = 85.0
    estimate = estimate_screen(
        target_area,
        lay_diameter;
        coverage_min = minimum_coverage
    )

    @test estimate isa WireEstimate{Float64}
    @test estimate.feasible
    for choice in estimate
        @test 1e6 * choice.total_area >= target_area
        @test minimum_coverage <= choice.coverage <= 100.0
        @test choice.radius ==
              (choice.lay_diameter + choice.wire_diameter) / 2
        separation = 2 * choice.radius * sinpi(1 / choice.wires)
        @test separation >= choice.wire_diameter
    end
    @test estimate[:fewest_wires].wires <= estimate[:smallest_diameter].wires
    @test estimate[:smallest_diameter].wire_diameter <= estimate[:fewest_wires].wire_diameter

    custom = estimate_screen(
        target_area,
        lay_diameter;
        wire_diameters = [1.2],
        awg_min = 40,
        awg_max = 40,
        max_area_overshoot = Inf
    )
    @test any(choice -> startswith(choice.awg, "custom"), custom)
    @test any(choice -> choice.wire_diameter == 0.0012, custom)

    screen = first(estimate)
    @test Wires.maxfill(
        Wires.ScreenPattern,
        screen.radius,
        screen.wire_diameter / 2
    ) >= screen.wires

    infeasible = estimate_screen(
        target_area,
        lay_diameter;
        coverage_min = 99,
        coverage_max = 99
    )
    @test !infeasible.feasible
    @test !isempty(infeasible.reasons)
    @test infeasible[:closest_area] isa Wires.ScreenPattern

    @test_throws DomainError estimate_screen(0.0, lay_diameter)
    @test_throws DomainError estimate_screen(target_area, 0.0)
    @test_throws DomainError estimate_screen(
        target_area,
        lay_diameter;
        coverage_min = 101.0
    )
end
