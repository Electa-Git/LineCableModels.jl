@testmodule ClearanceFixtures begin
    using LineCableModels

    function cable(radius = 0.01; offset = Pose2(0, 0))
        copper = Material(kind = :conductor, rho = 1.7e-8)
        return build(CableDesign, "clearance-test",
            Group(:core, Region(:metal, Disk(radius), copper); at = offset))
    end

    function pair(radius, separation; y = -1.0)
        design = cable(radius)
        return build(LineCableSystem, [design, design],
            [Pose2(-separation / 2, y), Pose2(separation / 2, y)];
            connections = [(core = 1,), (core = 2,)])
    end

    function gaps(system)
        shapes = [resolve(pose, boundary(design.geometry))
                  for (pose, design) in zip(system.positions, system.designs)]
        [hypot(centroid(shapes[i])[1] - centroid(shapes[j])[1],
            centroid(shapes[i])[2] - centroid(shapes[j])[2]) -
            shapes[i].r - shapes[j].r
         for i in eachindex(system.positions) for j in 1:(i - 1)]
    end
end

@testitem "DataModel / exterior clearance / touching, overlap and reconstruction" tags=[:unit] setup=[ClearanceFixtures] begin
    using LineCableModels
    import LineCableModels.ImportExport as IE
    fixture = ClearanceFixtures
    design = fixture.cable()
    declared = trefoil(design; center = at(0.0, -1.0), spacing = 0.02,
        connections = (core = (1, 2, 3),))
    system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem, declared)
    @test all(>=(1.0e-6), fixture.gaps(system))
    @test maximum(fixture.gaps(system)) - minimum(fixture.gaps(system)) < 1e-14
    @test sum(p -> p.x, system.positions) ≈ sum(p -> p.x, system.declared_positions) atol=1e-15
    @test sum(p -> p.y, system.positions) ≈ sum(p -> p.y, system.declared_positions)
    @test getproperty.(declared, :pose) == system.declared_positions
    @test all(==(design), system.designs)
    @test validate(system) === system
    rebuilt = @test_logs build(LineCableSystem, system.designs, system.positions;
        connections = system.connections)
    @test rebuilt.positions == system.positions
    decoded = @test_logs IE.deserialize_value(IE.serialize_value(system))
    @test decoded.positions == system.positions
    @test decoded.declared_positions == system.declared_positions
    @test decoded.clearances == system.clearances
    legacy = IE.serialize_value(system)
    delete!(legacy, "declared_positions")
    delete!(legacy, "clearances")
    @test IE.deserialize_value(legacy).positions == system.positions

    reversed = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem, reverse(declared))
    @test getproperty.(reverse(reversed.positions), :x) ≈ getproperty.(system.positions, :x) atol=1e-14
    @test getproperty.(reverse(reversed.positions), :y) ≈ getproperty.(system.positions, :y) atol=1e-14
    @test_logs fixture.pair(0.01, 0.03)
    @test_logs (:warn, r"Cable placements adjusted") fixture.pair(0.01, prevfloat(0.02))
    @test_throws DomainError fixture.pair(0.01, 0.02 - 1e-10)
    @test_throws DomainError fixture.pair(0.01, 0.0)

    # Several neighbours must be resolved together, not by a one-pass pair push.
    line = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
        fill(design, 5), [Pose2(0.02i, -1) for i in 1:5])
    @test minimum(fixture.gaps(line)) >= 1e-6
    line.positions[2] = line.positions[1]
    @test_throws DomainError validate(line)

    # A tangent interface is resolved at problem construction, where earth is known.
    neutral = build(LineCableSystem, design, Pose2(0, -0.01))
    problem = @test_logs (:warn, r"Cable placements adjusted") LineParametersProblem(
        neutral; earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    @test -only(problem.system.positions).y - outer_radius(design) >= 1e-6
    @test only(neutral.positions).y == -0.01
    @test_throws DomainError LineParametersProblem(
        build(LineCableSystem, design, Pose2(0, -0.009));
        earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    small_precision = @test_logs (:warn, r"Cable placements adjusted") fixture.pair(
        0.01f0, 0.02f0; y = -1.0f0)
    @test only(fixture.gaps(small_precision)) >= 1e-6
end

@testitem "DataModel / exterior clearance / propagated radii and correlated at coordinates" tags=[:extension] setup=[ClearanceFixtures] begin
    using LineCableModels
    using Measurements: measurement, derivative
    fixture = ClearanceFixtures
    radius = measurement(0.01, 1e-4)
    shared_shift = measurement(0.3, 0.004)
    design = fixture.cable(radius)
    declarations = [
        (@at design (shared_shift - radius, -1.0) core=1),
        (@at design (shared_shift + radius, -1.0) core=2)
    ]
    system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem, declarations)
    @test nominal(system.clearances[1, 2]) ≈ 1e-6 + uncertainty(radius)
    @test nominal(only(fixture.gaps(system))) >= nominal(system.clearances[1, 2])
    @test uncertainty(only(fixture.gaps(system))) < 1e-14
    @test derivative(sum(p -> p.x, system.positions) / 2, shared_shift) ≈ 1
    @test derivative(system.positions[2].x - system.positions[1].x, radius) ≈ 2

    independent = measurement(0.02, 0.003)
    positional = @test_logs (:warn, r"Cable placements adjusted") fixture.pair(0.01, independent)
    @test nominal(positional.clearances[1, 2]) ≈ 1e-6 + 0.003
    @test uncertainty(only(fixture.gaps(positional))) < 1e-14

    left = fixture.cable(measurement(0.01, 1e-4))
    right = fixture.cable(measurement(0.01, 2e-4))
    radii = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
        [left, right], [Pose2(-0.01, -1), Pose2(0.01, -1)])
    @test nominal(radii.clearances[1, 2]) ≈ 1e-6 + hypot(1e-4, 2e-4)

    # Local cross-section poses contribute to the exterior envelope as well.
    local_offset = measurement(0.002, 5e-4)
    eccentric = fixture.cable(0.01; offset = Pose2(local_offset, 0))
    @test uncertainty(outer_radius(eccentric)) ≈ 5e-4
    local_system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
        [eccentric, eccentric], [Pose2(-0.01, -1), Pose2(0.01, -1)])
    @test nominal(local_system.clearances[1, 2]) ≈ 1e-6 + 0.0005
    @test nominal(only(fixture.gaps(local_system))) >= nominal(local_system.clearances[1, 2])
    zero_offset = fixture.cable(0.01; offset = Pose2(measurement(0.0, 5e-4), 0))
    zero_system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
        [zero_offset, zero_offset], [Pose2(-0.01, -1), Pose2(0.01, -1)])
    @test all(isfinite, nominal.(zero_system.clearances))
    @test nominal(zero_system.clearances[1, 2]) >= 1e-6 + 0.0005
    @test isfinite(uncertainty(only(fixture.gaps(zero_system))))
    compensated = fixture.cable(radius;
        offset = Pose2(0.002 - (radius - nominal(radius)), 0))
    @test uncertainty(outer_radius(compensated)) < 1e-14
    compensated_system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
        [compensated, compensated], [Pose2(-radius, -1), Pose2(radius, -1)])
    @test nominal(compensated_system.clearances[1, 2]) ≈ 1e-6 + uncertainty(radius)
end

@testitem "UQ / exterior clearance / retained budgets, extreme draws and bounded warnings" tags=[:integration] setup=[ClearanceFixtures] begin
    using LineCableModels
    using Measurements
    using Random
    fixture = ClearanceFixtures
    DM = LineCableModels.DataModel
    space = Gridspace{LineCableSystem}(
        (radius, spacing) -> fixture.pair(radius, spacing),
        (Grid(0.01, AbsoluteError(1e-4)), Grid(0.021, AbsoluteError(0.004))))
    point = only(LineCableModels.points(space))
    context = DM.prepare_clearance(point)
    @test context.records[1].clearances[1, 2] > 0.004
    @test DM._CLEARANCE_CONTEXT[] === nothing
    for (radius, spacing) in ((0.0108, 0.015), (0.0112, 0.006), (0.009, 0.045), (0.01, 0.0))
        result = @test_logs DM.with_clearance(context) do
            LineCableModels.realize(point, (radius, spacing))
        end
        @test only(fixture.gaps(result)) >= context.records[1].clearances[1, 2]
        @test result.clearances == context.records[1].clearances
    end
    @test context.adjustments[] == 3
    @test DM._CLEARANCE_CONTEXT[] === nothing
    @test_throws DomainError DM.with_clearance(context) do
        LineCableModels.realize(point, (-0.01, 0.02))
    end
    @test DM._CLEARANCE_CONTEXT[] === nothing
    @test_logs (:warn, r"Sampled cable placements adjusted") rand(MersenneTwister(8), space)
    zipped = Gridspace{LineCableSystem}(space.build, space.grids; combine = :zip)
    @test_logs (:warn, r"Sampled cable placements adjusted") rand(MersenneTwister(8), zipped)

    # Even a draw that swaps air/earth sides must remain pairwise feasible.
    design = fixture.cable()
    reference = [Pose2(0, 0.1), Pose2(0, -0.1)]
    required = [0.001 0.1; 0.1 0.001]
    resolved, _, _ = DM.clearance_geometry([design, design],
        [Pose2(0, -0.1), Pose2(0, 0.1)]; reference, required,
        interface = true, sampling = true)
    @test resolved[1].y > 0 && resolved[2].y < 0
    @test abs(resolved[1].y - resolved[2].y) - 0.02 >= 0.1

    # Retain the side of the composed local @at, not the sampled local offset.
    # An unrelated second system with the same default ID must not replace it.
    local_space = Gridspace{LineParametersProblem}(offset -> begin
        local_design = fixture.cable(0.001; offset = Pose2(0, offset))
        buried = build(LineCableSystem, local_design, Pose2(0, -0.01))
        build(LineCableSystem, local_design, Pose2(0, 0.3))
        LineParametersProblem(buried; earth_props = homogeneous(rho = 100.0),
            frequencies = [50.0])
    end, (Grid(0.0, AbsoluteError(0.05)),))
    local_point = only(LineCableModels.points(local_space))
    local_context = DM.prepare_clearance(local_point)
    local_problem = @test_logs DM.with_clearance(local_context) do
        LineCableModels.realize(local_point, (0.15,))
    end
    exterior = resolve(only(local_problem.system.positions),
        boundary(only(local_problem.system.designs).geometry))
    @test -centroid(exterior)[2] - exterior.r >= 0.050001
    @test length(local_context.references) <= 3

    # Exercise the actual MC loop with owned engine results and forced inward draws.
    constructed = LineCableSystem[]
    problems = Gridspace{LineParametersProblem}(spacing -> begin
        system = fixture.pair(0.01, spacing)
        push!(constructed, system)
        LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
            frequencies = [50.0])
    end, (Grid(0.021, AbsoluteError(0.002)),))
    inner = Formulation(options = (kron_reduction = false,
        reduce_bundle = false, ideal_transposition = false))
    formulation = MonteCarlo(inner; trials = 3, seed = 8,
        distribution = (_rng, mean, sigma) -> mean - 4sigma,
        options = (retain_details = true,))
    sampled = @test_logs (:warn, r"Sampled cable placements adjusted") compute(
        ParametricProblem(problems), formulation)
    @test sampled.trial_counts == [3]
    @test isempty(only(sampled.details.failures))
    @test only(sampled.details.clearance).adjustments == 3
    @test length(constructed) == 4 # one uncertain preparation, three draws
    @test all(system -> nominal(only(fixture.gaps(system))) >= 0.002001, constructed)
    @test DM._CLEARANCE_CONTEXT[] === nothing
end
