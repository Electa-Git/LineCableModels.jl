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

@testitem "DataModel / exterior clearance / the 1 μm floor, overlaps and reconstruction" tags=[:unit, :importexport] setup=[ClearanceFixtures] begin
    using LineCableModels
    import LineCableModels.ImportExport as IE
    fixture=ClearanceFixtures
    design=fixture.cable()

    # A pair gap from zero, within roundoff, up to 1 μm is raised to 1 μm, silently and
    # symmetrically about the centroid.
    for separation in (0.02, prevfloat(0.02), 0.02 + 5e-7)
        system=@test_logs fixture.pair(0.01, separation)
        @test only(fixture.gaps(system)) ≈ 1e-6 atol=1e-14
        @test sum(p -> p.x, system.positions) ≈ 0 atol=1e-15
        @test all(p -> p.y == -1.0, system.positions)
    end
    # A gap of 1 μm or more stays as declared.
    @test fixture.pair(0.01, 0.03).positions == [Pose2(-0.015, -1.0), Pose2(0.015, -1.0)]
    overlap=try
        fixture.pair(0.01, 0.02 - 1e-10)
    catch error
        error
    end
    @test overlap isa OverlapError
    @test overlap.cables == (1, 2)
    @test overlap.gap ≈ -1e-10 rtol=1e-4
    @test occursin("overlap", sprint(showerror, overlap))
    @test_throws OverlapError fixture.pair(0.01, 0.0)
    small_precision=@test_logs fixture.pair(0.01f0, 0.02f0; y = -1.0f0)
    @test only(fixture.gaps(small_precision)) >= 1e-6

    # Several neighbors are resolved together, by one expansion.
    line=@test_logs build(LineCableSystem, fill(design, 5), [Pose2(0.02i, -1) for i in 1:5])
    @test minimum(fixture.gaps(line)) >= 1e-6
    @test validate(line) === line
    # A rebuild from the built positions keeps them, and so does a saved copy.
    rebuilt=@test_logs build(LineCableSystem, line.designs, line.positions;
        connections = line.connections)
    @test rebuilt.positions == line.positions
    encoded=IE.serialize_value(line)
    @test !haskey(encoded, "input_positions") && !haskey(encoded, "clearances")
    @test IE.deserialize_value(encoded).positions == line.positions
    # Called on its own, `validate` checks a system again and catches a changed vector.
    line.positions[2]=line.positions[1]
    @test_throws OverlapError validate(line)

    # The air-earth interface has the same rule. A line problem supplies the earth.
    earth=homogeneous(rho = 100.0)
    tangent=build(LineCableSystem, design, Pose2(0, -0.01))
    problem=@test_logs LineParametersProblem(tangent; earth_props = earth, frequencies = [50.0])
    @test -only(problem.system.positions).y - outer_radius(design) ≈ 1e-6 atol=1e-15
    @test only(tangent.positions).y == -0.01
    for crossing in (Pose2(0, -0.009), Pose2(0, 0))
        @test_throws OverlapError LineParametersProblem(build(LineCableSystem, design, crossing);
            earth_props = earth, frequencies = [50.0])
    end
    @test_throws OverlapError build(LineCableSystem, design, Pose2(0, -0.009); environment = earth)

    # The system checks its pair gaps, and the problem checks the interface gaps.
    placed=LineParametersProblem(fixture.pair(0.01, 0.03); earth_props = earth,
        frequencies = [50.0])
    placed.system.positions[2]=Pose2(-0.005, -1.0)
    @test validate(placed) === placed
    @test_throws OverlapError validate(placed.system)
    placed.system.positions[2]=Pose2(0.015, -0.005)
    @test_throws OverlapError validate(placed)
end

@testitem "DataModel / exterior clearance / uncertainty keeps the nominal geometry" tags=[:extension, :measurements] setup=[ClearanceFixtures] begin
    using LineCableModels
    using Measurements: measurement, derivative
    fixture=ClearanceFixtures
    radius=measurement(0.01, 1e-4)
    shared_shift=measurement(0.3, 0.004)
    design=fixture.cable(radius)

    # A nominal zero gap is raised to 1 μm. The positions keep their correlation with the
    # radius and the shared shift.
    touching=@test_logs build(LineCableSystem, [
        (@at design (shared_shift - radius, -1.0) core=1),
        (@at design (shared_shift + radius, -1.0) core=2)])
    @test nominal(only(fixture.gaps(touching))) ≈ 1e-6 atol=1e-14
    @test uncertainty(only(fixture.gaps(touching))) < 1e-14
    @test derivative(sum(p -> p.x, touching.positions) / 2, shared_shift) ≈ 1
    @test derivative(touching.positions[2].x - touching.positions[1].x, radius) ≈ 2

    # A declared gap of 1 μm or more stays as declared, whatever its uncertainty.
    left=fixture.cable(measurement(0.01, 1e-4))
    right=fixture.cable(measurement(0.01, 2e-4))
    for system in (fixture.pair(radius, 0.0205), fixture.pair(0.01, measurement(0.0205, 0.003)),
            build(LineCableSystem, [left, right], [Pose2(-0.01025, -1), Pose2(0.01025, -1)]))
        @test nominal.(getproperty.(system.positions, :x)) == [-0.01025, 0.01025]
    end
    independent=fixture.pair(0.01, measurement(0.0205, 0.003))
    @test uncertainty(independent.positions[2].x) ≈ 0.0015

    # A local cross-section offset moves the exterior circle, and the gap follows it.
    eccentric=fixture.cable(0.01; offset = Pose2(measurement(0.002, 5e-4), 0))
    offset=@test_logs build(LineCableSystem,
        [eccentric, eccentric], [Pose2(-0.01, -1), Pose2(0.01, -1)])
    @test nominal(only(fixture.gaps(offset))) ≈ 1e-6 atol=1e-14
    @test_throws OverlapError fixture.pair(radius, 0.019)
end

@testitem "UQ / exterior clearance / Monte Carlo discards overlapping draws" tags=[:integration, :measurements] setup=[ClearanceFixtures] begin
    using LineCableModels
    using Measurements
    fixture=ClearanceFixtures
    earth=homogeneous(rho = 100.0)
    drawn=Float64[]
    space=Gridspace{LineParametersProblem}(
        spacing->begin
            push!(drawn, spacing)
            LineParametersProblem(fixture.pair(0.01, spacing); earth_props = earth,
                frequencies = [50.0])
        end,
        (Grid(0.021, AbsoluteError(0.002)),))
    inner=Formulation(options = (kron_reduction = false,
        reduce_bundle = false, ideal_transposition = false))
    formulation=MonteCarlo(inner; trials = 8, seed = 8, return_samples = true,
        retain_details = true)

    # A draw whose cables overlap is discarded, as if it was never drawn, and drawn again.
    sampled=@test_logs compute(ParametricProblem(space), formulation)
    discarded=count(<(0.02), drawn)
    @test discarded > 0
    @test length(drawn) == 8 + discarded + 1 # accepted, discarded, the nominal record
    @test sampled.trial_counts == [8]
    @test keys(sampled.details.data) == (:trials,)
    accepted=filter(>=(0.02), drawn[1:(end - 1)])
    replay=compute(ParametricProblem(space), formulation)
    @test samples(replay) == samples(sampled)
    @test filter(>=(0.02), drawn[(8 + discarded + 2):(end - 1)]) == accepted

    # `max_failures` overlapping draws in a row stop a study that cannot draw a feasible
    # case. The distribution of the point has no feasible draw, a `DomainError` that a
    # Combinatorial study records.
    empty!(drawn)
    infeasible=MonteCarlo(inner; trials = 2, seed = 8, max_failures = 5,
        distribution = (_rng, mean, sigma)->mean-4sigma)
    failure=try
        compute(ParametricProblem(space), infeasible)
    catch error
        error
    end
    @test failure isa DomainError
    @test only(failure.val).standard_deviation == 0.002
    @test occursin("max_failures", sprint(showerror, failure))
    @test occursin("OverlapError", sprint(showerror, failure))
    @test LineCableModels.on_error(Combinatorial(inner), failure) === failure
    @test length(drawn) == 5

    # A grid point is a requested case, and the study records an overlapping one.
    requested=Gridspace{LineParametersProblem}(
        spacing->LineParametersProblem(fixture.pair(0.01, spacing); earth_props = earth,
            frequencies = [50.0]),
        (Grid((0.03, 0.015)),))
    recorded=compute(ParametricProblem(requested), Combinatorial(inner))
    @test length(recorded) == 1
    @test last(only(recorded.failures)) isa OverlapError
end

@testitem "DataModel / exterior clearance / a study records infeasible points" tags=[:integration, :measurements] setup=[ClearanceFixtures, SavedResults] begin
    using LineCableModels
    using Measurements
    fixture=ClearanceFixtures
    inner=Formulation(options = (kron_reduction = false,
        reduce_bundle = false, ideal_transposition = false))
    # The cables of point 2 overlap, and the earth of point 3 has no permittivity: both
    # builds throw. The earth of point 4 has a layer over the half space. The space builds two
    # problem types. `shift` is the uncertainty of the spacing.
    earth(eps_r, layered)=layered ? build(EarthModel,
        (layer(; rho = 100.0, eps_r, thickness = 0.5), layer(rho = 500.0, eps_r = 20.0))) :
        homogeneous(; rho = 100.0, eps_r)
    problem=(spacing, eps_r, layered, shift)->LineParametersProblem(
        fixture.pair(0.01, spacing + shift); earth_props = earth(eps_r, layered),
        frequencies = [50.0, 1e3])
    space(shift)=Gridspace{LineParametersProblem}(problem,
        (Grid((0.03, 0.015, 0.03, 0.04)), Grid((10.0, 10.0, 0.0, 10.0)),
            Grid((false, false, false, true)), shift);
        combine = :zip)
    formulations=Grid((inner, Formulation()))
    result=compute(ParametricProblem(space(Grid((0.0,)))), Combinatorial(formulations))
    gridpoints(value)=[(gridpoint.problem_index, gridpoint.formulation_index)
                       for (gridpoint, _) in value.failures]
    @test gridpoints(result) == [(2, 1), (3, 1), (2, 2), (3, 2)]
    @test [typeof(exception) for (_, exception) in result.failures] ==
          [OverlapError, DomainError, OverlapError, DomainError]
    # The values equal a run without the failed points.
    alone=compute(Gridspace{LineParametersProblem}(problem,
            (Grid((0.03, 0.04)), Grid((10.0, 10.0)), Grid((false, true)), Grid((0.0,)));
            combine = :zip),
        formulations)
    @test isempty(alone.failures)
    @test length(result) == length(alone) == 4
    for (computed, reference) in zip(result, alone)
        @test Z(computed) == Z(reference)
        @test Y(computed) == Y(reference)
    end
    @test Z(result[4, 2]) == Z(alone[2, 2])
    @test_throws ErrorException result[2, 1]
    # A comparison compares the computed points and keeps the failures.
    compared=LineCableModels.Engine.compare(first(alone), result, Z)
    @test length(compared) == 4
    @test compared.failures === result.failures
    # Linear propagation over the same space records the same failures.
    linear=compute(ParametricProblem(space(Grid(0.0, AbsoluteError(1e-4)))), LinearError(inner))
    @test gridpoints(linear) == [(2, 1), (3, 1)]
    @test [typeof(exception) for (_, exception) in linear.failures] == [OverlapError, DomainError]
    @test length(linear) == 2
    @test sprint(show, linear) == "LinearErrorResult(2 values, 2 failed)"
    # A saved result keeps its failures.
    for saved in (result, linear)
        restored=SavedResults.round_trip(saved)
        @test SavedResults.same(restored, saved)
        @test gridpoints(restored) == gridpoints(saved)
    end
end
