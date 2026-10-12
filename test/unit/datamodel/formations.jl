@testitem "DataModel / formations / a spacing below the diameter lays the cables touching" tags=[:unit, :datamodel] begin
    using LineCableModels
    copper = Material(kind = :conductor, rho = 1.7e-8)
    design = build(CableDesign, "formation", Group(:core, Region(:metal, Disk(0.01), copper)))
    connections = (core = (1, 2, 3),)
    distance(a, b) = hypot(a.x - b.x, a.y - b.y)
    # The side of the trefoil, or the adjacent spacing of a flat formation.
    side(poses) = minimum(distance(poses[i], poses[j]) for i in 1:3 for j in 1:(i - 1))
    for formation in (trefoil, hflat, vflat)
        for spacing in (0.02, 0.015, 0.001)
            records = formation(design; center = at(0.0, -1.0), spacing, connections)
            poses = getproperty.(records, :pose)
            @test side(poses) - 2outer_radius(design) ≈ 1e-6 atol = 1e-15
            # The build checks the formation and moves nothing.
            @test (@test_logs build(LineCableSystem, records)).positions == poses
        end
        wide = getproperty.(formation(design; spacing = 0.05, connections), :pose)
        @test side(wide) ≈ 0.05 rtol = 1e-15
    end
end

@testitem "DataModel / formations / each realization is laid with its own radius" tags=[:unit, :measurements] begin
    using LineCableModels
    using Measurements: measurement, derivative
    copper = Material(kind = :conductor, rho = 1.7e-8)
    cable(radius) = build(CableDesign, "formation",
        Group(:core, Region(:metal, Disk(radius), copper)))
    connections = (core = (1, 2, 3),)
    distance(a, b) = hypot(a.x - b.x, a.y - b.y)
    gaps(system) = [distance(system.positions[i], system.positions[j]) -
                    outer_radius(system.designs[i]) - outer_radius(system.designs[j])
                    for i in 1:3 for j in 1:(i - 1)]

    # First-order propagation: the spacing of a touching formation depends on the radius.
    radius = measurement(0.01, 1e-4)
    for formation in (trefoil, hflat, vflat)
        touching = getproperty.(formation(cable(radius); spacing = 0.02, connections), :pose)
        @test nominal(distance(touching[1], touching[2])) - 0.02 ≈ 1e-6 atol = 1e-15
        @test derivative(distance(touching[1], touching[2]), radius) ≈ 2
        wide = getproperty.(formation(cable(radius); spacing = 0.05, connections), :pose)
        @test derivative(distance(wide[1], wide[2]), radius) == 0
    end

    # Monte Carlo: a larger drawn cable stays touching, so no draw overlaps.
    built = LineCableSystem[]
    space = Gridspace{LineParametersProblem}(
        r -> begin
            system = build(LineCableSystem, trefoil(cable(r); center = at(0.0, -1.0),
                spacing = 0.02, connections))
            push!(built, system)
            LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
                frequencies = [50.0])
        end,
        (Grid(0.01, AbsoluteError(1e-4)),))
    inner = Formulation(options = (kron_reduction = false,
        reduce_bundle = false, ideal_transposition = false))
    sampled = compute(ParametricProblem(space), MonteCarlo(inner; trials = 8, seed = 5))
    @test sampled.trial_counts == [8]
    @test length(built) == 8 + 1 # every draw, and the nominal record
    smallest = [nominal(minimum(gaps(system))) for system in built]
    @test all(>=(1e-6 - 1e-15), smallest)
    @test count(gap -> isapprox(gap, 1e-6; atol = 1e-12), smallest[1:8]) > 0
    @test count(gap -> gap > 1e-6 + 1e-12, smallest[1:8]) > 0
end
