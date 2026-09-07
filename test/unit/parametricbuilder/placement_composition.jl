@testitem "ParametricBuilder / placement composition preserves geometry and connections" tags=[:unit] begin
    copper = Material(kind=:conductor, rho=1.72e-8)
    part = @terminal :phase begin
        core(copper; r=1e-3)
    end
    design = @cable "placed-wire" begin
        part
    end
    inner = at(0.03, -0.02; φ=0.2)
    outer = at(0.4, -1.0; φ=pi / 3)
    @test at(inner) === inner
    member = at(part, inner)
    local_placements = (member, inner)
    transformed = at(local_placements, outer)
    @test transformed isa Tuple
    @test transformed[1].item === part
    @test transformed[1].at == outer * inner
    @test transformed[2] == outer * inner
    @test at(collect(local_placements), outer) == collect(transformed)
    @test member.at == inner

    poses = (inner, at(-0.01, 0.02; φ=-0.4))
    tuple_space = Gridspace{Tuple}(pose -> (at(part, pose), pose), (Grid(poses),))
    shifted_tuples = at(tuple_space, outer)
    @test shifted_tuples isa Gridspace{Tuple}
    @test collect(shifted_tuples) == [at((at(part, pose), pose), outer) for pose in poses]

    connections = (phase=(1, 2, 3),)
    family = trefoil(design; spacing=Grid((0.02, 0.04)), connections)
    transformed_family = at(family, outer)
    @test transformed_family isa Gridspace{Vector}
    @test collect(transformed_family) == [at(formation, outer) for formation in family]
    offsets = at(Grid((0.1, 0.2)), -1.0)
    product = at(family, offsets)
    zipped = at(family, offsets; combine=:zip)
    @test length(product) == 4
    @test length(zipped) == 2
    @test collect(product) == [at(formation, pose) for pose in offsets for formation in family]
    @test collect(zipped) == [at(formation, pose) for (formation, pose) in zip(family, offsets)]
    for formation in product
        @test all(placed -> placed.design === design, formation)
        @test getproperty.(getproperty.(formation, :connections), :phase) == [1, 2, 3]
    end

    for formation in (trefoil, hflat, vflat)
        named = formation(design; spacing=0.02, connections)
        dictionary = formation(design; spacing=0.02, connections=Dict(:phase=>[1, 2, 3]))
        individual = formation(design; spacing=0.02,
            connections=((phase=1,), (phase=2,), (phase=3,)))
        @test getproperty.(dictionary, :pose) == getproperty.(named, :pose)
        @test [placed.connections[:phase] for placed in dictionary] == [1, 2, 3]
        @test individual == named
        for declaration in ((phase=2,), Dict(:phase=>2))
            common = formation(design; spacing=0.02, connections=declaration)
            @test [placed.connections[:phase] for placed in common] == [2, 2, 2]
        end
        for declaration in ((phase=(1, 2),), Dict(:phase=>[1, 2]), ((phase=1,), (phase=2,)))
            @test_throws DimensionMismatch formation(design; spacing=0.02, connections=declaration)
        end
    end
end

@testitem "ParametricBuilder / placed formations build systems and scalar problems consistently" tags=[:unit] begin
    copper = Material(kind=:conductor, rho=1.72e-8)
    design = @cable "placed-problem" begin
        @terminal :core begin
            core(copper; r=1e-3)
        end
    end
    soil = homogeneous(rho=100.0)
    family = at(hflat(design; spacing=Grid((0.02, 0.04)),
        connections=(core=(1, 2, 3),)), Pose2(0.0, -1.0))
    systems = build(LineCableSystem, family; environment=soil,
        system_id="placed-system", line_length=250.0)
    problems = LineParametersProblem(family; environment=soil,
        system_id="placed-system", line_length=250.0,
        earth_props=soil, frequencies=[50.0, 1000.0], temperature=40.0)
    @test systems isa Gridspace{LineCableSystem}
    @test problems isa Gridspace{LineParametersProblem}
    @test length(systems) == length(problems) == 2
    for (placed, system, problem) in zip(family, systems, problems)
        @test system.positions == getproperty.(placed, :pose)
        @test system.environment === soil
        @test system.line_length == 250.0
        @test problem.temperature == 40.0
        @test problem.frequencies == [50.0, 1000.0]
        @test problem.earth_props === soil
        # Grouped placement collections lower to the same physical ordering.
        grouped = (placed[1:2], placed[3:3])
        assembled = build(LineCableSystem, grouped; environment=soil,
            system_id="placed-system", line_length=250.0)
        for candidate in (problem.system, assembled), property in (
                :system_id, :line_length, :designs, :positions, :connections,
                :environment, :terminal_order, :terminal_map, :connection_order)
            @test getproperty(candidate, property) == getproperty(system, property)
        end
        direct = LineParametersProblem(assembled; earth_props=soil,
            frequencies=problem.frequencies, temperature=problem.temperature)
        actual = @inferred compute(problem)
        expected = @inferred compute(direct)
        @test actual.Z.values == expected.Z.values
        @test actual.Y.values == expected.Y.values
    end
    soils = homogeneous(rho=Grid((10.0, 100.0)))
    product = LineParametersProblem(systems, soils; frequencies=[50.0])
    zipped = LineParametersProblem(systems, soils; frequencies=[50.0], combine=:zip)
    @test length(product) == 4
    @test length(zipped) == 2
    @test [problem.earth_props.layers[2].rho for problem in zipped] == [10.0, 100.0]
    @test_throws ArgumentError build(LineCableSystem, ())
    @test_throws ArgumentError build(LineCableSystem, ((),))
    @test_throws ArgumentError build(LineCableSystem, ((design=design, pose=Pose2(0.0, -1.0)),))
end
