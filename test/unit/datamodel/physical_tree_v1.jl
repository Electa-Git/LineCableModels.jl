@testitem "DataModel / v1 physical tree / Region and Stack" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    const PB=LineCableModels.ParametricBuilder
    copper_props=LineCableModels.Material(kind = :conductor, rho = 1.7241e-8)
    insulator_props=LineCableModels.Material(
        kind = :insulator, rho = 1.0e14, eps_r = 2.3
    )

    core=DM.Region(:core, DM.Disk(0.01), copper_props)
    insulation=DM.Region(:insulation, DM.Shell(0.005), insulator_props)
    @test core isa DM.AbstractCablePart
    @test core.tag === :core
    @test core.material === copper_props
    @test (@inferred DM.resolve(DM.EmptyBoundary(), core)).regions[1].shape ==
          DM.DiskShape(0.01)
    @test_throws ArgumentError DM.Region(Symbol(""), DM.Disk(1), copper_props)
    @test_throws ArgumentError DM.Region(:bad, DM.Disk(1), :copper)

    stack=PB.Stack(core, insulation)
    @test stack.items == DM.AbstractCablePart[core, insulation]
    resolved=DM.resolve(DM.EmptyBoundary(), stack)
    @test length(resolved.regions) == 2
    @test resolved.regions[1].shape == DM.DiskShape(0.01)
    @test resolved.regions[2].shape == DM.AnnulusShape(0.01, 0.015)
    @test DM.boundary(resolved) == DM.DiskShape(0.015)
    @test_throws ArgumentError DM.Stack(DM.AbstractCablePart[])
    @test_throws MethodError DM.resolve(DM.EmptyBoundary(), DM.Stack(insulation))

    nested=PB.Stack(PB.Stack(core), insulation)
    nested_resolved=DM.resolve(DM.EmptyBoundary(), nested)
    @test DM.boundary(nested_resolved) == DM.DiskShape(0.015)

    radii=PB.Grid((0.01, 0.02))
    cores=PB.Region(:core, PB.Disk(radii), copper_props)
    spaces=PB.Stack(cores, insulation)
    @test spaces isa PB.Gridspace{DM.Stack}
    @test length(spaces) == 2
    @test all(space -> space isa DM.Stack, spaces)
    @test DM.r_ex(DM.boundary(DM.resolve(DM.EmptyBoundary(), first(spaces)))) == 0.015
end

@testitem "DataModel / v1 physical tree / Enclosure and class conveniences" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    copper=LineCableModels.Material(kind = :conductor, rho = 1.7241e-8)
    oil=LineCableModels.Material(kind = :insulator, rho = 1.0e12, eps_r = 2.2)
    semicon=LineCableModels.Material(kind = :semicon, rho = 1000.0, eps_r = 1000.0)

    core_region=LineCableModels.Conductor.Solid(:core_metal, copper; r = 1.0)
    core=LineCableModels.Group(:core, core_region)
    wall=LineCableModels.Conductor.Shell(:pipe_wall, copper; t = 0.5)
    pipe=LineCableModels.Enclosure(
        :pipe,
        core;
        shape = LineCableModels.Disk(3.0),
        fill = oil,
        wall
    )
    resolved=DM.resolve(DM.EmptyBoundary(), pipe)
    @test resolved.terminals == [:core]
    @test getproperty.(getproperty.(resolved.regions, :region), :tag) ==
          [:core_metal, :pipe_fill, :pipe_wall]
    @test DM.support(DM.boundary(resolved)) == 3.5
    @test resolved.regions[2].shape.shape == DM.AnnulusShape(1.0, 3.0)

    insulation=LineCableModels.Insulator.Shell(:insulation, oil; t = 0.2)
    screen=LineCableModels.Semiconductor.Shell(:screen, semicon; t = 0.1)
    @test insulation.material.kind === :insulator
    @test screen.material.kind === :semicon
    @test_throws ArgumentError LineCableModels.Conductor.Solid(:bad, oil; r = 1.0)
    @test_throws ArgumentError LineCableModels.Insulator.Shell(:bad, copper; t = 1.0)
    @test_throws ArgumentError LineCableModels.Semiconductor.Shell(:bad, oil; t = 1.0)

    fill=LineCableModels.Filler.Region(
        :fill,
        oil;
        primitive = LineCableModels.Annulus(1.0, 3.0)
    )
    explicit=LineCableModels.Enclosure(
        :duct,
        core;
        shape = LineCableModels.Disk(3.0),
        fill
    )
    @test getproperty.(
        getproperty.(DM.resolve(DM.EmptyBoundary(), explicit).regions, :region),
        :tag
    ) == [:core_metal, :fill]

    materials=LineCableModels.Grid((copper, oil))
    varying=LineCableModels.Conductor.Solid(:varying, materials; r = 1.0)
    @test varying isa LineCableModels.Gridspace{DM.Region}
    @test first(varying).material === copper
    @test_throws ArgumentError collect(varying)

    nested=LineCableModels.Enclosure(
        :outer,
        pipe;
        shape = LineCableModels.Disk(5.0),
        fill = oil
    )
    @test DM.support(DM.boundary(DM.resolve(DM.EmptyBoundary(), nested))) == 5.0
end

@testitem "DataModel / v1 physical tree / placement, Group, and Assembly" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    copper=LineCableModels.Material(kind = :conductor, rho = 1.7241e-8)
    wire=LineCableModels.Region(:wire, LineCableModels.Disk(0.5), copper)

    ring=LineCableModels.Ring(6; r = 2.0, φ0 = π / 6)
    poses=LineCableModels.placements(ring, wire, nothing)
    @test length(poses) == 6
    @test all(pose -> hypot(pose.x, pose.y) ≈ 2.0, poses)
    compacted=LineCableModels.placements(
        ring, wire, LineCableModels.DiameterFactor(0.9)
    )
    @test all(pose -> hypot(pose.x, pose.y) ≈ 1.8, compacted)

    helix=LineCableModels.Helix(LineCableModels.LayRatio(10); dir = -1)
    @test LineCableModels.pitch(helix, 2.0) == 40.0
    @test LineCableModels.angle(helix, 2.0) ≈ atan(4π, 40.0)
    @test LineCableModels.overlength(helix, 2.0) > 1
    fixed=LineCableModels.Helix(LineCableModels.Pitch(20.0))
    @test LineCableModels.pitch(fixed, 3.0) == 20.0
    angled=LineCableModels.Helix(LineCableModels.LayAngle(π / 4))
    @test LineCableModels.pitch(angled, 2.0) ≈ 4π

    strand=LineCableModels.Group(
        :strand,
        wire;
        pattern = LineCableModels.Ring(6; r = 1.0),
        path = helix
    )
    strand_resolution=DM.resolve(DM.EmptyBoundary(), strand)
    @test strand_resolution.terminals == [:strand]
    @test length(strand_resolution.regions) == 6
    @test all(region -> region.terminal === :strand, strand_resolution.regions)
    @test all(region -> region.overlength > 1, strand_resolution.regions)
    @test DM.support(DM.boundary(strand_resolution)) == 1.5

    rope=LineCableModels.Group(
        :phase,
        strand;
        pattern = LineCableModels.Ring(3; r = 4.0)
    )
    rope_resolution=DM.resolve(DM.EmptyBoundary(), rope)
    @test rope_resolution.terminals == [:phase]
    @test length(rope_resolution.regions) == 18
    @test all(region -> region.terminal === :phase, rope_resolution.regions)

    core=LineCableModels.Group(:core, wire)
    phases=LineCableModels.Assembly(
        core;
        pattern = LineCableModels.Ring(3; r = 5.0, φ0 = π / 2),
        names = (:a, :b, :c)
    )
    phase_resolution=DM.resolve(DM.EmptyBoundary(), phases)
    @test phase_resolution.terminals == [:a, :b, :c]
    @test getproperty.(phase_resolution.regions, :terminal) == [:a, :b, :c]

    indexed=LineCableModels.Assembly(
        core;
        pattern = LineCableModels.Lattice(nx = 2, ny = 1, dx = 3.0, dy = 0.0),
        names = :phase
    )
    @test DM.resolve(DM.EmptyBoundary(), indexed).terminals == [:phase_1, :phase_2]
    @test_throws DimensionMismatch DM.resolve(
        DM.EmptyBoundary(),
        LineCableModels.Assembly(core; pattern = ring, names = (:a, :b))
    )
end

@testitem "DataModel / v1 eager construction / concrete boundaries and allocations" tags=[:unit] begin
    using LineCableModels

    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    root=Stack(
        Group(:phase, Region(:core, Disk(0.01), conductor)),
        Region(:insulation, Shell(0.003), dielectric)
    )
    make_design() = CableDesign(root; cable_id = "allocation-baseline")
    design=make_design()
    make_system() = LineCableSystem(
        design;
        position = Pose2(0.0, -1.0),
        connections = (phase = 1,)
    )
    system=make_system()
    problem=LineCableModels.Engine.LineParametersProblem(
        system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0),
        frequencies = [50.0]
    )

    @test isconcretetype(typeof(design))
    @test isconcretetype(typeof(system))
    @test isconcretetype(typeof(problem))
    @test (@inferred LineCableModels.Engine.AnalyticalInput(problem)) isa
          LineCableModels.Engine.AnalyticalInput{Float64}

    # These bounds record the eager-construction scale without turning the
    # physical grammar into an allocation optimization exercise.
    make_design()
    make_system()
    @test @allocated(make_design()) <= 20_000
    @test @allocated(make_system()) <= 12_000
end
