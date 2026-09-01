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
    @test (@inferred DM.resolve(DM.EmptyBoundary(), core)).regions[1].primitive ==
          DM.Disk(0.01)
    @test_throws ArgumentError DM.Region(Symbol(""), DM.Disk(1), copper_props)
    @test_throws ArgumentError DM.Region(:bad, DM.Disk(1), :copper)

    stack=PB.Stack(core, insulation)
    @test stack.items == DM.AbstractCablePart[core, insulation]
    resolved=DM.resolve(DM.EmptyBoundary(), stack)
    @test length(resolved.regions) == 2
    @test resolved.regions[1].primitive == DM.Disk(0.01)
    @test resolved.regions[2].primitive == DM.Annulus(0.01, 0.015)
    @test DM.boundary(resolved) == DM.Disk(0.015)
    @test_throws ArgumentError DM.Stack(DM.AbstractCablePart[])
    @test_throws MethodError DM.resolve(DM.EmptyBoundary(), DM.Stack(insulation))

    nested=PB.Stack(PB.Stack(core), insulation)
    nested_resolved=DM.resolve(DM.EmptyBoundary(), nested)
    @test DM.boundary(nested_resolved) == DM.Disk(0.015)

    radii=PB.Grid((0.01, 0.02))
    cores=PB.Region(:core, PB.Disk(radii), copper_props)
    spaces=PB.Stack(cores, insulation)
    @test spaces isa PB.Gridspace{DM.Stack}
    @test length(spaces) == 2
    @test all(space -> space isa DM.Stack, spaces)
    @test DM.r_ex(DM.boundary(DM.resolve(DM.EmptyBoundary(), first(spaces)))) == 0.015
end

@testitem "DataModel / v1 physical tree / primitive extension is local dispatch" tags=[:unit] begin
    import LineCableModels.DataModel as DM

    struct CapsuleDefinition{T <: Real} <: DM.AbstractPrimitive{T}
        radius::T
    end
    struct Capsule{T <: Real} <: DM.AbstractPrimitive{T}
        radius::T
        at::DM.Pose2{T}
    end

    DM.resolve(::DM.EmptyBoundary, definition::CapsuleDefinition) = Capsule(
        definition.radius,
        DM.Pose2(0, 0, 0)
    )
    DM.resolve(at::DM.Pose2, primitive::Capsule) = Capsule(primitive.radius, at *
                                                                             primitive.at)
    DM.boundary(primitive::Capsule) = primitive
    DM.area(primitive::Capsule) = pi * primitive.radius^2
    DM.centroid(primitive::Capsule) = (primitive.at.x, primitive.at.y)
    DM.support(primitive::Capsule, φ::Real) = primitive.at.x * cos(φ) +
                                              primitive.at.y * sin(φ) + primitive.radius
    DM.support(primitive::Capsule) = hypot(primitive.at.x, primitive.at.y) +
                                     primitive.radius

    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    design=build(
        CableDesign,
        "local-extension",
        Group(:core, Region(:capsule, CapsuleDefinition(0.01), conductor))
    )
    @test only(design.geometry.regions).primitive ==
          Capsule(0.01, DM.Pose2(0, 0, 0))
    @test design.terminal_order == [:core]
end

@testitem "DataModel / v1 physical tree / circular and rectangular strands" tags=[:unit] begin
    copper=Material(kind = :conductor, rho = 1.7241e-8)
    disk=Disk(0.5e-3)
    rectangle=Rectangle(0.35e-3, 0.8e-3)

    circular=stranded(
        copper;
        shape = disk,
        layers = 2,
        n = (6, 12),
        lay = (LayRatio(15), LayRatio(11))
    )
    rectangular=stranded(
        copper;
        shape = rectangle,
        layers = 2,
        n = (6, 12),
        lay = (LayRatio(15), LayRatio(11))
    )
    strand_space=stranded(
        copper;
        shape = Grid((disk, rectangle)),
        layers = 0
    )
    @test circular isa Stack
    @test rectangular isa Stack
    @test strand_space isa Gridspace{Stack}
    resolved_strands=collect(strand_space)
    @test typeof.(getproperty.(first.(getproperty.(resolved_strands, :items)), :item)) ==
          [Region{typeof(disk), typeof(copper)},
        Region{typeof(rectangle), typeof(copper)}]
    @test all(item -> item isa Group, circular.items)
    @test all(item -> item isa Group, rectangular.items)
    @test getproperty.(circular.items, :name) == fill(:strand, 3)
    @test getproperty.(rectangular.items, :name) == fill(:strand, 3)
    @test circular.items[1].path === nothing
    @test circular.items[2].path.lay == LayRatio(15)
    @test circular.items[3].path.lay == LayRatio(11)
    @test getproperty.(circular.items[2:3], :pattern) == [
        Hexa(1),
        Hexa(2)
    ]
    @test all(item -> item.pattern isa Ring, rectangular.items[2:3])

    circular_design=build(
        CableDesign, "circular-strands", terminal(:core, circular)
    )
    rectangular_design=build(
        CableDesign, "rectangular-strands", terminal(:core, rectangular)
    )
    @test length(circular_design.geometry.regions) == 19
    @test length(rectangular_design.geometry.regions) == 19
    @test all(
        region -> region.primitive isa LineCableModels.DataModel.Rectangle,
        rectangular_design.geometry.regions
    )
    @test sum(area, getproperty.(circular_design.geometry.regions, :primitive)) ≈
          19pi * disk.r^2
    @test sum(area, getproperty.(rectangular_design.geometry.regions, :primitive)) ≈
          19rectangle.w * rectangle.h

    strand_centres=centroid.(getproperty.(
        circular_design.geometry.regions,
        :primitive
    ))
    course_two=strand_centres[8:19]
    strand_spacing=2disk.r
    adjacent_distances=[hypot(
                            course_two[mod1(index + 1, length(course_two))][1] -
                            course_two[index][1],
                            course_two[mod1(index + 1, length(course_two))][2] -
                            course_two[index][2]
                        )
                        for index in eachindex(course_two)]
    @test all(distance -> distance ≈ strand_spacing, adjacent_distances)
    all_pair_distances=[hypot(
                            strand_centres[left][1] - strand_centres[right][1],
                            strand_centres[left][2] - strand_centres[right][2]
                        )
                        for left in eachindex(strand_centres)
                        for right in (left + 1):length(strand_centres)]
    @test minimum(all_pair_distances) ≈ strand_spacing
    course_two_radii=hypot.(first.(course_two), last.(course_two))
    @test count(radius -> radius ≈ 2strand_spacing, course_two_radii) == 6
    @test count(radius -> radius ≈ sqrt(3)*strand_spacing, course_two_radii) == 6
    @test [last(region.placement.patterns).member
           for region in circular_design.geometry.regions[8:19]] == collect(1:12)
    flattened_circular=only(LineCableModels.DataModel.flatten(
        circular_design,
        50.0
    )).conductor
    @test flattened_circular.gmr ≈
          LineCableModels.DataModel.BaseParams.strand_gmr(
        strand_centres,
        disk.r,
        copper.mu_r
    )

    initial_radius=hypot(rectangle.w, rectangle.h) / 2
    expected_centre=initial_radius + rectangle.h / 2
    @test hypot(centroid(rectangular_design.geometry.regions[2].primitive)...) ≈
          expected_centre
    @test getproperty.(rectangular_design.geometry.regions, :terminal) ==
          fill(:core, 19)
    @test all(
        region -> region.source.tag === :wire,
        rectangular_design.geometry.regions
    )
    @test all(isempty, getproperty.(rectangular_design.geometry.regions[1:1], :paths))
    @test all(!isempty, getproperty.(rectangular_design.geometry.regions[2:end], :paths))

    for (member, placed) in enumerate(rectangular_design.geometry.regions[2:7])
        x, y = centroid(placed.primitive)
        radial = atan(y, x)
        @test mod(placed.primitive.at.φ - radial, 2pi) ≈ pi / 2
        @test placed.source.primitive === rectangle
    end
    @test_throws DomainError build(
        CableDesign,
        "overpacked",
        terminal(
            :core,
            stranded(
                copper;
                shape = Rectangle(10e-3, 1e-3),
                layers = 1,
                n = (6,)
            )
        )
    )
    @test_throws DimensionMismatch stranded(
        copper; shape = disk, layers = 2, n = (6,)
    )
    @test_throws DimensionMismatch stranded(
        copper; shape = disk, layers = 2, n = (6, 12), lay = (LayRatio(11),)
    )
    elliptical=stranded(
        copper;
        shape = Ellipse(0.5e-3, 0.25e-3),
        layers = 1,
        n = 6
    )
    @test build(CableDesign, "elliptical", terminal(:core, elliptical)) isa CableDesign
    @test_throws ArgumentError LineCableModels.DataModel.flatten(
        rectangular_design, 50.0
    )
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
        primitive = LineCableModels.Disk(3.0),
        fill = oil,
        wall
    )
    resolved=DM.resolve(DM.EmptyBoundary(), pipe)
    @test unique(filter(!isnothing, getproperty.(resolved.regions, :terminal))) == [:core]
    @test getproperty.(getproperty.(resolved.regions, :source), :tag) ==
          [:core_metal, :pipe_fill, :pipe_wall]
    @test DM.support(DM.boundary(resolved)) == 3.5
    @test resolved.regions[2].primitive == DM.Annulus(1.0, 3.0)

    insulation=LineCableModels.Insulator.Shell(:insulation, oil; t = 0.2)
    screen=LineCableModels.Semiconductor.Shell(:screen, semicon; t = 0.1)
    @test insulation.material.kind === :insulator
    @test screen.material.kind === :semicon
    @test_throws ArgumentError LineCableModels.Conductor.Solid(:bad, oil; r = 1.0)
    @test_throws ArgumentError LineCableModels.Insulator.Shell(:bad, copper; t = 1.0)
    @test_throws ArgumentError LineCableModels.Semiconductor.Shell(:bad, oil; t = 1.0)

    fill=LineCableModels.filler(
        oil,
        LineCableModels.Annulus(1.0, 3.0);
        tag = :fill
    )
    explicit=LineCableModels.Enclosure(
        :duct,
        core;
        primitive = LineCableModels.Disk(3.0),
        fill
    )
    @test getproperty.(
        getproperty.(DM.resolve(DM.EmptyBoundary(), explicit).regions, :source),
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
        primitive = LineCableModels.Disk(5.0),
        fill = oil
    )
    @test DM.support(DM.boundary(DM.resolve(DM.EmptyBoundary(), nested))) == 5.0
end

@testitem "DataModel / v1 physical tree / Enclosure material occupancy" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    copper=Material(kind = :conductor, rho = 1.7241e-8)
    matrix=Material(kind = :insulator, rho = Inf, eps_r = 2.3, mu_r = 1.0)
    insulation_material=Material(kind = :insulator, rho = 1.0e12, eps_r = 2.4)
    bundle=terminal(
        :core,
        stranded(copper; shape = Disk(0.5), layers = 1, n = 6)
    )
    packed=Enclosure(
        :matrix,
        bundle;
        primitive = Disk(1.5),
        fill = matrix
    )
    resolved=DM.resolve(DM.EmptyBoundary(), packed)
    fill=last(resolved.regions)

    @test length(resolved.regions) == 8
    @test fill.primitive isa DM._DifferencePrimitive
    @test length(fill.primitive.holes) == 7
    @test DM.area(fill.primitive) ≈ π * 1.5^2 - 7π * 0.5^2
    @test all(region -> region.terminal === :core, resolved.regions[1:7])
    @test fill.terminal === nothing
    @test DM.boundary(resolved) == Disk(1.5)

    annular=Enclosure(
        :annular_matrix,
        Region(:occupied_annulus, Annulus(0.5, 1.0), copper);
        primitive = Disk(1.5),
        fill = matrix
    )
    annular_fill=last(DM.resolve(DM.EmptyBoundary(), annular).regions)
    @test annular_fill.primitive isa DM._DifferencePrimitive
    @test DM.area(annular_fill.primitive) ≈
          π * 1.5^2 - π * (1.0^2 - 0.5^2)

    armor=Group(
        :armor,
        Region(:armor_wires, Disk(0.1), copper);
        pattern = Ring(6; r = 1.1)
    )
    annular_stage=Enclosure(
        :armor_matrix,
        armor;
        primitive = Annulus(1.0, 1.2),
        fill = matrix
    )
    staged_design=build(
        CableDesign,
        "annular-filled-stage",
        Group(:core, Region(:inner_core, Disk(1.0), copper)),
        annular_stage
    )
    staged_fill=last(staged_design.geometry.regions)
    @test staged_fill.primitive isa DM._DifferencePrimitive
    @test staged_fill.primitive.outer isa Annulus
    @test length(staged_fill.primitive.holes) == 6
    @test sum(DM.area, staged_design.geometry.regions) ≈
          DM.area(staged_design.geometry.outer)
    @test DM.boundary(staged_design.geometry) == Disk(1.2)

    discontinuous_stage=Enclosure(
        :discontinuous_matrix,
        armor;
        primitive = Annulus(1.01, 1.21),
        fill = matrix
    )
    @test_throws DomainError build(
        CableDesign,
        "discontinuous-annular-stage",
        Group(:core, Region(:inner_core, Disk(1.0), copper)),
        discontinuous_stage
    )

    bare_design=build(CableDesign, "bare-bundle", bundle)
    packed_design=build(CableDesign, "packed-bundle", packed)
    bare=only(DM.flatten(bare_design, 50.0))
    reduced=only(DM.flatten(packed_design, 50.0))
    @test reduced.conductor == bare.conductor
    @test reduced.dielectric == bare.dielectric

    IE=LineCableModels.ImportExport
    restored=IE.deserialize_value(IE.serialize_value(packed_design))
    @test restored == packed_design
    restored_fill=last(restored.geometry.regions).primitive
    @test restored_fill isa DM._DifferencePrimitive
    @test length(restored_fill.holes) == 7

    insulated=build(
        CableDesign,
        "packed-insulated",
        packed,
        Region(:insulation, Shell(0.2), insulation_material)
    )
    component=only(DM.flatten(insulated, 50.0))
    @test length(component.dielectric.layers) == 1
    @test only(component.dielectric.layers).r_in ≈ 1.5
    @test only(component.dielectric.layers).r_ex ≈ 1.7

    nested=Group(:core, packed; pattern = Ring(2; r = 4.0))
    nested_design=build(CableDesign, "nested-packed", nested)
    @test only(DM.flatten(nested_design, 50.0)).conductor.num_wires == 14

    lossy=Material(kind = :insulator, rho = 1.0e12, eps_r = 2.3, mu_r = 1.0)
    magnetic=Material(kind = :insulator, rho = Inf, eps_r = 2.3, mu_r = 2.0)
    semiconducting=Material(kind = :semicon, rho = 1.0e3, eps_r = 10.0)
    lossy_design=build(
        CableDesign,
        "lossy-matrix",
        Enclosure(:matrix, bundle; primitive = Disk(1.5), fill = lossy)
    )
    magnetic_design=build(
        CableDesign,
        "magnetic-matrix",
        Enclosure(:matrix, bundle; primitive = Disk(1.5), fill = magnetic)
    )
    semiconducting_design=build(
        CableDesign,
        "semiconducting-matrix",
        Enclosure(:matrix, bundle; primitive = Disk(1.5), fill = semiconducting)
    )
    wound_matrix_design=build(
        CableDesign,
        "longitudinal-matrix-path",
        Group(:wound_core, packed; path = Helix(LayRatio(10.0)))
    )
    @test_throws ArgumentError DM.flatten(lossy_design, 50.0)
    @test_throws ArgumentError DM.flatten(magnetic_design, 50.0)
    @test_throws ArgumentError DM.flatten(semiconducting_design, 50.0)
    @test_throws ArgumentError DM.flatten(wound_matrix_design, 50.0)
end

@testitem "DataModel / v1 physical tree / placement, Group, and Assembly" tags=[:unit] begin
    const DM=LineCableModels.DataModel
    copper=LineCableModels.Material(kind = :conductor, rho = 1.7241e-8)
    wire=LineCableModels.Region(:wire, LineCableModels.Disk(0.5), copper)

    ring=LineCableModels.Ring(6; r = 2.0, φ0 = π / 6)
    poses=LineCableModels.placements(ring, wire.primitive, nothing)
    @test length(poses) == 6
    @test all(pose -> hypot(pose.x, pose.y) ≈ 2.0, poses)
    @test_throws MethodError LineCableModels.placements(
        ring, wire.primitive, LineCableModels.DiameterFactor(0.9)
    )

    for course in 1:4
        hexagonal=Hexa(course; φ0 = π / 6)
        hexagonal_poses=placements(hexagonal, wire.primitive, nothing)
        @test capacity(hexagonal, wire.primitive, nothing) == 6course
        @test length(hexagonal_poses) == 6course
        @test all(eachindex(hexagonal_poses)) do index
            next=hexagonal_poses[mod1(index + 1, length(hexagonal_poses))]
            current=hexagonal_poses[index]
            hypot(next.x - current.x, next.y - current.y) ≈ 1.0
        end
    end
    hexagonal=Hexa(2; φ0 = π / 6)
    hexagonal_poses=placements(hexagonal, wire.primitive, nothing)
    @test length(unique(round.(
        hypot.(getfield.(hexagonal_poses, :x), getfield.(hexagonal_poses, :y));
        digits = 12
    ))) == 2
    @test_throws ArgumentError Hexa(0)
    @test_throws DomainError Hexa(1; gap_frac = -0.1)
    @test_throws MethodError placements(hexagonal, Rectangle(1.0, 1.0), nothing)
    course_space=Hexa(Grid((1, 2)))
    @test course_space isa Gridspace{Hexa}
    @test collect(course_space) == [Hexa(1), Hexa(2)]

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
    @test unique(getproperty.(strand_resolution.regions, :terminal)) == [:strand]
    @test length(strand_resolution.regions) == 6
    @test all(region -> region.terminal === :strand, strand_resolution.regions)
    @test all(
        region -> prod(
            entry -> LineCableModels.overlength(entry.path, entry.radius),
            region.paths;
            init = 1.0
        ) > 1,
        strand_resolution.regions
    )
    @test DM.support(DM.boundary(strand_resolution)) == 1.5

    rope=LineCableModels.Group(
        :phase,
        strand;
        pattern = LineCableModels.Ring(3; r = 4.0)
    )
    rope_resolution=DM.resolve(DM.EmptyBoundary(), rope)
    @test unique(getproperty.(rope_resolution.regions, :terminal)) == [:phase]
    @test length(rope_resolution.regions) == 18
    @test all(region -> region.terminal === :phase, rope_resolution.regions)

    core=LineCableModels.Group(:core, wire)
    phases=LineCableModels.Assembly(
        core;
        pattern = LineCableModels.Ring(3; r = 5.0, φ0 = π / 2),
        names = (:a, :b, :c)
    )
    phase_resolution=DM.resolve(DM.EmptyBoundary(), phases)
    @test getproperty.(phase_resolution.regions, :terminal) == [:a, :b, :c]

    indexed=LineCableModels.Assembly(
        core;
        pattern = LineCableModels.Lattice(nx = 2, ny = 1, dx = 3.0, dy = 0.0),
        names = (:phase_1, :phase_2)
    )
    @test getproperty.(
        DM.resolve(DM.EmptyBoundary(), indexed).regions,
        :terminal
    ) == [:phase_1, :phase_2]
    @test_throws DimensionMismatch DM.resolve(
        DM.EmptyBoundary(),
        LineCableModels.Assembly(core; pattern = ring, names = (:a, :b))
    )
end

@testitem "DataModel / v1 construction / concrete boundaries and allocations" tags=[:unit] begin
    using LineCableModels

    conductor=Material(kind = :conductor, rho = 1.7241e-8)
    dielectric=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    root=Stack(
        Group(:phase, Region(:core, Disk(0.01), conductor)),
        Region(:insulation, Shell(0.003), dielectric)
    )
    make_design() = build(CableDesign, "allocation-baseline", root)
    design=make_design()
    make_system() = build(
        LineCableSystem,
        design,
        Pose2(0.0, -1.0);
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
    execution=computation_options(Val(LineCableModelsCoaxial), (;))
    blueprints=LineCableModels.Engine.CableBlueprint{eltype(problem)}[
        LineCableModels.Engine.flatten(LineCableModelsCoaxial(), source, eltype(problem))
        for source in problem.system.designs
    ]
    @test (@inferred LineCableModels.Engine.LineParametersWorkspace(
        problem,
        Formulation(),
        execution,
        blueprints
    )) isa LineCableModels.Engine.LineParametersWorkspace{Float64}

    # These bounds record the construction scale without turning the
    # physical grammar into an allocation optimization exercise.
    make_design()
    make_system()
    @test @allocated(make_design()) <= 20_000
    @test @allocated(make_system()) <= 12_000
end
