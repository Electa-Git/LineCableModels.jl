@testitem "DataModel / assembly placement preserves local conductor lowering" tags=[:unit] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    const IE = LineCableModels.ImportExport
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
    dielectric = Material(:insulator, 1e14, 2.3, 1.0)
    engine = LineCableModelsCoaxial()
    wire = Disk(0.5e-3)
    single = terminal(:core, solid(copper, wire))
    placed = at(single, Pose2(0.01, 0.0))
    for constructor in (assembly, cores)
        @test only(constructor(single).item).item == single
        @test only(constructor(placed).item) == placed
        gridded = constructor(at(single, Pose2(Grid((0.01, 0.02)), 0.0)))
        @test gridded isa Gridspace{Assembly}
        @test length(gridded) == 2
        @test [only(item.item).at.x for item in gridded] == [0.01, 0.02]
    end
    @test_throws ArgumentError assembly(single, single; pattern=Ring(2; r=0.01))
    @test_throws ArgumentError assembly(single; names=(:a,))
    @test_throws ArgumentError cores(single; r=0.01)
    @test_throws ArgumentError cores(single; n=3, names=(:a, :b, :c))
    bodies = (
        solid(copper, Disk(1e-3)),
        stranded(copper; shape=wire, boundary=Disk(5wire.r), lay=LayRatio(12)),
        stranded(copper; shape=wire, boundary=Disk(sqrt(19) * wire.r), compact=true),
        stranded(copper; shape=Rectangle(0.3e-3, 0.1e-3), center=Disk(0.2e-3),
            boundary=Disk(0.59e-3)),
        stranded(copper; shape=Disk(0.35e-3),
            boundary=Sector(span=2pi/3, r_base=0.6e-3, r_back=4e-3, fillet=0.2e-3)),
    )
    for body in bodies, count in (1, 3)
        names = ntuple(index -> Symbol(:phase_, index), count)
        member = terminal(:core, body, insulation(dielectric; t=0.5e-3))
        pattern = Ring(count; r=0.01, φ0=0.27)
        repeated = build(CableDesign, "repeated", assembly(member; pattern, names))
        poses = DM.placements(pattern, DM.resolve(EmptyBoundary(), member), nothing)
        explicit = build(CableDesign, "explicit", assembly((
            at(terminal(name, body, insulation(dielectric; t=0.5e-3)), pose)
            for (name, pose) in zip(names, poses))...))
        actual = @inferred EN.flatten(engine, repeated, Float64)
        expected = @inferred EN.flatten(engine, explicit, Float64)
        @test actual.assembly_ranges == expected.assembly_ranges == [i:i for i in 1:count]
        for (left, right) in zip(actual.conductors, expected.conductors)
            @test left.terminal === right.terminal
            for field in (:r_in, :r_ex, :cross_section, :num_wires, :num_turns,
                    :resistance, :alpha, :gmr, :position)
                @test all(isapprox.(getproperty(left, field), getproperty(right, field)))
            end
        end
        @test actual.dielectrics == expected.dielectrics
        @test IE.deserialize_value(IE.serialize_value(repeated)) == repeated
        @test all(region -> last(region.placement.patterns).owner === DM.Assembly,
            repeated.geometry.regions)
        relocated = DM.resolve(Pose2(0.004, -0.006, 0.31), first(repeated.geometry.regions))
        @test last(relocated.placement.patterns).owner === DM.Assembly
        @test last(relocated.placement.patterns).member == 1
    end

    # The same Assembly becomes a strand course when a Group coalesces it.
    centre = terminal(:core, solid(copper, wire))
    member = terminal(:wire, solid(copper, wire))
    pattern = Ring(6; r=2wire.r)
    outer = assembly(member; pattern, names=ntuple(index -> Symbol(:wire_, index), 6))
    merged = terminal(:core, assembly(centre, outer))
    direct = Group(:core, Region(:solid, wire, copper); pattern)
    merged_design = build(CableDesign, "coalesced", merged)
    direct_design = build(CableDesign, "direct", centre, direct)
    @test DM.radial_components(merged_design, Float64) == DM.radial_components(direct_design, Float64)
    @test all(region -> all(entry -> entry.owner !== DM.Assembly, region.placement.patterns),
        merged_design.geometry.regions)

    # A helical placement retains the physical overlength independently of
    # whether its child is a strand or a complete coaxial core.
    helix = Helix(LayRatio(12))
    helical = build(CableDesign, "helical-assemblies", assembly(member;
        pattern=Ring(3; r=0.01), names=(:a, :b, :c), path=helix))
    rows = EN.flatten(engine, helical, Float64).conductors
    @test all(row -> row.num_wires == 0, rows)
    @test all(row -> row.resistance ≈ copper.rho / area(wire) * overlength(helix, 0.01), rows)
    @test all(region -> only(region.paths).path === helix, helical.geometry.regions)
    @test all(region -> only(region.paths).radius ≈ 0.01, helical.geometry.regions)
    prototype = terminal(:core, bodies[2])
    isolated = build(CableDesign, "stranded", prototype)
    helical_stranded = build(CableDesign, "helical-stranded", assembly(prototype;
        pattern=Ring(3; r=0.01), names=(:a, :b, :c), path=helix))
    base = only(EN.flatten(engine, isolated, Float64).conductors)
    wound = EN.flatten(engine, helical_stranded, Float64).conductors
    @test all(row -> row.resistance ≈ base.resistance * overlength(helix, 0.01), wound)
end

@testitem "Engine / repeated coaxials agree with explicit placements and remain Gridable" tags=[:integration] begin
    const DM = LineCableModels.DataModel
    const EN = LineCableModels.Engine
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
    dielectric = Material(:insulator, 1e14, 2.3, 1.0)
    parts = (solid(copper, Disk(1e-3)), insulation(dielectric; t=1e-3))
    member = terminal(:core, parts...)
    names = (:a, :b, :c)
    spaces = build(CableDesign, "assembly-grid", cores(member;
        n=3, r=Grid((0.01, 0.02)), names))
    @test spaces isa Gridspace{CableDesign}
    @test length(spaces) == 2
    for (design, radius) in zip(spaces, (0.01, 0.02))
        poses = DM.placements(Ring(3; r=radius), DM.resolve(EmptyBoundary(), member), nothing)
        explicit = build(CableDesign, "explicit", assembly((
            at(terminal(name, parts...), pose) for (name, pose) in zip(names, poses))...))
        systems = map((design, explicit)) do source
            build(LineCableSystem, [source], [Pose2(0.0, -1.0)];
                connections=[Dict(:a=>1, :b=>2, :c=>3)])
        end
        problems = map(system -> LineParametersProblem(system;
            earth_props=homogeneous(rho=100.0), frequencies=[50.0, 1000.0]), systems)
        actual = @inferred compute(first(problems))
        expected = @inferred compute(last(problems))
        @test actual.Z.values == expected.Z.values
        @test actual.Y.values == expected.Y.values
        @test EN.flatten(LineCableModelsCoaxial(), design, Float64).assembly_ranges == [1:1, 2:2, 3:3]
    end
    # Combining the core declaration and placement axes uses the same
    # Gridspace rule, including a zip across the nested Ring declaration.
    radii = (0.01, 0.02)
    angles = (0.0, 0.3)
    wires = (0.5e-3, 1e-3)
    members = terminal(:core, solid(copper, Disk(Grid(wires))))
    zipped = cores(members; n=3, r=Grid(radii), φ0=Grid(angles), names, combine=:zip)
    @test zipped isa Gridspace{Assembly}
    @test collect(zipped) == [cores(terminal(:core, solid(copper, Disk(wire)));
        n=3, r=radius, φ0=angle, names) for (wire, radius, angle) in zip(wires, radii, angles)]
    product = cores(members; n=3, r=Grid(radii), φ0=Grid(angles), names)
    @test length(product) == 8
    @test all(cores(terminal(:core, solid(copper, Disk(wire)));
        n=3, r=radius, φ0=angle, names) in product
        for wire in wires, radius in radii, angle in angles)
end

@testitem "Gmsh FEM / assembly placement does not coarsen solid-core mesh size" tags=[:extension] begin
    using Gmsh
    const DM = LineCableModels.DataModel
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
    for primitive in (Disk(1e-3), Rectangle(1e-3, 0.5e-3), Ellipse(1e-3, 0.5e-3))
        member = terminal(:core, solid(copper, primitive))
        isolated = build(CableDesign, "isolated", member)
        repeated = build(CableDesign, "repeated", cores(member; n=3, r=0.01, names=(:a, :b, :c)))
        expected = extension._fem_region_mesh_size(only(isolated.geometry.regions))
        @test all(region -> extension._fem_region_mesh_size(region) == expected, repeated.geometry.regions)
    end
end
