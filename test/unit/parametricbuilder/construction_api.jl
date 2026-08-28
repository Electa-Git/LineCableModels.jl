@testitem "ParametricBuilder / construction API / lowering and terminal ownership" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.DataModel as DM

    copper = Material(kind = :conductor, rho = 1.72e-8)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    semicon = Material(kind = :semicon, rho = 100.0, eps_r = 20.0)

    phase = terminal(
        :phase,
        core(copper; r = 5e-3),
        screen(semicon; t = 0.5e-3),
        insulation(xlpe; t = 3e-3)
    )
    explicit = Group(
        :phase,
        layers(
            Region(:core, Disk(5e-3), copper),
            Region(:screen, Shell(0.5e-3), semicon),
            Region(:insulation, Shell(3e-3), xlpe)
        )
    )
    @test phase == explicit

    functional = build(CableDesign, "coaxial", phase)
    notated = @cable "coaxial" begin
        @terminal :phase begin
            core(copper; r = 5e-3)
            screen(semicon; t = 0.5e-3)
            insulation(xlpe; t = 3e-3)
        end
    end
    @test notated == functional
    @test functional.terminal_order == [:phase]
    @test functional.root isa Group
    @test all(==(1), filter(!iszero, functional.terminal_map))
    @test_throws ArgumentError build(
        CableDesign,
        "undeclared-terminal",
        core(copper; r = 5e-3)
    )

    insulating_tapes = tape(
        xlpe;
        section = Sector(8.5e-3, 9e-3, -0.08, 0.16),
        n = 12,
        lay = LayRatio(10)
    )
    insulating_resolution = DM.resolve(DM.EmptyBoundary(), insulating_tapes)
    @test all(isnothing, getproperty.(insulating_resolution.regions, :terminal))
    @test_throws ArgumentError terminal(:insulating_tapes, insulating_tapes)
    with_tapes = build(CableDesign, "taped", phase, insulating_tapes)
    @test with_tapes.terminal_order == [:phase]
    @test count(iszero, with_tapes.terminal_map) > 0

    distributed_tapes = @distribute tape(
        xlpe;
        section = Sector(8.5e-3, 9e-3, -0.08, 0.16),
        gap_frac = 0.02,
        lay = LayRatio(10)
    )
    @test distributed_tapes.pattern.n == capacity()
    distributed_design = build(
        CableDesign,
        "distributed-tapes",
        phase,
        distributed_tapes
    )
    tape_regions = filter(
        region -> region.source.tag === :tape,
        distributed_design.geometry.regions
    )
    @test length(tape_regions) == capacity(
        Ring(capacity(); r = 0, gap_frac = 0.02),
        Sector(8.5e-3, 9e-3, -0.08, 0.16),
        nothing
    )
    @test all(region -> length(region.paths) == 1, tape_regions)

    nested = build(
        CableDesign,
        "nested-terminal",
        terminal(:outer, terminal(:inner, core(copper; r = 2e-3)))
    )
    @test nested.terminal_order == [:outer]

    sector_a = terminal(
        :a,
        core(copper, Sector(0, 3e-3, -pi / 6, pi / 3)),
        insulation(xlpe; t = 0.5e-3)
    )
    sector_b = terminal(
        :b,
        core(copper, Sector(0, 2.5e-3, -pi / 8, pi / 4)),
        insulation(xlpe; t = 0.4e-3)
    )
    explicit_assembly = assembly(
        at(sector_a, -5e-3, 0; φ = 0.2),
        at(sector_b, 5e-3, 0; φ = -0.3)
    )
    notated_assembly = @assembly begin
        @at sector_a (-5e-3, 0, 0.2)
        @at sector_b (5e-3, 0, -0.3)
    end
    @test notated_assembly == explicit_assembly
    multicore = build(CableDesign, "sectorized", explicit_assembly)
    @test multicore.terminal_order == [:a, :b]
    @test multicore.geometry.regions[1].primitive isa DM.Sector

    separate_screen = build(
        CableDesign,
        "screened",
        terminal(:phase, core(copper; r = 2e-3), insulation(xlpe; t = 1e-3)),
        terminal(:screen, screen(copper; t = 0.2e-3))
    )
    @test separate_screen.terminal_order == [:phase, :screen]
end

@testitem "ParametricBuilder / construction API / public surface and macro lowering" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.DataModel as DM
    import LineCableModels.ParametricBuilder as PB

    required = (
        :build, :terminal, :core, :strand, :rope, :cores, :tape,
        :insulation, :screen, :sheath, :armor, :bedding, :jacket, :filler,
        :pipe, :duct, :at, :trefoil, :hflat, :vflat, :capacity,
        :solid, :shell, :wires, :layers, :assembly,
        :Region, :Stack, :Group, :Assembly, :Enclosure, :Pose2,
        :Disk, :Rectangle, :Ellipse,
        :Sector, :Annulus, :Shell,
        :Polygon, :Ring, :Polar, :Fill, :Lattice,
        :Helix, :LayRatio, :Pitch, :LayAngle, :FillFactor,
        :DiameterFactor, :TabulatedCompaction, :AffineCompaction,
        Symbol("@cable"), Symbol("@terminal"), Symbol("@assembly"),
        Symbol("@duct"), Symbol("@at"), Symbol("@hflat"), Symbol("@vflat"),
        Symbol("@trefoil"), Symbol("@distribute")
    )
    public_names = Set(names(LineCableModels))
    @test all(name -> name in public_names, required)

    forbidden = (
        :DuctBank, :Pipe, :Duct, :Core, :Cable, :Trefoil, :distribute,
        Symbol("@design"), Symbol("@pipe"), Symbol("@core"),
        Symbol("@insulation"), Symbol("@strand"), Symbol("@rope"),
        Symbol("@armor"), Symbol("@screen"), Symbol("@shell"),
        Symbol("@wire"), Symbol("@fill")
    )
    @test all(name -> !(name in public_names), forbidden)
    @test LineCableModels.build === DM.build === PB.build
    @test all(name -> name in public_names,
        (:Disk, :Rectangle, :Ellipse, :Sector, :Annulus, :Polygon))

    expansions = (
        :(@cable "x" begin
            part
        end) => :build,
        :(@terminal :phase begin
            part
        end) => :terminal,
        :(@assembly begin
            left
            right
        end) => :assembly,
        :(@duct shape=shape fill=fill begin
            left
            right
        end) => :duct,
        :(@at subject (1, 2) φ=3) => :at,
        :(@trefoil design spacing=spacing center=(0, 0) phase=(1, 2, 3)) =>
            :trefoil,
        :(@hflat design spacing=spacing phase=(1, 2, 3)) => :hflat,
        :(@vflat design spacing=spacing phase=(1, 2, 3)) => :vflat
    )
    for (surface, callee) in expansions
        expanded = Base.remove_linenums!(macroexpand(@__MODULE__, surface))
        @test expanded isa Expr && expanded.head === :call
        @test expanded.args[1] == GlobalRef(PB, callee)
    end

    distributed = Base.remove_linenums!(macroexpand(
        @__MODULE__,
        :(@distribute wires(material; wire = wire, r = radius))
    ))
    @test distributed.head === :call
    @test distributed.args[1] === :wires
    parameters = only(filter(
        argument -> argument isa Expr && argument.head === :parameters,
        distributed.args
    ))
    deferred = only(filter(
        argument -> argument isa Expr && argument.head === :kw &&
                    argument.args[1] === :n,
        parameters.args
    ))
    @test deferred.args[2] == Expr(:call, GlobalRef(PB, :capacity))

    copper = Material(kind = :conductor, rho = 1.72e-8)
    evaluations = Dict{Symbol, Int}()
    once(name, value) = (evaluations[name] = get(evaluations, name, 0) + 1; value)
    nested = @cable once(:identifier, "single-evaluation") begin
        @terminal once(:terminal, :phase) begin
            core(
                once(:material, copper);
                r = once(:radius, 1e-3)
            )
        end
    end
    @test nested isa CableDesign
    @test evaluations == Dict(
        :identifier => 1,
        :terminal => 1,
        :material => 1,
        :radius => 1
    )
end

@testitem "ParametricBuilder / construction API / macro arguments evaluate once" tags=[:unit] begin
    using LineCableModels

    copper = Material(kind = :conductor, rho = 1.72e-8)
    air = Material(kind = :insulator, rho = 1.0e16, eps_r = 1.0)
    counts = Dict{Symbol, Int}()
    once(name, value) = (counts[name] = get(counts, name, 0) + 1; value)

    design = @cable once(:cable_id, "macro-evaluation") begin
        @terminal once(:terminal_name, :phase) begin
            core(once(:core_material, copper); r = once(:core_radius, 1e-3))
        end
    end
    member = terminal(:aux, core(copper; r = 0.5e-3))
    placed_member = at(member, 0.0, 0.0)

    @assembly begin
        once(:assembly_first, placed_member)
        once(:assembly_second, at(member, 2e-3, 0.0))
    end
    @duct shape=once(:duct_shape, Disk(5e-3)) fill=once(:duct_fill, air) begin
        once(:duct_member, placed_member)
    end
    @at once(:at_subject, design) (
        once(:at_x, 0.0),
        once(:at_y, -1.0),
        once(:at_angle, 0.2)
    ) connections = once(:at_connections, (phase = 1,))

    @trefoil once(:trefoil_design, design) spacing=once(:trefoil_spacing, 0.1) center=(
        once(:trefoil_x, 0.0), once(:trefoil_y, -1.0)
    ) phase=once(:trefoil_connections, (1, 2, 3))
    @hflat once(:hflat_design, design) spacing=once(:hflat_spacing, 0.1) phase=once(
        :hflat_connections, (1, 2, 3))
    @vflat once(:vflat_design, design) spacing=once(:vflat_spacing, 0.1) phase=once(
        :vflat_connections, (1, 2, 3))

    @distribute wires(
        once(:wire_material, copper);
        wire = once(:wire_definition, Disk(0.2e-3)),
        r = once(:wire_radius, 2e-3)
    )

    @test !isempty(counts)
    @test all(==(1), values(counts))
end

@testitem "ParametricBuilder / construction API / enclosures and macro transparency" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.DataModel as DM

    copper = Material(kind = :conductor, rho = 1.72e-8)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    air = Material(kind = :insulator, rho = 1.0e16, eps_r = 1.0)
    concrete = Material(kind = :insulator, rho = 100.0, eps_r = 4.0)

    cell_a = duct(
        terminal(:a, core(copper; r = 1e-3), insulation(xlpe; t = 0.5e-3));
        shape = Disk(4e-3),
        fill = air
    )
    cell_b = duct(
        terminal(:b, core(copper; r = 1.2e-3), insulation(xlpe; t = 0.4e-3));
        shape = Disk(4e-3),
        fill = air
    )
    asymmetric = duct(
        at(cell_a, -6e-3, 0),
        at(cell_b, 6e-3, 1e-3; φ = 0.2);
        shape = Rectangle(20e-3, 12e-3),
        fill = concrete
    )
    notated = @duct shape=Rectangle(20e-3, 12e-3) fill=concrete begin
        @at cell_a (-6e-3, 0)
        @at cell_b (6e-3, 1e-3, 0.2)
    end
    @test notated == asymmetric
    @test asymmetric.item isa Assembly
    design = build(CableDesign, "asymmetric-ducts", asymmetric)
    @test design.terminal_order == [:a, :b]
    @test any(region -> region.source.tag === :duct_fill, design.geometry.regions)

    repeated = duct(
        cell_a;
        formation = Lattice(nx = 2, ny = 1, dx = 12e-3, dy = 0),
        shape = Rectangle(24e-3, 12e-3),
        fill = concrete
    )
    @test repeated.item isa Assembly
    @test repeated.item.item === cell_a
    @test repeated.item.pattern == Lattice(nx = 2, ny = 1, dx = 12e-3, dy = 0)

    evaluations = Symbol[]
    identifier() = (push!(evaluations, :identifier); "single-evaluation")
    phase() = (push!(evaluations, :phase); terminal(:phase, core(copper; r = 1e-3)))
    built = @cable identifier() begin
        "This literal is documentation, not a cable part."
        nothing
        phase()
    end
    @test built isa CableDesign
    @test evaluations == [:identifier, :phase]

    expanded = Base.remove_linenums!(macroexpand(
        @__MODULE__,
        :(@terminal :phase begin
            core(copper; r = 1e-3)
            insulation(xlpe; t = 0.5e-3)
        end)
    ))
    @test expanded isa Expr && expanded.head === :call
    @test expanded.args[1] == GlobalRef(LineCableModels.ParametricBuilder, :terminal)
    @test expanded.args[2] == QuoteNode(:phase)
    @test length(expanded.args) == 4
end

@testitem "ParametricBuilder / construction API / repetition and compaction" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.DataModel as DM

    copper = Material(kind = :conductor, rho = 1.72e-8)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)

    round = strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 2,
        n = 6,
        lay = (LayRatio(12), Pitch(0.2)),
        dir = (1, -1),
        φ0 = (0.0, 0.1)
    )
    @test round isa Stack
    @test length(round.items) == 3
    @test round.items[1].pattern === nothing
    @test round.items[1].path === nothing
    @test [item.pattern.n for item in round.items[2:3]] == [6, 12]
    @test getfield.(round.items[2:3], :path) == [
        Helix(LayRatio(12); dir = 1, φ0 = 0.0),
        Helix(Pitch(0.2); dir = -1, φ0 = 0.1)
    ]

    exact = strand(
        copper;
        wire = Rectangle(0.6e-3, 0.8e-3),
        layers = 2,
        n = (5, 11),
        lay = nothing
    )
    @test [item.pattern.n for item in exact.items[2:3]] == [5, 11]
    rectangular = build(
        CableDesign,
        "rectangular-strand",
        terminal(:phase, exact),
        insulation(xlpe; t = 1e-3)
    )
    rectangles = filter(
        region -> region.primitive isa DM.Rectangle,
        rectangular.geometry.regions
    )
    @test length(rectangles) == 17
    @test all(region -> DM.area(region) == 0.6e-3 * 0.8e-3, rectangles)

    compacted = strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 1,
        n = 6,
        compact = FillFactor(0.9),
        lay = LayRatio(11)
    )
    compacted_design = build(
        CableDesign,
        "compacted",
        terminal(:phase, compacted)
    )
    sectors = filter(
        region -> region.primitive isa DM.Sector,
        compacted_design.geometry.regions
    )
    @test length(sectors) == 6
    @test sum(DM.area, sectors) ≈ 6pi * (0.5e-3)^2
    @test all(region -> length(region.paths) == 1, sectors)
    expected_outer = sqrt((0.5e-3)^2 + 6 * (0.5e-3)^2 / 0.9)
    @test all(region -> region.primitive.ro ≈ expected_outer, sectors)
    @test all(
        region -> only(region.paths).radius ≈
                  (0.5e-3 + expected_outer) / 2,
        sectors
    )

    mixed = strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 2,
        n = (6, capacity()),
        compact = (nothing, FillFactor(0.9))
    )
    mixed_design = build(CableDesign, "mixed-capacity", terminal(:phase, mixed))
    resolved_counts = unique([last(region.placement.patterns).pattern.n
                              for region in mixed_design.geometry.regions
                              if !isempty(region.placement.patterns)])
    @test first(resolved_counts) == 6
    @test last(resolved_counts) > 6

    direct_capacity = capacity(Ring, 10e-3, 0.5e-3; gap_frac = 0.03)
    pattern = Ring(capacity(); r = 10e-3, gap_frac = 0.03)
    @test capacity(pattern, Disk(0.5e-3), nothing) == direct_capacity
    distributed = @distribute wires(
        copper;
        wire = Disk(0.5e-3),
        r = 10e-3,
        gap_frac = 0.03
    )
    @test distributed.pattern.n == capacity()
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@distribute wires(copper; wire = Disk(1e-3), n = 6, r = 2e-3))
    )

    fill_pattern = Fill(r = 2.6e-3, φ = pi / 6)
    filled_poses = placements(fill_pattern, Disk(0.5e-3), nothing)
    @test length(filled_poses) == 19
    @test capacity(fill_pattern, Disk(0.5e-3), nothing) == 19
    @test first(filled_poses) == Pose2(0, 0, 0)
    @test maximum(pose -> hypot(pose.x, pose.y), filled_poses) + 0.5e-3 <=
          fill_pattern.r + eps(fill_pattern.r)

    nested = rope(
        round;
        layers = 1,
        n = 6,
        lay = LayRatio(9)
    )
    nested_design = build(
        CableDesign,
        "nested-rope",
        terminal(:phase, nested)
    )
    @test any(region -> length(region.paths) == 2, nested_design.geometry.regions)
    @test any(region -> isempty(region.paths), nested_design.geometry.regions)
    leaf_resistances = map(nested_design.geometry.regions) do region
        primitive = region.primitive
        base = DM.tubular_resistance(
            0.0,
            primitive.r,
            region.source.material.rho,
            region.source.material.alpha,
            region.source.material.T0,
            region.source.material.T0
        )
        foldl(region.paths; init = base) do resistance, entry
            resistance * overlength(entry.path, entry.radius)
        end
    end
    expected_resistance = inv(sum(inv, leaf_resistances))
    component = only(LineCableModels.Engine.homogeneous_components(
        Formulation(),
        nested_design,
        50.0
    ))
    @test component.conductor.resistance ≈ expected_resistance

    @test_throws DimensionMismatch strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 2,
        lay = (LayRatio(10),)
    )
    @test_throws MethodError DM.placements(
        Ring(6; r = 2e-3),
        Disk(0.5e-3),
        DiameterFactor(0.9)
    )
end

@testitem "ParametricBuilder / construction API / placement, macros, and Gridspace" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.DataModel as DM

    copper = Material(kind = :conductor, rho = 1.72e-8)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    design = @cable "formation-cable" begin
        @terminal :phase begin
            core(copper; r = 3e-3)
            insulation(xlpe; t = 1e-3)
        end
    end

    explicit = trefoil(
        design;
        center = at(0.2, -1.0),
        spacing = 0.1,
        connections = (phase = (1, 2, 3),)
    )
    notated = @trefoil design spacing=0.1 center=(0.2, -1.0) phase=(1, 2, 3)
    @test notated == explicit
    @test getproperty.(getproperty.(explicit, :connections), :phase) == [1, 2, 3]
    @test all(placement -> placement.design === design, explicit)

    moved = at(explicit, 1.0, 2.0; φ = 0.25)
    outer = at(1.0, 2.0; φ = 0.25)
    @test getproperty.(moved, :pose) == [outer * item.pose for item in explicit]
    @test all(placement -> placement.design === design, moved)
    @test design.geometry == first(moved).design.geometry

    explicit_connections = (phase = (1, 2, 3),)
    @test (@hflat design spacing=0.1 center=(0.2, -1.0) connections=explicit_connections) ==
          hflat(
        design;
        spacing = 0.1,
        center = at(0.2, -1.0),
        connections = explicit_connections
    )
    @test (@vflat design spacing=0.1 center=(0.2, -1.0) phase=(1, 2, 3)) ==
          vflat(
        design;
        spacing = 0.1,
        center = at(0.2, -1.0),
        connections = explicit_connections
    )
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@trefoil design spacing=0.1 connections=explicit_connections phase=(1, 2, 3))
    )

    placed = @at design (0.0, -1.0) connections = (phase = 1,)
    @test placed == at(design, 0.0, -1.0; connections = (phase = 1,))
    @test (@at design (0.0, -1.0) φ=0.2 connections=(phase = 1,)) ==
          at(design, 0.0, -1.0; φ = 0.2, connections = (phase = 1,))
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@at design (0.0, -1.0, 0.2) φ = 0.3)
    )

    root_space = terminal(
        :phase,
        core(copper; r = Grid((2e-3, 3e-3))),
        insulation(xlpe; t = 1e-3)
    )
    @test root_space isa Gridspace{Group}
    designs = build(CableDesign, "grid-cable", root_space)
    @test designs isa Gridspace{CableDesign}
    @test collect(designs) == [build(
               CableDesign,
               "grid-cable",
               terminal(:phase, core(copper; r), insulation(xlpe; t = 1e-3))
           ) for r in (2e-3, 3e-3)]
    @test rand(designs) isa CableDesign

    macro_designs = @cable "macro-grid" begin
        @terminal :phase begin
            core(copper; r = Grid((2e-3, 3e-3)))
            insulation(xlpe; t = 1e-3)
        end
    end
    @test macro_designs isa Gridspace{CableDesign}
    @test collect(macro_designs) == [build(
               CableDesign,
               "macro-grid",
               terminal(:phase, core(copper; r), insulation(xlpe; t = 1e-3))
           ) for r in (2e-3, 3e-3)]

    schedules = Grid((
        (LayRatio(10), LayRatio(11)),
        (LayRatio(12), LayRatio(13))
    ))
    strands = strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = schedules
    )
    @test strands isa Gridspace{Stack}
    @test length(strands) == 2
    @test strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = (LayRatio(10), LayRatio(11))
    ) isa Stack

    repeated_cores = cores(
        terminal(:core, core(copper; r = 0.5e-3));
        n = 2,
        r = Grid((2e-3, 3e-3)),
        names = (:a, :b)
    )
    @test repeated_cores isa Gridspace{Assembly}
    @test collect(repeated_cores) == [cores(
               terminal(:core, core(copper; r = 0.5e-3));
               n = 2,
               r,
               names = (:a, :b)
           ) for r in (2e-3, 3e-3)]

    radii = Grid((2e-3, 3e-3))
    advanced_spaces = (
        Pose2(radii, 0, 0),
        Ring(6; r = radii),
        Polar(nr = 1, nφ = 6, dr = radii),
        Fill(r = radii, φ = pi / 6),
        Lattice(nx = 2, ny = 1, dx = radii, dy = 0),
        LayRatio(radii),
        Pitch(radii),
        LayAngle(Grid((0.1, 0.2))),
        Helix(LayRatio(10); φ0 = Grid((0.0, 0.1))),
        FillFactor(Grid((0.8, 0.9))),
        DiameterFactor(Grid((0.8, 0.9))),
        TabulatedCompaction(Grid(((course = 1,), (course = 2,)))),
        AffineCompaction(Grid((([1.0 0.0; 0.0 1.0]), ([0.9 0.0; 0.0 1.0]))))
    )
    @test all(space -> space isa Gridspace, advanced_spaces)
    @test all(==(2), length.(advanced_spaces))
    @test collect(advanced_spaces[2]) == [Ring(6; r) for r in (2e-3, 3e-3)]
    @test collect(advanced_spaces[6]) == collect(LayRatio.((2e-3, 3e-3)))
    @test collect(advanced_spaces[9]) == [Helix(LayRatio(10); φ0) for φ0 in (0.0, 0.1)]

    local_member = @at terminal(:aux, core(copper; r = 1e-3)) (0.01, 0.02, 0.3)
    @test local_member == at(
        terminal(:aux, core(copper; r = 1e-3)),
        0.01,
        0.02;
        φ = 0.3
    )
    @test local_member.at == DM.Pose2(0.01, 0.02, 0.3)

    @test (@distribute strand(
        copper;
        wire = Disk(0.5e-3),
        layers = 1
    )).items[2].pattern.n == capacity()
    @test (@distribute rope(
        terminal(:aux, core(copper; r = 0.5e-3));
        layers = 1
    )).items[2].pattern.n == capacity()

    distributed_space = @distribute wires(
        copper;
        wire = Disk(Grid((0.25e-3, 0.5e-3))),
        r = 3e-3
    )
    @test distributed_space isa Gridspace{Group}
    resolved_counts = [length(DM.resolve(DM.EmptyBoundary(), group).regions)
                       for group in distributed_space]
    @test resolved_counts == [capacity(Ring, 3e-3, member_radius)
           for member_radius in (0.25e-3, 0.5e-3)]
end
