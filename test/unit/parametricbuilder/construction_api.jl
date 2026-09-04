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
        section = Rectangle(1.0e-3, 0.5e-3),
        n = 12,
        lay = LayRatio(10)
    )
    @test_throws ArgumentError terminal(:insulating_tapes, insulating_tapes)
    with_tapes = build(CableDesign, "taped", phase, insulating_tapes)
    @test with_tapes.terminal_order == [:phase]
    @test count(iszero, with_tapes.terminal_map) > 0

    distributed_tapes = @distribute tape(
        xlpe;
        section = Rectangle(1.0e-3, 0.5e-3),
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
        Ring(capacity(); r = 8.75e-3, gap_frac = 0.02),
        Rectangle(1.0e-3, 0.5e-3),
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
        core(copper, Sector(span = pi / 3, r_base = 0.5e-3, r_back = 3e-3)),
        insulation(xlpe; t = 0.5e-3)
    )
    sector_b = terminal(
        :b,
        core(copper, Sector(span = pi / 4, r_base = 0.4e-3, r_back = 2.5e-3)),
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
    @test multicore.geometry.regions[1].primitive isa DM.SectorShape

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
        :build, :terminal, :core, :stranded, :rope, :cores, :tape,
        :insulation, :screen, :sheath, :armor, :bedding, :jacket, :filler,
        :pipe, :duct, :at, :trefoil, :hflat, :vflat, :capacity,
        :solid, :shell, :wires, :layers, :assembly, :layer, :homogeneous,
        :EarthLayer, :EarthModel,
        :Region, :Stack, :Group, :Assembly, :Enclosure, :Pose2,
        :Disk, :Rectangle, :Ellipse,
        :Sector, :Annulus, :Shell,
        :Polygon, :Ring, :Polar, :Fill, :Lattice,
        :Helix, :LayRatio, :Pitch, :LayAngle, :FillFactor,
        Symbol("@cable"), Symbol("@system"), Symbol("@earth"),
        Symbol("@terminal"), Symbol("@assembly"),
        Symbol("@pipe"), Symbol("@duct"), Symbol("@at"),
        Symbol("@hflat"), Symbol("@vflat"),
        Symbol("@trefoil"), Symbol("@distribute")
    )
    public_names = Set(names(LineCableModels))
    @test all(name -> name in public_names, required)

    forbidden = (
        :DuctBank, :Pipe, :Duct, :Core, :Cable, :Trefoil, :Earth,
        :RoundedSector, :Hexa, :DiameterFactor, :TabulatedCompaction,
        :AffineCompaction, :distribute, :strand,
        Symbol("@design"), Symbol("@core"),
        Symbol("@insulation"), Symbol("@strand"), Symbol("@rope"),
        Symbol("@armor"), Symbol("@screen"), Symbol("@shell"),
        Symbol("@wire"), Symbol("@fill")
    )
    @test all(name -> !(name in public_names), forbidden)
    @test Base.ispublic(DM, :BentStrip)
    @test Base.ispublic(DM, :BoundedPlacement)
    @test LineCableModels.build === DM.build === PB.build
    @test all(name -> name in public_names,
        (:Disk, :Rectangle, :Ellipse, :Sector, :Annulus, :Polygon))

    expansions = (
        :(@cable "x" nominal_data=metadata combine=:zip begin
            part
        end) => :build,
        :(@system "x" line_length=1 begin
            placed
        end) => :build,
        :(@earth begin
            layer
        end) => :build,
        :(@terminal :phase combine=:zip begin
            part
        end) => :terminal,
        :(@assembly combine=:zip begin
            left
            right
        end) => :assembly,
        :(@pipe shape=shape fill=fill combine=:zip begin
            left
            right
        end) => :pipe,
        :(@duct shape=shape fill=fill combine=:zip begin
            left
            right
        end) => :duct,
        :(@at subject (1, 2) φ=3) => :at,
        :(@trefoil design spacing=spacing center=(0, 0) combine=:zip phase=(1, 2, 3)) =>
            :trefoil,
        :(@hflat design spacing=spacing combine=:zip phase=(1, 2, 3)) => :hflat,
        :(@vflat design spacing=spacing combine=:zip phase=(1, 2, 3)) => :vflat
    )
    for (surface, callee) in expansions
        expanded = Base.remove_linenums!(macroexpand(@__MODULE__, surface))
        @test expanded isa Expr && expanded.head === :call
        @test expanded.args[1] == GlobalRef(PB, callee)
    end

    combined_surfaces = (
        :(@cable "x" combine=:zip begin
            part
        end),
        :(@terminal :phase combine=:zip begin
            part
        end),
        :(@assembly combine=:zip begin
            left
            right
        end),
        :(@pipe shape=shape fill=fill combine=:zip begin
            part
        end),
        :(@duct shape=shape fill=fill combine=:zip begin
            part
        end),
        :(@trefoil design spacing=spacing combine=:zip phase=(1, 2, 3)),
        :(@hflat design spacing=spacing combine=:zip phase=(1, 2, 3)),
        :(@vflat design spacing=spacing combine=:zip phase=(1, 2, 3))
    )
    for surface in combined_surfaces
        expanded = Base.remove_linenums!(macroexpand(@__MODULE__, surface))
        parameters = only(filter(
            argument -> argument isa Expr && argument.head === :parameters,
            expanded.args
        ))
        @test any(parameters.args) do keyword
            keyword isa Expr && keyword.head === :kw &&
                keyword.args == [:combine, QuoteNode(:zip)]
        end
    end

    distributed = Base.remove_linenums!(macroexpand(
        @__MODULE__,
        :(@distribute wires(material; shape = wire, r = radius))
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
    nested = @cable once(:identifier, "single-evaluation") nominal_data=(
        rated_voltage = once(:rated_voltage, 132e3),
    ) begin
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
        :rated_voltage => 1,
        :terminal => 1,
        :material => 1,
        :radius => 1
    )
    @test nested.nominal_data == (rated_voltage = 132e3,)

    system = @system once(:system_id, "single-system-evaluation") line_length=once(
        :line_length, 1.0
    ) begin
        @at nested (0.0, -1.0) phase=1
    end
    @test system isa LineCableSystem
    @test system.system_id == "single-system-evaluation"
    @test evaluations[:system_id] == 1
    @test evaluations[:line_length] == 1
end

@testitem "ParametricBuilder / construction API / earth grammar" tags=[:unit] begin
    using LineCableModels

    earth = @earth begin
        layer(rho = 100.0, eps_r = 10.0, thickness = 5.0)
        layer(rho = 500.0, eps_r = 20.0)
    end
    @test earth isa EarthModel{Float64}
    @test !earth.vertical_layers
    @test length(earth.layers) == 3
    @test isinf(first(earth.layers).rho)
    @test getproperty.(earth.layers[2:end], :rho) == (100.0, 500.0)
    @test getproperty.(earth.layers[2:end], :thickness) == (5.0, Inf)
    @test layer(rho = 100.0f0) isa EarthLayer{Float32}

    air = layer(rho = Inf, eps_r = 1.0006, mu_r = 1.0)
    explicit_air = @earth air_layer=air begin
        layer(rho = 100.0)
    end
    @test first(explicit_air.layers) === air

    earth_model = @inferred homogeneous(rho = 100.0, eps_r = 10.0)
    @test earth_model isa EarthModel{Float64}
    @test length(earth_model.layers) == 2
    @test earth_model.layers[2].rho == 100.0

    spaces = @earth begin
        layer(rho = Grid((10.0, 100.0)), eps_r = 5.0, thickness = 2.0)
        layer(rho = 500.0, eps_r = 20.0)
    end
    @test spaces isa Gridspace{EarthModel}
    @test length(spaces) == 2
    @test [model.layers[2].rho for model in spaces] == [10.0, 100.0]
    @test all(model -> model.layers[3].rho == 500.0, spaces)

    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@earth unsupported=true begin
            layer(rho = 100.0)
        end)
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
        shape = once(:wire_definition, Disk(0.2e-3)),
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

    piped = pipe(
        at(cell_a, -5e-3, 0),
        at(cell_b, 5e-3, 0);
        shape = Disk(15e-3),
        fill = air
    )
    notated_pipe = @pipe shape=Disk(15e-3) fill=air begin
        @at cell_a (-5e-3, 0)
        @at cell_b (5e-3, 0)
    end
    @test notated_pipe == piped

    pipe_space = @pipe shape=Disk(Grid((15e-3, 16e-3))) fill=Grid((
        air,
        concrete
    )) combine=:zip begin
        cell_a
    end
    @test pipe_space isa Gridspace{Enclosure}
    @test length(pipe_space) == 2

    repeated = duct(
        cell_a;
        formation = Lattice(nx = 2, ny = 1, dx = 12e-3, dy = 0),
        shape = Rectangle(24e-3, 12e-3),
        fill = concrete
    )
    @test repeated.item isa Assembly
    @test repeated.item.item === cell_a
    @test repeated.item.pattern == Lattice(nx = 2, ny = 1, dx = 12e-3, dy = 0)

    elliptical = duct(
        at(cell_a, -6e-3, 0),
        at(cell_b, 6e-3, 0);
        shape = Ellipse(15e-3, 8e-3),
        fill = concrete
    )
    elliptical_design = build(CableDesign, "elliptical-duct", elliptical)
    @test elliptical_design.geometry.outer isa Ellipse
    @test area(elliptical_design.geometry.outer) ≈ pi * 15e-3 * 8e-3
    elliptical_fill = only(filter(elliptical_design.geometry.regions) do region
        region.primitive isa DM.DifferenceShape &&
            region.primitive.outer isa Ellipse
    end)
    @test elliptical_fill.primitive isa DM.DifferenceShape
    @test length(elliptical_fill.primitive.holes) == 2
    @test all(hole -> hole isa Disk, elliptical_fill.primitive.holes)
    @test_throws DomainError build(
        CableDesign,
        "overfilled-elliptical-duct",
        duct(
            at(cell_a, -14e-3, 0),
            at(cell_b, 14e-3, 0);
            shape = Ellipse(15e-3, 8e-3),
            fill = concrete
        )
    )

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
    using Measurements: Measurement, uncertainty
    import LineCableModels.DataModel as DM

    copper = Material(kind = :conductor, rho = 1.72e-8)
    xlpe = Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)

    @test LayRatio(13, 12, 11) ==
          (LayRatio(13), LayRatio(12), LayRatio(11))
    @test LayRatio((13, 12, 11)) == LayRatio(13, 12, 11)
    @test Pitch(0.13, 0.12) == (Pitch(0.13), Pitch(0.12))
    @test LayAngle([0.13, 0.12]) == (LayAngle(0.13), LayAngle(0.12))
    @test FillFactor(0.88, 0.9) == (FillFactor(0.88), FillFactor(0.9))

    round_boundary = Disk(3.0e-3)
    round = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = 6,
        lay = (LayRatio(12), Pitch(0.2)),
        dir = (1, -1),
        φ0 = (0.0, 0.1),
        boundary = round_boundary
    )
    @test round isa Stack
    @test length(round.items) == 1
    round_group = only(round.items)
    @test round_group.boundary == round_boundary
    @test round_group.pattern === nothing
    round_courses = round_group.item.items
    @test round_courses[1].pattern === nothing
    @test round_courses[1].path === nothing
    @test [item.pattern.n for item in round_courses[2:3]] == [6, 12]
    @test getfield.(round_courses[2:3], :path) == [
        Helix(LayRatio(12); dir = 1, φ0 = 0.0),
        Helix(Pitch(0.2); dir = -1, φ0 = 0.1)
    ]

    compact_boundary = Disk(sqrt(19 / 0.9) * 0.5e-3)
    normalized = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = LayRatio(12, 11),
        compact = FillFactor(0.9),
        boundary = compact_boundary
    )
    normalized_group = only(normalized.items)
    @test normalized_group.compact == FillFactor(0.9)
    @test getfield.(normalized_group.item.items[2:3], :path) == [
        Helix(LayRatio(12)),
        Helix(LayRatio(11); φ0 = pi / 6)
    ]

    sector_boundary = Sector(
        span = deg2rad(119.0),
        r_base = 1.10e-3,
        r_back = 10.24e-3,
        fillet = 1.02e-3
    )
    resolved_sector = DM.resolve(DM.EmptyBoundary(), sector_boundary)
    sector_fill = 18pi * (0.5e-3)^2 / DM.area(resolved_sector)
    bounded = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = (LayRatio(12), LayRatio(11)),
        compact = FillFactor(sector_fill),
        boundary = sector_boundary
    )
    @test bounded isa Stack
    @test length(bounded.items) == 1
    bounded_group = only(bounded.items)
    @test bounded_group isa Group
    @test bounded_group.pattern === nothing
    @test bounded_group.path === nothing
    @test bounded_group.compact == FillFactor(sector_fill)
    @test bounded_group.boundary == sector_boundary
    @test bounded_group.item isa Stack
    @test length(bounded_group.item.items) == 2
    @test getproperty.(bounded_group.item.items, :path) == [
        Helix(LayRatio(12)),
        Helix(LayRatio(11); φ0 = pi / 6)
    ]
    bounded_design = build(
        CableDesign,
        "bounded-stranded-core",
        terminal(:core, bounded)
    )
    @test length(bounded_design.geometry.regions) == 18
    @test DM.boundary(bounded_design.geometry).primitive == sector_boundary
    @test all(region -> region.primitive isa DM.Polygon,
        bounded_design.geometry.regions)
    @test sum(DM.area, bounded_design.geometry.regions) ≈ 18pi * (0.5e-3)^2
    sector_centre = DM.centroid(resolved_sector)
    course_counts = [first(region.placement.patterns).pattern.n
                     for region in bounded_design.geometry.regions]
    path_radii = [only(region.paths).radius for region in bounded_design.geometry.regions]
    @test all(zip(bounded_design.geometry.regions, path_radii)) do (region, radius)
        centre = DM.centroid(region.primitive)
        radius ≈ hypot(
            centre[1] - sector_centre[1],
            centre[2] - sector_centre[2]
        )
    end
    @test maximum(path_radii[course_counts .== 6]) <
          minimum(path_radii[course_counts .== 12])
    bounded_component = only(DM.flatten(bounded_design, 50.0))
    @test bounded_component.conductor.r_ex ≈ sqrt(DM.area(resolved_sector) / pi)
    @test all(isapprox.(
        bounded_component.conductor.position,
        DM.centroid(resolved_sector)
    ))
    member_resistances = map(bounded_design.geometry.regions) do region
        path = only(region.paths)
        copper.rho / DM.area(region.source.primitive) *
        overlength(path.path, path.radius)
    end
    @test bounded_component.conductor.resistance ≈
          inv(sum(inv, member_resistances))

    natural_sector = stranded(
        copper;
        shape = Disk(0.5e-3),
        boundary = sector_boundary
    )
    natural_group = only(natural_sector.items)
    @test length(natural_group.item.items) == 1
    @test only(natural_group.item.items).pattern.n == capacity()
    natural_sector_design = build(
        CableDesign,
        "natural-sector-strands",
        terminal(:core, natural_sector)
    )
    @test all(region -> region.primitive isa Disk,
        natural_sector_design.geometry.regions)
    @test length(natural_sector_design.geometry.regions) ==
          length(DM.sector_wire_sites(resolved_sector, Disk(0.5e-3)))
    natural_centres = DM.centroid.(getproperty.(
        natural_sector_design.geometry.regions,
        :primitive
    ))
    @test all(
        hypot(
            natural_centres[left][1] - natural_centres[right][1],
            natural_centres[left][2] - natural_centres[right][2]
        ) + 1e-12 >= 1e-3
    for left in 1:(length(natural_centres) - 1)
    for right in (left + 1):length(natural_centres)
    )
    @test_throws ArgumentError stranded(
        copper;
        center = Disk(0.5e-3),
        shape = Disk(0.5e-3),
        layers = 1,
        n = 6,
        boundary = sector_boundary
    )

    second_boundary = Sector(
        span = deg2rad(118.0),
        r_base = 1.0e-3,
        r_back = 10.0e-3,
        fillet = 0.9e-3
    )
    boundary_space = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 1,
        n = 6,
        boundary = Grid((sector_boundary, second_boundary))
    )
    @test boundary_space isa Gridspace{Stack}
    @test collect(boundary_space) == [stranded(
               copper;
               shape = Disk(0.5e-3),
               layers = 1,
               n = 6,
               boundary
           ) for boundary in (sector_boundary, second_boundary)]

    product_space = stranded(
        copper;
        shape = Disk(Grid((0.4e-3, 0.5e-3))),
        layers = 1,
        n = 6,
        dir = Grid((1, -1)),
        boundary = Disk(2.0e-3)
    )
    @test product_space isa Gridspace{Stack}
    @test length(product_space) == 4

    radii = (0.4e-3, 0.5e-3)
    counts = (6, 7)
    fills = (0.9, 0.8)
    zipped_space = stranded(
        copper;
        shape = Grid(Disk.(radii)),
        layers = 1,
        n = Grid(counts),
        lay = Grid((LayRatio(12), LayRatio(11))),
        dir = Grid((1, -1)),
        compact = Grid(FillFactor.(fills)),
        boundary = Grid(Disk.(sqrt.((1 .+ counts) ./ fills) .* radii)),
        combine = :zip
    )
    zipped_designs = [build(CableDesign, "zipped-bounded-$index", terminal(:phase, value))
                      for (index, value) in enumerate(zipped_space)]
    @test length(zipped_designs) == 2
    @test [length(design.geometry.regions) for design in zipped_designs] == [7, 8]

    uncertain_space = stranded(
        copper;
        shape = Disk(Grid(0.5e-3, 1.0)),
        layers = 1,
        n = 6,
        compact = FillFactor(0.9),
        boundary = Disk(sqrt(7 / 0.9) * 0.5e-3)
    )
    uncertain_design = build(
        CableDesign,
        "uncertain-bounded-core",
        terminal(:phase, first(uncertain_space))
    )
    @test eltype(uncertain_design) <: Measurement
    @test any(
        region -> uncertainty(DM.area(region.primitive)) > 0,
        uncertain_design.geometry.regions
    )

    @test_throws ArgumentError stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 0,
        boundary = "rounded"
    )

    @test_throws ArgumentError stranded(
        copper;
        center = Disk(0.5e-3),
        shape = Rectangle(0.6e-3, 0.8e-3),
        layers = 2,
        n = (5, 11),
        lay = nothing,
        compact = FillFactor(1),
        boundary = Disk(sqrt((0.5e-3)^2 + 16 * 0.6e-3 * 0.8e-3 / pi))
    )

    compacted = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 1,
        n = 6,
        compact = FillFactor(0.9),
        lay = LayRatio(11),
        boundary = Disk(sqrt(7 / 0.9) * 0.5e-3)
    )
    compacted_design = build(
        CableDesign,
        "compacted",
        terminal(:phase, compacted)
    )
    deformed = filter(
        region -> region.primitive isa DM.Polygon,
        compacted_design.geometry.regions
    )
    @test length(deformed) == 7
    @test sum(DM.area, deformed) ≈ 7pi * (0.5e-3)^2
    outer_course = filter(region -> length(region.paths) == 1, deformed)
    @test length(outer_course) == 6
    @test all(
        region -> only(region.paths).radius ≈ hypot(DM.centroid(region)...),
        outer_course
    )

    saturated_boundary = Disk(sqrt(7) * 0.5e-3)
    saturated = build(
        CableDesign,
        "saturated-compaction",
        terminal(
            :phase,
            stranded(
                copper;
                shape = Disk(0.5e-3),
                layers = 1,
                n = 6,
                compact = FillFactor(1),
                boundary = saturated_boundary
            )
        )
    )
    @test sum(DM.area, saturated.geometry.regions) ≈ DM.area(saturated_boundary)
    @test all(
        region -> DM.support(region.primitive) <=
                  saturated_boundary.r * (1 + 5e-6),
        saturated.geometry.regions
    )
    @test_throws DomainError build(
        CableDesign,
        "inconsistent-fill-factor",
        terminal(
            :phase,
            stranded(
                copper;
                shape = Disk(0.5e-3),
                layers = 1,
                n = 6,
                compact = FillFactor(0.5),
                boundary = saturated_boundary
            )
        )
    )

    @test_throws ArgumentError stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = (6, capacity()),
        boundary = compact_boundary
    )

    direct_capacity = capacity(Ring, 10e-3, 0.5e-3; gap_frac = 0.03)
    pattern = Ring(capacity(); r = 10e-3, gap_frac = 0.03)
    @test capacity(pattern, Disk(0.5e-3), nothing) == direct_capacity
    distributed = @distribute wires(
        copper;
        shape = Disk(0.5e-3),
        r = 10e-3,
        gap_frac = 0.03
    )
    @test distributed.pattern.n == capacity()
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@distribute wires(copper; shape = Disk(1e-3), n = 6, r = 2e-3))
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
    scheduled_rope = rope(
        round;
        layers = 2,
        n = (6, 12),
        lay = LayRatio(10, 9)
    )
    @test getfield.(scheduled_rope.items[2:3], :path) == [
        Helix(LayRatio(10)),
        Helix(LayRatio(9))
    ]
    nested_design = build(
        CableDesign,
        "nested-rope",
        terminal(:phase, nested)
    )
    @test any(region -> length(region.paths) == 2, nested_design.geometry.regions)
    @test any(region -> isempty(region.paths), nested_design.geometry.regions)
    leaf_resistances = map(nested_design.geometry.regions) do region
        primitive = region.source.primitive
        base = DM.tubular_resistance(
            0.0,
            primitive.r,
            region.source.material.rho
        )
        foldl(region.paths; init = base) do resistance, entry
            resistance * overlength(entry.path, entry.radius)
        end
    end
    expected_resistance = inv(sum(inv, leaf_resistances))
    component = only(LineCableModels.DataModel.flatten(
        nested_design,
        50.0
    ))
    @test component.conductor.resistance ≈ expected_resistance

    @test_throws DimensionMismatch stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        lay = (LayRatio(10),),
        boundary = round_boundary
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
    @test (@at design (0.0, -1.0) phase=1) == placed
    @test (@at design (0.0, -1.0) φ=0.2 connections=(phase = 1,)) ==
          at(design, 0.0, -1.0; φ = 0.2, connections = (phase = 1,))
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@at design (0.0, -1.0, 0.2) φ = 0.3)
    )
    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@at design (0.0, -1.0) connections=(phase = 1,) phase=2)
    )

    system = @system "placed-system" line_length=2.0 begin
        @at design (0.0, -1.0) phase=1
        @at design (0.1, -1.0) phase=2
    end
    expected_system = build(
        LineCableSystem,
        [
            at(design, 0.0, -1.0; connections = (phase = 1,)),
            at(design, 0.1, -1.0; connections = (phase = 2,))
        ];
        system_id = "placed-system",
        line_length = 2.0
    )
    @test system.system_id == expected_system.system_id
    @test system.line_length == expected_system.line_length
    @test system.designs == expected_system.designs
    @test system.positions == expected_system.positions
    @test system.connections == expected_system.connections
    @test system.terminal_order == expected_system.terminal_order
    @test system.terminal_map == expected_system.terminal_map
    @test system.connection_order == expected_system.connection_order

    formation_system = @system "formation-system" begin
        @hflat design spacing=0.1 phase=(1, 2, 3)
    end
    @test length(formation_system.designs) == 3
    @test formation_system.connection_order == [1, 2, 3]

    @test_throws ArgumentError macroexpand(
        @__MODULE__,
        :(@system "invalid" system_id="duplicate" begin
            placed
        end)
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

    zipped_terminal = @terminal :phase combine=:zip begin
        core(copper; r = Grid((2e-3, 3e-3)))
        insulation(xlpe; t = Grid((1e-3, 2e-3)))
    end
    @test zipped_terminal isa Gridspace{Group}
    @test length(zipped_terminal) == 2

    zipped_designs = @cable Grid(("zip-a", "zip-b")) nominal_data=Grid((
        (rating = 1,),
        (rating = 2,)
    )) combine=:zip begin
        zipped_terminal
    end
    @test zipped_designs isa Gridspace{CableDesign}
    @test length(zipped_designs) == 2
    @test getproperty.(collect(zipped_designs), :cable_id) == ["zip-a", "zip-b"]
    @test getproperty.(collect(zipped_designs), :nominal_data) == [
        (rating = 1,),
        (rating = 2,)
    ]

    left_members = terminal(:left, core(copper; r = Grid((1e-3, 1.2e-3))))
    right_members = terminal(:right, core(copper; r = Grid((0.8e-3, 0.9e-3))))
    zipped_assembly = @assembly combine=:zip begin
        left_members
        right_members
    end
    @test zipped_assembly isa Gridspace{Assembly}
    @test length(zipped_assembly) == 2

    zipped_formation = @trefoil first(zipped_designs) spacing=Grid((
        0.08,
        0.1
    )) center=at(x = Grid((0.0, 0.2)), y = -1.0) combine=:zip phase=(1, 2, 3)
    @test zipped_formation isa Gridspace{Vector}
    @test length(zipped_formation) == 2

    systems = @system "macro-system" begin
        @at macro_designs (0.0, -1.0) phase=1
    end
    @test systems isa Gridspace{LineCableSystem}
    @test length(systems) == 2
    @test all(system -> system.system_id == "macro-system", systems)
    @test all(system -> system.connection_order == [1], systems)

    schedules = Grid((
        LayRatio(10, 11),
        LayRatio(12, 13)
    ))
    strands = stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = schedules,
        boundary = Disk(3e-3)
    )
    @test strands isa Gridspace{Stack}
    @test length(strands) == 2
    @test stranded(
        copper;
        shape = Disk(0.5e-3),
        layers = 2,
        n = (6, 12),
        lay = LayRatio(10, 11),
        boundary = Disk(3e-3)
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
        FillFactor(Grid((0.8, 0.9)))
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

    @test (@distribute rope(
        terminal(:aux, core(copper; r = 0.5e-3));
        layers = 1
    )).items[2].pattern.n == capacity()

    distributed_space = @distribute wires(
        copper;
        shape = Disk(Grid((0.25e-3, 0.5e-3))),
        r = 3e-3
    )
    @test distributed_space isa Gridspace{Group}
    resolved_counts = [length(DM.resolve(DM.EmptyBoundary(), group).regions)
                       for group in distributed_space]
    @test resolved_counts == [capacity(Ring, 3e-3, member_radius)
           for member_radius in (0.25e-3, 0.5e-3)]
end
