# Reference order: first the buried bare and insulated wire pairs and the two coaxial
# cases, then the overhead and mixed bare-wire pairs. Cables with sector-shaped conductors
# follow, then a pipe-type cable. The first two entries keep their positions.
const REFERENCES = let
    copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
    # The wire cases take the copper of the materials library that their scripts read.
    library_copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393)
    insulation = Material(:insulator, 1.97e14, 2.3)
    polymer = Material(:insulator, 1.0e14, 2.3)
    semicon = Material(:semicon, 1000.0, 1000.0)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    air = material(kind = :insulator, rho = Inf)
    frequencies = 10.0 .^ range(-1, 7; length = 101)
    metal = Terminal(:core, Region(:core_metal, Circle(0.0425), library_copper))
    wires = [
        "two_bare_wires" => build(CableDesign, "bare_wire", metal),
        "two_insulated_wires" => build(CableDesign, "insulated_wire", metal,
            Region(:core_insulation, Offset(1.0e-3), insulation))]
    cases = [case => LineParametersProblem(
                 build(CableSystem, [design, design],
                     [Placement(0.0, -1.0), Placement(1.0, -1.0)];
                     connections = [(core = 1,), (core = 2,)], id = case, line_length = 1.0),
                 homogeneous(rho = 0.1); frequencies)
             for (case, design) in wires]
    for (case, strand) in (("coaxial_circular_strands", Circle(1.5e-3)),
        ("coaxial_rectangular_strands", Rectangle(3.0e-3, 3.0e-3)))
        courses = [Lay(Helix(LayRatio(12)),
                       Repeat(CircularPattern(n; r = 3.0e-3k),
                           Region(Symbol(:core_strands_, k), strand, copper)))
                   for (k, n) in enumerate((6, 12, 18))]
        armour = Lay(Helix(LayRatio(10)),
            Repeat(CircularPattern(40; r = 26.5e-3),
                Region(:armour_wires, Circle(2.0e-3), steel)))
        design = build(CableDesign, "coaxial",
            Terminal(:core,
                Enclosure(Region(:core_fill, Circle(10.5e-3), air),
                    Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
            Region(:core_semicon_inner, Offset(1.0e-3), semicon),
            Region(:core_insulation, Offset(8.0e-3), polymer),
            Region(:core_semicon_outer, Offset(1.0e-3), semicon),
            Terminal(:sheath, Region(:sheath_copper, Annulus(20.5e-3, 21.5e-3), copper)),
            Region(:sheath_bedding, Offset(3.0e-3), polymer),
            Terminal(:armour, Enclosure(Region(:armour_fill, Offset(4.0e-3), air), armour)),
            Region(:armour_jacket, Offset(4.0e-3), polymer))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, sheath = 2, armour = 3)], id = case,
            line_length = 1000.0)
        push!(cases,
            case => LineParametersProblem(system, homogeneous(rho = 100, eps_r = 10);
                frequencies))
    end
    design = build(CableDesign, "bare_wire", metal)
    for (case, placements) in (
        ("overhead_bare_wires", [Placement(0.0, 1.0), Placement(1.0, 1.0)]),
        ("mixed_bare_wires", [Placement(0.0, 1.0), Placement(1.0, -1.0)]))
        system = build(CableSystem, [design, design], placements;
            connections = [(core = 1,), (core = 2,)], id = case, line_length = 1.0)
        push!(cases,
            case => LineParametersProblem(system, homogeneous(rho = 0.1); frequencies))
    end
    aluminum = Material(:conductor, 2.8264e-8, 1, 1, 20, 0.00429)
    pvc = Material(:insulator, Inf, 8.0)
    # The sector of a 1 kV cable of three cores. Its sides meet 1.1 mm / sin(60°) from its
    # center, where the insulations of three such cores touch.
    span, t = deg2rad(119.0), 1.1e-3
    sector = LineCableModels.Sector(
        0, 10.24e-3, span, t / sind(60), 1.02e-3, Placement(0, 0))
    single = build(CableDesign, "insulated_sector",
        Terminal(:core, Region(:core, sector, aluminum)),
        Region(:core_insulation, Offset(t), pvc),
        Region(:jacket, Offset(3.0e-3), pvc))
    system = build(CableSystem, [single], [Placement(0.0, -1.0)];
        connections = [(core = 1,)], id = "one_insulated_sector", line_length = 1.0)
    push!(cases,
        "one_insulated_sector" => LineParametersProblem(system,
            homogeneous(rho = 100, eps_r = 10); frequencies))
    pe = material(kind = :insulator, rho = 1.0e15, eps_r = 2.3)
    xlpe = material(kind = :insulator, rho = 1.97e14, eps_r = 2.5, tan_delta = 1.0e-3)
    semiconductive = Material(:semicon, 0.06, 1000.0)
    # A 220 kV single core whose conductor is six solid aluminum segments in a PE fill.
    segment = Region(:segment,
        LineCableModels.Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, Placement(0, 0)), aluminum)
    milliken = build(CableDesign, "milliken_220kv",
        Terminal(:core,
            Enclosure(Region(:core_fill, Circle(31.15e-3), pe),
                Repeat(CircularPattern(6; r = 0.0), segment))),
        Region(:binder, Offset(0.5e-3), pe),
        Region(:semicon_inner, Offset(1.6e-3), semiconductive),
        Region(:xlpe, Offset(20.3e-3), xlpe),
        Region(:semicon_outer, Offset(1.3e-3), semiconductive),
        Region(:bedding_inner, Offset(0.3e-3), pe),
        # The package's builder `wires`, since the local `wires` above is a list of designs.
        Terminal(:sheath,
            LineCableModels.wires(copper; shape = Circle(1.11e-3), n = 65,
                lay = LayRatio(5.5012442232492), fill = pe, tag = :screen)),
        Region(:bedding_outer, Offset(0.73e-3), pe),
        Terminal(:foil, Region(:foil, Offset(0.2e-3), aluminum)),
        Region(:jacket, Offset(4.8e-3), pe))
    # A 1 kV cable of three such sectors, each with its insulation, and a neutral of 30 wires.
    function insulated(name, θ)
        member = Stack(Terminal(name, Region(name, sector, aluminum)),
            Region(Symbol(name, :_insulation), Offset(t), pvc))
        return Placement(0.0, 0.0, θ) => member
    end
    r_wire, r_center, r_outer = 0.79e-3, 14.36e-3, 17.25e-3
    three = build(CableDesign, "NAYCWY_3x95_30x2_5",
        Enclosure(Region(:filler, Circle(r_center - r_wire), pvc),
            Assembly((insulated(:a, 0.0), insulated(:b, 2π / 3),
                insulated(:c, 4π / 3)))),
        Terminal(:neutral,
            Enclosure(Region(:neutral_fill, Offset(2r_wire), pvc),
                Repeat(CircularPattern(30; r = Offset(0)),
                    Region(:neutral_wires, Circle(r_wire), copper)))),
        Region(:jacket, Offset(r_outer - r_center - r_wire), pvc))
    # A 36 kV pipe-type cable of three round cores, after Høidalen 2013, Table 1. Each core
    # lies 2 × 28 mm / √3 from the center of the pipe, so that the insulations touch.
    pipe_aluminum = Material(:conductor, 3.365e-8, 1, 1, 20, 0.00429)
    pipe_steel = Material(:conductor, 7.0e-8, 1, 10, 20, 0.0045)
    offset = 2 * 28.0e-3 / sqrt(3)
    function round_core(name, θ)
        conductor = Region(name, Circle(19.5e-3, Placement(offset, 0.0)), pipe_aluminum)
        member = Stack(Terminal(name, conductor),
            Region(Symbol(name, :_insulation), Offset(8.5e-3), xlpe))
        return Placement(0.0, 0.0, θ) => member
    end
    pipe = build(CableDesign, "pipe_type_3x1000",
        Enclosure(Region(:filler, Circle(62.0e-3), xlpe),
            Assembly((round_core(:a, 0.0), round_core(:b, 2π / 3),
                round_core(:c, 4π / 3)))),
        Terminal(:pipe, Region(:pipe, Offset(3.0e-3), pipe_steel)),
        Region(:jacket, Offset(5.0e-3), pe))
    for (case, design, connections) in (
        ("milliken_solid_segments", milliken, [(core = 1, sheath = 2, foil = 3)]),
        ("three_insulated_sectors", three, [(a = 1, b = 2, c = 3, neutral = 4)]),
        ("three_cores_steel_pipe", pipe, [(a = 1, b = 2, c = 3, pipe = 4)]))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)]; connections,
            id = case, line_length = 1.0)
        push!(cases,
            case => LineParametersProblem(system, homogeneous(rho = 100, eps_r = 10);
                frequencies))
    end
    cases
end
