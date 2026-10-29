# Reference order: the buried bare and insulated wire pairs, the two coaxial cases,
# then the overhead and mixed bare-wire pairs. The first two entries keep their positions.
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
    cases
end
