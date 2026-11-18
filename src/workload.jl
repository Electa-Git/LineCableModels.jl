# The precompile workload. Julia compiles the methods that the cables below call and stores
# them with the package. The loop at the end flattens each cable and computes its line
# parameters at one frequency.
@setup_workload begin
    copper = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.0, 20.0, 0.00429)
    steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045)
    semicon = Material(:semicon, 1000.0, 1000.0)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    pvc = Material(:insulator, Inf, 8.0)
    air = Material(:insulator, Inf)
    earth = homogeneous(rho = 100.0)
    @compile_workload begin
        cables = Pair{<:CableDesign, <:NamedTuple}[]
        # A bare wire, then the wire with its insulation.
        wire = Terminal(:core, Region(:core, Circle(0.0425), copper))
        push!(cables, build(CableDesign, "bare", wire) => (core = 1,))
        push!(cables,
            build(CableDesign, "insulated", wire,
                Region(:insulation, Offset(1.0e-3), xlpe)) => (core = 1,))
        # A stranded core under a sheath and an armour of wires. The core takes round wires
        # first and rectangular ones second.
        for strand in (Circle(1.5e-3), Rectangle(3.0e-3, 3.0e-3))
            courses = [Lay(Helix(LayRatio(12)),
                           Repeat(CircularPattern(n; r = 3.0e-3k),
                               Region(:core_strands, strand, copper)))
                       for (k, n) in enumerate((6, 12, 18))]
            armour = Lay(Helix(LayRatio(10)),
                Repeat(CircularPattern(40; r = 26.5e-3),
                    Region(:armour_wires, Circle(2.0e-3), steel)))
            push!(cables,
                build(CableDesign, "stranded",
                    Terminal(:core,
                        Group(Region(:core_fill, Circle(10.5e-3), air),
                            Stack(Region(:core_centre, Circle(1.5e-3), copper), courses...))),
                    Region(:semicon_inner, Offset(1.0e-3), semicon),
                    Region(:insulation, Offset(8.0e-3), xlpe),
                    Region(:semicon_outer, Offset(1.0e-3), semicon),
                    Terminal(:sheath, Region(:sheath, Annulus(20.5e-3, 21.5e-3), copper)),
                    Region(:bedding, Offset(3.0e-3), xlpe),
                    Terminal(:armour, Group(Region(:armour_fill, Offset(4.0e-3), air), armour)),
                    Region(:jacket, Offset(4.0e-3), xlpe)) => (
                    core = 1, sheath = 2, armour = 3))
        end
        # A conductor of solid segments in a fill, under a screen of wires and a foil.
        segment = Region(:segment,
            Sector(0, 31.15e-3, π / 3, 0, 1.25e-3, Placement(0, 0)), aluminum)
        push!(cables,
            build(CableDesign, "segments",
                Terminal(:core,
                    Group(Region(:core_fill, Circle(31.15e-3), xlpe),
                        Repeat(CircularPattern(6; r = 0.0), segment))),
                Region(:semicon_inner, Offset(1.6e-3), semicon),
                Region(:insulation, Offset(20.3e-3), xlpe),
                Region(:semicon_outer, Offset(1.3e-3), semicon),
                Terminal(:sheath,
                    wires(copper; shape = Circle(1.11e-3), n = 65,
                        lay = LayRatio(5.5), fill = xlpe, tag = :screen)),
                Region(:bedding, Offset(0.73e-3), xlpe),
                Terminal(:foil, Region(:foil, Offset(0.2e-3), aluminum)),
                Region(:jacket, Offset(4.8e-3), xlpe)) => (core = 1, sheath = 2, foil = 3))
        # An insulated sector alone, then such sectors in a filler, inside a neutral of wires.
        t = 1.1e-3
        sector = Sector(0, 10.24e-3, deg2rad(119.0), t / sind(60), 1.02e-3, Placement(0, 0))
        push!(cables,
            build(CableDesign, "sector",
                Terminal(:core, Region(:core, sector, aluminum)),
                Region(:insulation, Offset(t), pvc),
                Region(:jacket, Offset(3.0e-3), pvc)) => (core = 1,))
        sectors = map((:a, :b, :c), (0.0, 2π / 3, 4π / 3)) do name, θ
            member = Stack(Terminal(name, Region(name, sector, aluminum)),
                Region(:insulation, Offset(t), pvc))
            return Placement(0.0, 0.0, θ) => member
        end
        push!(cables,
            build(CableDesign, "sectors",
                Enclosure(Region(:filler, Circle(13.57e-3), pvc), Assembly(sectors)),
                Terminal(:neutral,
                    Group(Region(:neutral_fill, Offset(1.58e-3), pvc),
                        Repeat(CircularPattern(30; r = Offset(0)),
                            Region(:neutral_wires, Circle(0.79e-3), copper)))),
                Region(:jacket, Offset(2.1e-3), pvc)) => (a = 1, b = 2, c = 3, neutral = 4))
        # Round cores in a steel pipe.
        at = Placement(2 * 28.0e-3 / sqrt(3), 0.0)
        cores = map((:a, :b, :c), (0.0, 2π / 3, 4π / 3)) do name, θ
            member = Stack(Terminal(name, Region(name, Circle(19.5e-3, at), aluminum)),
                Region(:insulation, Offset(8.5e-3), xlpe))
            return Placement(0.0, 0.0, θ) => member
        end
        push!(cables,
            build(CableDesign, "pipe",
                Enclosure(Region(:filler, Circle(62.0e-3), xlpe), Assembly(cores)),
                Terminal(:pipe, Region(:pipe, Offset(3.0e-3), steel)),
                Region(:jacket, Offset(5.0e-3), xlpe)) => (a = 1, b = 2, c = 3, pipe = 4))
        # A cable that the builders declare, with a screen of wires under a tape.
        r = 1.8e-3
        push!(cables,
            (@design "builders" begin
                @terminal :core begin
                    stranded(copper; center = Circle(r), shape = Circle(r), n = (6, 12),
                        lay = LayRatio(13.0), boundary = Circle(5r), tag = :core)
                end
                screen(semicon; t = 0.8e-3, tag = :semicon_inner)
                insulation(xlpe; t = 8.0e-3)
                screen(semicon; t = 0.8e-3, tag = :semicon_outer)
                @terminal :sheath begin
                    wires(copper; shape = Circle(0.6e-3), n = 40, lay = LayRatio(10.0),
                        tag = :screen)
                    tape(copper; w = 10.0e-3, t = 0.1e-3, lay = LayRatio(10.0))
                    sheath(aluminum; t = 0.3e-3)
                end
                bedding(xlpe; t = 0.5e-3)
                jacket(xlpe; t = 2.5e-3)
            end) => (core = 1, sheath = 2))
        for (design, connection) in cables
            flatten(design)
            system = build(CableSystem, [design], [Placement(0.0, -1.0)];
                connections = [connection], id = design.id, line_length = 1.0)
            compute(LineParametersProblem(system, earth; frequencies = [50.0]),
                CoaxialFormulation())
        end
    end
end
