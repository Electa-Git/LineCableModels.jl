# Three 18 kV cables with a stranded aluminum core of 1000 mm² in a trefoil formation 1 m
# deep. The core and the aluminum tape under the jacket are terminals of each cable. So is
# the screen of copper wires with its copper tape. The problem has 101 frequencies from
# 0.1 Hz to 10 MHz.
let
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429)
    copper = Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    pe = Material(:insulator, 1.97e14, 2.3)
    semicon1 = Material(:semicon, 1000.0, 1000.0)
    semicon2 = Material(:semicon, 500.0, 1000.0)
    polyacrylate = Material(:semicon, 5.3e3, 32.3)
    matrix = Material(:insulator, Inf)

    r = 4.19e-3 / 2 # radius of a core wire [m]
    design = @design "18kV_1000mm2" begin
        @terminal :core begin
            # The central wire is a strand at the radius 0. The courses start at the angle 0,
            # and the even ones half a pitch further.
            Group(Region(:core_fill, Circle(9r), matrix),
                Stack(
                    Repeat(CircularPattern(1; r = 0.0), Region(:core, Circle(r), aluminum)),
                    (Lay(Helix(LayRatio(q)),
                         Repeat(CircularPattern(n; r = Offset(0), φ0 = φ0),
                             Region(:core, Circle(r), aluminum)))
                    for (n, q, φ0) in zip((6, 12, 18, 24), (15.0, 13.5, 12.5, 11.0),
                        (0.0, π / 12, 0.0, π / 24)))...))
        end
        screen(polyacrylate; t = 0.3e-3, tag = :core_semicon_tape_inner)
        screen(semicon1; t = 0.768e-3, tag = :core_semicon_inner)
        insulation(xlpe; t = 8.3e-3, tag = :core_insulation)
        screen(semicon2; t = 0.472e-3, tag = :core_semicon_outer)
        screen(polyacrylate; t = 0.3e-3, tag = :core_semicon_tape_outer)
        @terminal :sheath begin
            Group(Region(:screen_matrix, Offset(0.94e-3 + 0.1e-3), matrix),
                Stack(
                    Lay(Helix(LayRatio(10.0)),
                        Repeat(CircularPattern(49; r = Offset(0)),
                            Region(:sheath_wires, Circle(0.47e-3), copper))),
                    Lay(Helix(LayRatio(10.0)),
                        Repeat(CircularPattern(1; r = Offset(0)),
                            Region(:sheath_copper_tape, Rectangle(10.0e-3, 0.1e-3), copper)))))
        end
        screen(polyacrylate; t = 0.94e-3, tag = :sheath_water_blocking)
        @terminal :jacket begin
            sheath(aluminum; t = 0.15e-3, tag = :jacket_aluminum_tape)
        end
        jacket(pe; t = 5.0e-5, tag = :jacket_pe_face)
        jacket(pe; t = 3.4e-3, tag = :jacket_insulation)
    end

    placements = trefoil(design; spacing = 2.2r_ex(design.outer), at = Placement(0.0, -1.0))
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = 3k - 2, sheath = 3k - 1, jacket = 3k) for k in 1:3], id = "cable_18kv_1000mm2_trefoil",
        line_length = 1000.0)
    LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
        temperature = 20.0, frequencies = 10.0 .^ range(-1, 7; length = 101))
end
