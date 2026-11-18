# Three 18/30 kV NA2XS2Y cables with a stranded aluminum core of 630 mm² in a touching
# trefoil formation 1 m deep. Each cable has two terminals: the core and the screen of
# copper wires with its copper tape. The problem has 101 frequencies from 0.1 Hz to 10 MHz.
let
    aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429)
    copper = Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
    xlpe = Material(:insulator, 1.97e14, 2.5)
    pe = Material(:insulator, 1.97e14, 2.3)
    semicon1 = Material(:semicon, 1000.0, 1000.0)
    semicon2 = Material(:semicon, 500.0, 1000.0)
    air = Material(:insulator, Inf)

    r = 3.65e-3 / 2 # radius of a core wire [m]
    r_screen = 0.85e-3 / 2 # radius of a screen wire [m]
    design = @design "na2xs2y_1x630_35_18_30kv" begin
        @terminal :core begin
            # The central wire is a strand at the radius 0. The courses start at the angle 0,
            # and the even ones half a pitch further.
            Group(Region(:core_fill, Circle(9r), air),
                Stack(
                    Repeat(CircularPattern(1; r = 0.0), Region(:core, Circle(r), aluminum)),
                    (Lay(Helix(LayRatio(15.0)),
                         Repeat(CircularPattern(n; r = Offset(0), φ0 = φ0),
                             Region(:core, Circle(r), aluminum)))
                    for (n, φ0) in zip((6, 12, 18, 24), (0.0, π / 12, 0.0, π / 24)))...))
        end
        screen(semicon1; t = 0.5e-3, tag = :inner_semiconductor)
        insulation(xlpe; t = 8.0e-3)
        screen(semicon2; t = 0.4e-3, tag = :outer_semiconductor)
        @terminal :sheath begin
            Group(Region(:screen_matrix, Offset(2r_screen + 0.1e-3), air),
                Stack(
                    Lay(Helix(LayRatio(15.0)),
                        Repeat(CircularPattern(60; r = Offset(0)),
                            Region(:sheath_wires, Circle(r_screen), copper))),
                    Lay(Helix(LayRatio(15.0)),
                        Repeat(CircularPattern(1; r = Offset(0)),
                            Region(:sheath_copper_tape, Rectangle(15.0e-3, 0.1e-3), copper)))))
        end
        jacket(pe; t = 2.7e-3)
    end

    placements = trefoil(design; spacing = 2r_ex(design.outer), at = Placement(0.0, -1.0))
    system = build(CableSystem, fill(design, 3), placements;
        connections = [(core = 2k - 1, sheath = 2k) for k in 1:3], id = "cable_30kv_na2xs2y_630mm2_trefoil",
        line_length = 1000.0)
    LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
        temperature = 20.0, frequencies = 10.0 .^ range(-1, 7; length = 101))
end
