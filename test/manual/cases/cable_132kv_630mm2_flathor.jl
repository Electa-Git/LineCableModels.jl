# Three 132 kV cables with a stranded copper core of 630 mm² in a horizontal flat formation
# 1 m deep, 0.5 m apart. The core and the aluminum tape under the jacket are terminals of
# each cable, and so is the screen of copper wires with its copper tape. The script computes
# the line parameters with the coaxial formulation at 101 frequencies from 0.1 Hz to 10 MHz
# and with FEM at 11 of them. It writes both to `matrices.csv` in the format of
# `test/references/`.
#
# Run from the repository root:
#     `julia --project=docs test/manual/cases/cable_132kv_630mm2_flathor.jl`
# Results: `.linecablemodels/cases/cable_132kv_630mm2_flathor/<yyyymmdd-HHMMSS>/matrices.csv`

using LineCableModels, Gmsh, getdp_jll, Dates, LinearAlgebra, Printf

case = "cable_132kv_630mm2_flathor"
description = "Three 132 kV cables of 630 mm² in a horizontal flat formation 1 m deep, " *
              "the core, the wire screen with its tape and the aluminum tape each " *
              "a terminal"

aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429)
copper = Material(:conductor, 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393)
xlpe = Material(:insulator, 1.97e14, 2.5)
pe = Material(:insulator, 1.97e14, 2.3)
semicon1 = Material(:semicon, 1000.0, 1000.0)
semicon2 = Material(:semicon, 500.0, 1000.0)
polyacrylate = Material(:semicon, 5.3e3, 32.3)
matrix = Material(:insulator, Inf)

r = 0.0036648500338805653 / 2 # radius of a core wire [m]
r_screen = 0.0025881867280128637 / 2 # radius of a screen wire [m]
design = @design "132kV_630mm2" begin
    @terminal :core begin
        # The central wire is a strand at the radius 0. The courses start at the angle 0,
        # and the even ones half a pitch further.
        Enclosure(Region(:core_fill, Circle(9r), matrix),
            Stack(Repeat(CircularPattern(1; r = 0.0), Region(:core, Circle(r), copper)),
                (Lay(Helix(LayRatio(13.0)),
                     Repeat(CircularPattern(n; r = Offset(0), φ0 = φ0),
                         Region(:core, Circle(r), copper)))
                for (n, φ0) in zip((6, 12, 18, 24), (0.0, π / 12, 0.0, π / 24)))...))
    end
    screen(polyacrylate; t = 0.3e-3, tag = :core_semicon_tape_inner)
    screen(semicon1; t = 0.768e-3, tag = :core_semicon_inner)
    insulation(xlpe; t = 16.0e-3, tag = :core_insulation)
    screen(semicon2; t = 0.472e-3, tag = :core_semicon_outer)
    screen(polyacrylate; t = 0.3e-3, tag = :core_semicon_tape_outer)
    @terminal :sheath begin
        Enclosure(Region(:screen_matrix, Offset(2r_screen + 0.1e-3), matrix),
            Stack(
                Lay(Helix(LayRatio(10.0)),
                    Repeat(CircularPattern(19; r = Offset(0)),
                        Region(:sheath_wires, Circle(r_screen), copper))),
                tape(copper; w = 10.0e-3, t = 0.1e-3, lay = LayRatio(10.0),
                    tag = :sheath_copper_tape)))
    end
    screen(polyacrylate; t = 0.94e-3, tag = :sheath_water_blocking)
    @terminal :jacket begin
        sheath(aluminum; t = 0.15e-3, tag = :jacket_aluminum_tape)
    end
    jacket(pe; t = 5.05e-3, tag = :jacket_insulation)
end

placements = flat(design; n = 3, spacing = 0.5, at = Placement(0.0, -1.0))
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = 3k - 2, sheath = 3k - 1, jacket = 3k) for k in 1:3], id = case,
    line_length = 1000.0)
problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
    temperature = 20.0, frequencies = 10.0 .^ range(-1, 7; length = 101))
analytical = compute(problem,
    CoaxialFormulation(bundle = false, kron = false, transposition = false))
indices = 1:10:101

fem_problem = LineParametersProblem(system, problem.earth;
    temperature = problem.temperature,
    frequencies = problem.frequencies[indices])
# With `workers = 2`, the run meshes and solves two frequencies at a time.
fem = compute(fem_problem, FEM(bundle = false, kron = false, transposition = false);
    keep_run_directory = true, workers = 2)

folder = joinpath(pwd(), ".linecablemodels", "cases", case,
    Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)
open(joinpath(folder, "matrices.csv"), "w") do io
    println(io, "# Case: ", description, " (", case, ")")
    println(io, "# Conductors in row/column order: ",
        join(("$i = $label" for (i, label) in enumerate(analytical.terminals)), ", "))
    println(io, "# Units: frequency in Hz; Z in ohm/m and Y in S/m, per unit length")
    println(io, "# Sources: analytical = CoaxialFormulation(bundle = false, ",
        "kron = false, transposition = false), 101 frequencies")
    println(io, "#          fem = FEM(bundle = false, kron = false, ",
        "transposition = false), 11 frequencies")
    println(io, "source,frequency_hz,quantity,row,column,real,imaginary")
    for (source, result) in (("analytical", analytical), ("fem", fem)),
        k in eachindex(result.frequencies),
        (quantity, matrices) in (("Z", result.Z), ("Y", result.Y)),
        row in axes(matrices, 1), column in axes(matrices, 2)
        value = matrices[row, column, k]
        @printf(io, "%s,%.17g,%s,%d,%d,%.17g,%.17g\n", source, result.frequencies[k],
            quantity, row, column, real(value), imag(value))
    end
end

# FEM against the coaxial formulation, scaled by the largest diagonal entry.
for (k, j) in enumerate(indices)
    e_Z = maximum(abs, fem.Z[:, :, k] - analytical.Z[:, :, j]) /
          maximum(abs, diag(analytical.Z[:, :, j]))
    e_Y = maximum(abs, fem.Y[:, :, k] - analytical.Y[:, :, j]) /
          maximum(abs, diag(analytical.Y[:, :, j]))
    @printf("%11.4g Hz  e_Z %8.4f %%  e_Y %8.4f %%\n",
        fem.frequencies[k], 100e_Z, 100e_Y)
end
println("Results: ", folder)
