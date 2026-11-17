# The problem of `test/references/cable_132kv_630mm2_flathor/problem.jl`: three 132 kV
# cables in a horizontal flat formation. The script computes the line parameters with the
# coaxial formulation at the 101 frequencies of the problem and with FEM at 11 of them. It
# writes both to `matrices.csv` in the format of `test/references/`.
#
# Run from the repository root:
#     `julia --project=docs test/manual/cases/cable_132kv_630mm2_flathor.jl`
# Results: `.linecablemodels/cases/cable_132kv_630mm2_flathor/<yyyymmdd-HHMMSS>/matrices.csv`

using LineCableModels, Gmsh, getdp_jll, Dates, LinearAlgebra, Printf

case = "cable_132kv_630mm2_flathor"
description = "Three 132 kV cables of 630 mm² in a horizontal flat formation 1 m deep, " *
              "the core, the wire screen with its tape and the aluminum tape each " *
              "a terminal"

problem = include(joinpath(@__DIR__, "..", "..", "references", case, "problem.jl"))
analytical = compute(problem,
    CoaxialFormulation(bundle = false, kron = false, transposition = false))
indices = 1:10:101

fem_problem = LineParametersProblem(problem.system, problem.earth;
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
