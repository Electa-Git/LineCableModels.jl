# The problem of `test/references/cable_380kv_2000mm2_flatver/problem.jl`: three 380 kV
# cables in a vertical flat formation. The script computes the line parameters with the
# coaxial formulation at the 101 frequencies of the problem and with FEM at 11 of them. It
# writes both to `matrices.csv` in the format of `test/references/`.
#
# Run from the repository root:
#     `julia --project=docs test/manual/cases/cable_380kv_2000mm2_flatver.jl`
# Results: `.linecablemodels/cases/cable_380kv_2000mm2_flatver/<yyyymmdd-HHMMSS>/matrices.csv`

using LineCableModels, Gmsh, getdp_jll, Dates, LinearAlgebra, Printf

case = "cable_380kv_2000mm2_flatver"
description = "Three 380 kV cables of 2000 mm² in a vertical flat formation from 1 m to " *
              "2 m deep, the core, the wire screen with its tape and the aluminum tape " *
              "each a terminal"

problem = include(joinpath(@__DIR__, "..", "..", "references", case, "problem.jl"))
analytical = compute(problem,
    CoaxialFormulation(bundle = false, kron = false, transposition = false))
indices = 1:10:101

fem_problem = LineParametersProblem(problem.system, problem.earth;
    temperature = problem.temperature,
    frequencies = problem.frequencies[indices])
# With `workers = 1`, the run meshes and solves one frequency at a time. Two at a time
# exceeded 16 GB of memory.
fem = compute(fem_problem, FEM(bundle = false, kron = false, transposition = false);
    keep_run_directory = true, workers = 1)

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
