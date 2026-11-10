# The cables of the example `examples/cable_30kv_630mm2_trefoil.jl`, computed with FEM at
# the 50 frequencies of the example, with the screens eliminated and the phases transposed
# as there. The script writes the coaxial and FEM matrices to `matrices.csv` in the format
# of `test/references/`, and prints the deviations of the positive-sequence parameters.
#
# Run from the repository root:
#     `julia --project=docs test/manual/cases/cable_30kv_630mm2_trefoil.jl`
# Results: `matrices.csv` in
#     `.linecablemodels/cases/cable_30kv_630mm2_trefoil/<yyyymmdd-HHMMSS>/`

using LineCableModels, Gmsh, getdp_jll, Dates, LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "..", "..", "examples", "cable_30kv_630mm2_trefoil.jl"))
case = "cable_30kv_630mm2_trefoil"
description = "Three 18/30 kV cables of 630 mm² in a touching trefoil, the lower two 1 m " *
              "deep, the screens grounded and eliminated, the phases transposed"

fem = compute(problem, FEM(transposition = true); keep_run_directory = true)

folder = joinpath(pwd(), ".linecablemodels", "cases", case,
    Dates.format(now(), "yyyymmdd-HHMMSS"))
mkpath(folder)
open(joinpath(folder, "matrices.csv"), "w") do io
    println(io, "# Case: ", description, " (", case, ")")
    println(io, "# Conductors in row/column order: ",
        join(("$i = $label" for (i, label) in enumerate(result.terminals)), ", "))
    println(io, "# Units: frequency in Hz; Z in ohm/m and Y in S/m, per unit length")
    println(io, "# Sources: analytical = CoaxialFormulation(transposition = true), ",
        "50 frequencies")
    println(io, "#          fem = FEM(transposition = true), 50 frequencies")
    println(io, "source,frequency_hz,quantity,row,column,real,imaginary")
    for (source, parameters) in (("analytical", result), ("fem", fem)),
        k in eachindex(parameters.frequencies),
        (quantity, matrices) in (("Z", parameters.Z), ("Y", parameters.Y)),
        row in axes(matrices, 1), column in axes(matrices, 2)
        value = matrices[row, column, k]
        @printf(io, "%s,%.17g,%s,%d,%d,%.17g,%.17g\n", source, parameters.frequencies[k],
            quantity, row, column, real(value), imag(value))
    end
end

# FEM against the coaxial formulation. The comparison scales the matrices by their largest
# diagonal entry, and takes the positive-sequence resistance, reactance and susceptance
# relative to the coaxial.
z1_fem = fem.Z[1, 1, :] .- fem.Z[1, 2, :]
y1_fem = fem.Y[1, 1, :] .- fem.Y[1, 2, :]
for k in eachindex(result.frequencies)
    e_Z = maximum(abs, fem.Z[:, :, k] - result.Z[:, :, k]) /
          maximum(abs, diag(result.Z[:, :, k]))
    e_Y = maximum(abs, fem.Y[:, :, k] - result.Y[:, :, k]) /
          maximum(abs, diag(result.Y[:, :, k]))
    @printf("%7.1f Hz  e_Z %7.4f %%  e_Y %7.4f %%  R1 %+8.4f %%  X1 %+8.4f %%  B1 %+8.4f %%\n",
        result.frequencies[k], 100e_Z, 100e_Y,
        100(real(z1_fem[k]) / real(z1[k]) - 1), 100(imag(z1_fem[k]) / imag(z1[k]) - 1),
        100(imag(y1_fem[k]) / imag(y1[k]) - 1))
end
println("Results: ", folder)
