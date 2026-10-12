# CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil
# Adapted from cigre_tb880_case0_132kv_630cu.json, CableBuilderSpec v1,
# and its neighbouring SeraingThermal/specs/materials.json.
# Preserve the declared solid-core radius, 0.01515 m, rather than rescaling it
# to the nominal 630 mm² name: its represented metal area is about 721.06 mm².
# Layer dimensions are literal source values, in metres. The trefoil, earth,
# temperature and frequency grid follow existing gauntlet operating defaults.
# This is a case derived from the supplied design, not a CIGRE reference result.
# Run from the repository root:
#     julia --project=test test/manual/cases/cable_132kv_cigre_tb880_case0_630cu_trefoil.jl
# Results: .linecablemodels/fem/cable_132kv_cigre_tb880_case0_630cu_trefoil/yyyymmdd-HHMMSS/
# Files: analytical.json, fem.json (no formulation record), problem.json,
#        matrices.csv (primary reference) and comparison.txt.

using LineCableModels
using Gmsh
using LinearAlgebra
using Printf
using Dates

# Including a script from another checkout does not switch the loaded package.
script_checkout = normpath(joinpath(@__DIR__, "../../.."))
realpath(pkgdir(LineCableModels)) == realpath(script_checkout) || error(
    "This script belongs to $script_checkout, but LineCableModels was loaded from " *
    "$(pkgdir(LineCableModels)). Restart Julia with --project=$script_checkout/test.")
println("LineCableModels source: ", pkgdir(LineCableModels))

const core_radius = 0.01515
const inner_semicon_thickness = 0.0015
const insulation_thickness = 0.0155
const outer_semicon_thickness = 0.0013
const sheath_thickness = 0.0008
const jacket_thickness = 0.0035
const semicon_rho = 0.06
const semicon_eps_r = 1000.0
const xlpe_tan_delta = 0.001
const formation_x = 0.0
const formation_y = -1.0
const formation_clearance_ratio = 0.2
const line_length = 1000.0
const temperature = 20.0
const earth_rho = 100.0
const earth_eps_r = 10.0
const frequencies = 10.0 .^ range(-1, 7; length = 101)
# 11 points, one per 0.8 decade, 0.1 Hz to 10 MHz; all shared with the analytical grid.
const fem_frequency_indices = 1:10:101

materials = LineCableModels.MaterialsLibrary(add_defaults = true)
copper = LineCableModels.Material(materials, :copper)
aluminum = LineCableModels.Material(materials, :aluminum)
hdpe = LineCableModels.Material(materials, :pe)
xlpe = LineCableModels.Material(
    kind = :insulator, rho = 1.97e14, eps_r = 2.5,
    tan_delta = xlpe_tan_delta
)
carbon_pe = LineCableModels.Material(
    kind = :semicon, rho = semicon_rho, eps_r = semicon_eps_r
)
design = LineCableModels.build(
    LineCableModels.CableDesign,
    "cigre_tb880_case0_132kv_630cu",
    LineCableModels.Stack(
        LineCableModels.terminal(
            :core,
            LineCableModels.core(copper; r = core_radius, tag = :compacted_copper_conductor),
            LineCableModels.screen(carbon_pe; t = inner_semicon_thickness, tag = :inner_semiconductor),
            LineCableModels.insulation(xlpe; t = insulation_thickness, tag = :xlpe_insulation),
            LineCableModels.screen(carbon_pe; t = outer_semicon_thickness, tag = :outer_semiconductor)
        ),
        LineCableModels.terminal(
            :sheath,
            LineCableModels.sheath(aluminum; t = sheath_thickness, tag = :aluminum_sheath),
            LineCableModels.jacket(hdpe; t = jacket_thickness, tag = :hdpe_oversheath)
        )
    )
)
spacing = (2 + formation_clearance_ratio) * LineCableModels.outer_radius(design)
placements = LineCableModels.trefoil(
    design;
    center = LineCableModels.at(formation_x, formation_y),
    spacing,
    connections = (core = (1, 3, 5), sheath = (2, 4, 6))
)
system = LineCableModels.build(
    LineCableModels.LineCableSystem, placements;
    system_id = "cable_132kv_cigre_tb880_case0_630cu_trefoil",
    line_length = line_length
)
problem = LineCableModels.Engine.LineParametersProblem(
    system;
    temperature = temperature,
    earth_props = LineCableModels.homogeneous(rho = earth_rho, eps_r = earth_eps_r),
    frequencies = frequencies
)
earth = problem.earth_props
fem_problem = LineCableModels.Engine.LineParametersProblem(
    system;
    temperature = temperature,
    earth_props = earth,
    frequencies = frequencies[fem_frequency_indices]
)

reductions = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
# Default coaxial-annulus shunt model.
analytical_formulation = Formulation(
    earth_impedance = :unified, earth_admittance = :unified, options = reductions
)
fem_formulation = Formulation(:fem; options = (; physics = :helmholtz, reductions...))
fem_options = (
    keep_run_directory = true,
    verbosity = (default = 0,),
    gmsh_verbosity = 0,
    getdp_verbosity = 0
)

analytical_wall_seconds = @elapsed analytical_result = compute(problem, analytical_formulation)
fem_wall_seconds = @elapsed fem_result = compute(fem_problem, fem_formulation; options = fem_options)
observe(fem_result, LineCableModels.frequencies) == frequencies[fem_frequency_indices] ||
    error("FEM frequencies differ from the selected analytical grid points")

output_directory = joinpath(
    pwd(), ".linecablemodels", "fem", "cable_132kv_cigre_tb880_case0_630cu_trefoil", Dates.format(now(), "yyyymmdd-HHMMSS")
)
mkpath(output_directory)
save(ObservedResult(analytical_result, (Z, Y); length_unit = :base, clip = false),
    joinpath(output_directory, "analytical.json"))
# The FEM formulation type cannot be archived, so this archive keeps the matrices,
# frequencies and conductor labels without the formulation record.
fem_matrices = LineParameters(PhaseDomain(), observe(fem_result, Z), observe(fem_result, Y),
    observe(fem_result, LineCableModels.frequencies); basis = :pul,
    details = ComputationDetails(coordinates = details(fem_result).data.coordinates))
save(ObservedResult(fem_matrices, (Z, Y);
        gridpoint = (id = nothing, case = "cable_132kv_cigre_tb880_case0_630cu_trefoil", model = "fem"),
        length_unit = :base, clip = false),
    joinpath(output_directory, "fem.json"))
export_data(:json, problem; file_name = joinpath(output_directory, "problem.json"))

const mu0 = 4pi * 1e-7
const eps0 = 8.8541878128e-12
Z_analytical = observe(analytical_result, Z)
Y_analytical = observe(analytical_result, Y)
Z_fem = observe(fem_result, Z)
Y_fem = observe(fem_result, Y)
coordinates = details(analytical_result).data.coordinates
coordinates == details(fem_result).data.coordinates ||
    error("Analytical and FEM conductor row/column orders differ")
open(joinpath(output_directory, "matrices.csv"), "w") do io
    println(io, "# Case: CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil (cable_132kv_cigre_tb880_case0_630cu_trefoil)")
    println(io, "# Conductors in row/column order: ",
        join(("$i = $label" for (i, label) in enumerate(coordinates)), ", "))
    println(io, "# Units: frequency in Hz; Z in ohm/m and Y in S/m, per unit length")
    println(io, "# Sources: analytical = Unified earth, default shunt model, 101 frequencies")
    println(io, "#          fem = Formulation(:fem), production defaults, 11 frequencies")
    println(io, "source,frequency_hz,quantity,row,column,real,imaginary")
    for (source, result) in (("analytical", analytical_result), ("fem", fem_result))
        source_matrices = (("Z", observe(result, Z)), ("Y", observe(result, Y)))
        for (k, frequency) in enumerate(observe(result, LineCableModels.frequencies))
            for (quantity, matrices) in source_matrices
                for row in axes(matrices, 1), column in axes(matrices, 2)
                    value = matrices[row, column, k]
                    @printf(io, "%s,%.17g,%s,%d,%d,%.17g,%.17g\n",
                        source, frequency, quantity, row, column, real(value), imag(value))
                end
            end
        end
    end
end
fem_run_directory = details(fem_result).data.fem.run.run_directory
comparison_io = IOBuffer()
println(comparison_io, "Case: CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil (cable_132kv_cigre_tb880_case0_630cu_trefoil)")
println(comparison_io, "Analytical formulation: Unified earth, default shunt model")
println(comparison_io,
    "FEM settings: Helmholtz, Gamma = 0; production defaults; full conductor matrices")
@printf(comparison_io, "Analytical wall seconds: %.1f\n", analytical_wall_seconds)
@printf(comparison_io, "FEM wall seconds: %.1f\n", fem_wall_seconds)
println(comparison_io, "FEM run directory: ", fem_run_directory)
println(comparison_io, "Acceptance level: 1 %")
println(comparison_io, "Primary reference: matrices.csv (frequency in Hz; Z in ohm/m; Y in S/m)")
println(comparison_io,
    "Archive note: fem.json has matrices, frequencies and labels, no formulation record")
@printf(comparison_io, "%10s %9s %12s %9s %12s %13s %11s\n",
    "f_Hz", "e_Z_%", "Z_worst(i,j)", "e_Y_%", "Y_worst(i,j)", "G_sign_canary", "kappa_r_p")
for (k, analytical_index) in enumerate(fem_frequency_indices)
    Z_ana = Z_analytical[:, :, analytical_index]
    Y_ana = Y_analytical[:, :, analytical_index]
    Z_difference, Z_worst = findmax(abs.(Z_fem[:, :, k] - Z_ana))
    Y_difference, Y_worst = findmax(abs.(Y_fem[:, :, k] - Y_ana))
    e_Z = 100 * Z_difference / maximum(abs, diag(Z_ana))
    e_Y = 100 * Y_difference / maximum(abs, diag(Y_ana))
    conductance_threshold = 1e-6 * maximum(abs, diag(Y_ana))
    G_sign_canary = count(
        (abs.(real.(Y_ana)) .>= conductance_threshold) .&
        (sign.(real.(Y_fem[:, :, k])) .!= sign.(real.(Y_ana)))
    )
    w = 2pi * frequencies[analytical_index]
    kappa = sqrt(im * w * mu0 * (1 / earth_rho + im * w * eps0 * earth_eps_r))
    kappa_r_p = abs(kappa * LineCableModels.outer_radius(design))
    @printf(comparison_io, "%10.4g %9.3f %12s %9.3f %12s %13d %11.3g\n",
        frequencies[analytical_index], e_Z, @sprintf("(%d,%d)", Tuple(Z_worst)...),
        e_Y, @sprintf("(%d,%d)", Tuple(Y_worst)...), G_sign_canary, kappa_r_p)
end
comparison = String(take!(comparison_io))
print(comparison)
write(joinpath(output_directory, "comparison.txt"), comparison)

println("Output directory: ", output_directory)
println("FEM run directory: ", fem_run_directory)
