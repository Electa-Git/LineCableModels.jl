# cable 220kv milliken 1x2500 252 trefoil
# Run from the repository root:
#     julia --project=test test/manual/cases/cable_220kv_milliken_1x2500_252_trefoil.jl
# Results: .linecablemodels/fem/cable_220kv_milliken_1x2500_252_trefoil/yyyymmdd-HHMMSS/
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

const core_outer_radius = 0.03115
const core_wire_radius = 0.001475
const core_sectors = 6
const core_fillet_factor = 0.04
const binder_thickness = 0.0005
const inner_semicon_thickness = 0.0016
const insulation_thickness = 0.0203
const outer_semicon_thickness = 0.0013
const screen_inner_bedding_thickness = 0.0003
const screen_wires = 65
const screen_wire_radius = 0.00111
const screen_lay_ratio = 5.5012442232492
const screen_outer_bedding_thickness = 0.00073
const aluminum_foil_thickness = 0.0002
const jacket_thickness = 0.0048
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
const fem_frequency_indices = 1:10:101

# Construction only; verification stops before fem_problem.
problem = let
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    copper = LineCableModels.Material(materials, :copper)
    pe = LineCableModels.Material(materials, :pe)
    # HDPE and PE have identical electrical properties in the source library.
    hdpe = pe
    xlpe = LineCableModels.Material(
        kind = :insulator, rho = 1.97e14, eps_r = 2.5,
        tan_delta = xlpe_tan_delta
    )
    carbon_pe = LineCableModels.Material(
        kind = :semicon, rho = semicon_rho, eps_r = semicon_eps_r
    )

    span = 2pi / core_sectors
    fillet = core_fillet_factor * core_outer_radius
    segment = LineCableModels.Sector(
        span = span,
        r_base = fillet * (inv(sin(span / 2)) - 1),
        r_back = core_outer_radius,
        fillet = fillet
    )
    core = LineCableModels.terminal(
        :core,
        LineCableModels.milliken(
            aluminum;
            shape = LineCableModels.Disk(core_wire_radius),
            segment,
            segments = core_sectors,
            fill = pe
        ),
        LineCableModels.bedding(pe; t = binder_thickness, tag = :conductor_binder),
        LineCableModels.screen(carbon_pe; t = inner_semicon_thickness, tag = :inner_semicon),
        LineCableModels.insulation(xlpe; t = insulation_thickness, tag = :xlpe),
        LineCableModels.screen(carbon_pe; t = outer_semicon_thickness, tag = :outer_semicon),
        LineCableModels.bedding(pe; t = screen_inner_bedding_thickness, tag = :screen_bedding_inner)
    )
    screen_inner = core_outer_radius + binder_thickness +
                   inner_semicon_thickness + insulation_thickness +
                   outer_semicon_thickness + screen_inner_bedding_thickness
    screen_outer = screen_inner + 2screen_wire_radius
    wires = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(:copper_screen_wires, LineCableModels.Disk(screen_wire_radius), copper);
        pattern = LineCableModels.Ring(screen_wires; r = screen_inner +
                                                           screen_wire_radius),
        path = LineCableModels.Helix(LineCableModels.LayRatio(screen_lay_ratio))
    )
    metallic_screen = LineCableModels.Stack(
        LineCableModels.pipe(
            wires; shape = LineCableModels.Annulus(screen_inner, screen_outer), fill = pe
        ),
        LineCableModels.bedding(pe; t = screen_outer_bedding_thickness, tag = :screen_bedding_outer),
        LineCableModels.terminal(
            :foil,
            LineCableModels.sheath(aluminum; t = aluminum_foil_thickness, tag = :aluminium_foil),
            LineCableModels.jacket(hdpe; t = jacket_thickness, tag = :hdpe_oversheath)
        )
    )
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "220kv_eaxecew_1x2500_252",
        LineCableModels.Stack(core, metallic_screen)
    )
    spacing = (2 + formation_clearance_ratio) * LineCableModels.outer_radius(design)
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(formation_x, formation_y),
        spacing,
        connections = (core = (1, 4, 7), sheath = (2, 5, 8), foil = (3, 6, 9))
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem, placements;
        system_id = "cable_220kv_milliken_1x2500_252_trefoil",
        line_length = line_length
    )
    LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = temperature,
        earth_props = LineCableModels.homogeneous(rho = earth_rho, eps_r = earth_eps_r),
        frequencies = frequencies
    )
end
fem_problem = LineCableModels.Engine.LineParametersProblem(
    problem.system;
    temperature = problem.temperature,
    earth_props = problem.earth_props,
    frequencies = frequencies[fem_frequency_indices]
)

reductions = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
# Default coaxial-annulus shunt model.
analytical_formulation = Formulation(
    earth_impedance = :unified, earth_admittance = :unified, options = reductions
)
boundary_formulation = Formulation(
    earth_impedance = :unified, earth_admittance = :unified,
    shunt_model = :boundary, options = reductions
)
fem_formulation = Formulation(:fem; options = (; physics = :helmholtz, reductions...))
fem_options = (
    keep_run_directory = true,
    verbosity = (default = 0,),
    gmsh_verbosity = 0,
    getdp_verbosity = 0
)

analytical_wall_seconds = @elapsed analytical_result = compute(problem, analytical_formulation)
boundary_wall_seconds = @elapsed boundary_result = compute(problem, boundary_formulation)
fem_wall_seconds = @elapsed fem_result = compute(fem_problem, fem_formulation; options = fem_options)
observe(fem_result, LineCableModels.frequencies) == frequencies[fem_frequency_indices] ||
    error("FEM frequencies differ from the selected analytical grid points")

output_directory = joinpath(
    pwd(), ".linecablemodels", "fem", "cable_220kv_milliken_1x2500_252_trefoil", Dates.format(now(), "yyyymmdd-HHMMSS")
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
        gridpoint = (id = nothing, case = "cable_220kv_milliken_1x2500_252_trefoil", model = "fem"),
        length_unit = :base, clip = false),
    joinpath(output_directory, "fem.json"))
export_data(:json, problem; file_name = joinpath(output_directory, "problem.json"))

coordinates = details(analytical_result).data.coordinates
coordinates == details(fem_result).data.coordinates ||
    error("Analytical and FEM conductor row/column orders differ")
coordinates == details(boundary_result).data.coordinates ||
    error("Boundary and default analytical conductor orders differ")
open(joinpath(output_directory, "matrices.csv"), "w") do io
    println(io, "# Case: cable_220kv_milliken_1x2500_252_trefoil")
    println(io, "# Conductors in row/column order: ",
        join(("$i = $label" for (i, label) in enumerate(coordinates)), ", "))
    println(io, "# Units: frequency in Hz; Z in ohm/m and Y in S/m, per unit length")
    println(io, "# analytical: Unified earth, default shunt, 101 points; fem: Helmholtz, production defaults, 11 points")
    println(io, "# analytical_boundary: Unified earth, boundary shunt, same 101 points")
    println(io, "source,frequency_hz,quantity,row,column,real,imaginary")
    for (source, result) in (("analytical", analytical_result), ("analytical_boundary", boundary_result), ("fem", fem_result))
        for (k, frequency) in enumerate(observe(result, LineCableModels.frequencies))
            for (quantity, matrices) in (("Z", observe(result, Z)), ("Y", observe(result, Y)))
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
println(comparison_io, "Case: cable_220kv_milliken_1x2500_252_trefoil")
println(comparison_io, "FEM: Helmholtz, Gamma = 0; production defaults; full conductor matrices")
println(comparison_io, "Metric (%): 100 * max|FEM - analytical| / max|analytical diagonal|; no acceptance filtering")
@printf(comparison_io, "Analytical wall seconds: %.1f\n", analytical_wall_seconds)
@printf(comparison_io, "Boundary analytical wall seconds: %.1f\n", boundary_wall_seconds)
@printf(comparison_io, "FEM wall seconds: %.1f\n", fem_wall_seconds)
println(comparison_io, "FEM run directory: ", fem_run_directory)
println(comparison_io, "Primary reference: matrices.csv; fem.json omits the formulation record")
for (source, result) in (("analytical (default shunt)", analytical_result), ("analytical_boundary", boundary_result))
    println(comparison_io, "FEM vs ", source, " (Unified earth)")
    println(comparison_io, "f_Hz e_Z_% Z_worst(i,j) e_Y_% Y_worst(i,j)")
    for (k, analytical_index) in enumerate(fem_frequency_indices)
        errors = map((Z, Y)) do quantity
            reference = observe(result, quantity)[:, :, analytical_index]
            difference, worst = findmax(abs.(observe(fem_result, quantity)[:, :, k] - reference))
            (100 * difference / maximum(abs, diag(reference)), Tuple(worst))
        end
        @printf(comparison_io, "%.17g %.17g (%d,%d) %.17g (%d,%d)\n",
            frequencies[analytical_index], errors[1][1], errors[1][2]...,
            errors[2][1], errors[2][2]...)
    end
end
comparison = String(take!(comparison_io))
print(comparison)
write(joinpath(output_directory, "comparison.txt"), comparison)
println("Output directory: ", output_directory)
println("FEM run directory: ", fem_run_directory)
