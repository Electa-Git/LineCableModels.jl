# cable 525kv 1600mm2 bipole
# Run from the repository root:
#     julia --project=test test/manual/cases/cable_525kv_1600mm2_bipole.jl
# Results: .linecablemodels/fem/cable_525kv_1600mm2_bipole/yyyymmdd-HHMMSS/
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

const strand_layers = 6
const strands_per_layer = 6
const strand_diameter = 0.0036649
const core_lay_ratio = 11.0
const inner_semicon_thickness = 0.002
const insulation_thickness = 0.026
const outer_semicon_thickness = 0.0018
const water_blocking_thickness = 0.0003
const lead_screen_thickness = 0.0033
const inner_sheath_thickness = 0.003
const bedding_thickness = 0.003
const armor_wires = 68
const armor_wire_diameter = 0.005827
const armor_lay_ratio = 10.0
const armor_packing_clearance_ratio = 0.0
const jacket_thickness = 0.01
const cable_x = (-0.5, 0.5)
const cable_y = -1.0
const line_length = 1000.0
const temperature = 20.0
const earth_rho = 100.0
const earth_eps_r = 10.0
const frequencies = 10.0 .^ range(-1, 7; length = 101)
const fem_frequency_indices = 1:10:101

# Construction only; verification stops before fem_problem.
problem = let
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    copper = LineCableModels.Material(materials, :copper)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    pe = LineCableModels.Material(materials, :pe)
    polyacrylate = LineCableModels.Material(materials, :polyacrylate)
    lead = LineCableModels.Material(materials, :lead)
    pp = LineCableModels.Material(materials, :pp)
    steel = LineCableModels.Material(materials, :steel)
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )
    radius_before_bedding = (strand_layers + 0.5) * strand_diameter +
                            inner_semicon_thickness +
                            insulation_thickness +
                            outer_semicon_thickness +
                            water_blocking_thickness +
                            lead_screen_thickness +
                            inner_sheath_thickness
    armor_wire_radius = armor_wire_diameter / 2
    unbuffered_armor_radius = radius_before_bedding + bedding_thickness
    unbuffered_outer_radius = unbuffered_armor_radius + armor_wire_diameter +
                              jacket_thickness
    packing_buffer = armor_packing_clearance_ratio * unbuffered_outer_radius
    # A nonzero ratio is an explicit design override, never uncertainty repair.
    effective_bedding_thickness = bedding_thickness + packing_buffer
    strand_radius = strand_diameter / 2
    radius = (2strand_layers + 1) * strand_radius
    stranded_core = LineCableModels.stranded(
        copper;
        shape = LineCableModels.Disk(strand_radius),
        lay = LineCableModels.LayRatio(core_lay_ratio),
        boundary = LineCableModels.Disk(radius),
        fill = matrix
    )
    core = stranded_core
    parts = LineCableModels.AbstractCablePart[core]
    for (tag, thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon1),
        (:core_insulation, insulation_thickness, pe),
        (:core_semicon_outer, outer_semicon_thickness, semicon2),
        (:core_water_blocking, water_blocking_thickness, polyacrylate)
    )
        push!(parts, LineCableModels.Region(
            tag, LineCableModels.Shell(thickness), material
        ))
        radius += thickness
    end

    lead_outer = radius + lead_screen_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_lead_screen,
                LineCableModels.Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_inner,
            LineCableModels.Shell(inner_sheath_thickness),
            pe
        ))
    radius += inner_sheath_thickness
    push!(parts,
        LineCableModels.Region(
            :sheath_bedding,
            LineCableModels.Shell(effective_bedding_thickness),
            pp
        ))
    radius += effective_bedding_thickness

    armor_outer = radius + 2armor_wire_radius
    armor_centre = radius + armor_wire_radius
    armor = LineCableModels.Group(
        :armor,
        LineCableModels.Region(
            :armor_wires, LineCableModels.Disk(armor_wire_radius), steel
        );
        pattern = LineCableModels.Ring(armor_wires; r = armor_centre),
        path = LineCableModels.Helix(LineCableModels.LayRatio(armor_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :armor_matrix,
            armor;
            primitive = LineCableModels.Annulus(radius, armor_outer),
            fill = matrix
        ))
    push!(parts, LineCableModels.Region(
        :armor_jacket,
        LineCableModels.Shell(jacket_thickness),
        pp
    ))
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "525kV_1600mm2",
        LineCableModels.Stack(parts)
    )
    earth = LineCableModels.homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:2]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 2),
        [LineCableModels.Pose2(cable_x[index], cable_y) for index in 1:2];
        connections,
        system_id = "cable_525kv_1600mm2_bipole",
        line_length = line_length
    )
    LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = temperature,
        earth_props = earth,
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
fem_formulation = Formulation(:fem; options = (; physics = :helmholtz, reductions...))
fem_options = (
    frequency_workers = 8,
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
    pwd(), ".linecablemodels", "fem", "cable_525kv_1600mm2_bipole", Dates.format(now(), "yyyymmdd-HHMMSS")
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
        gridpoint = (id = nothing, case = "cable_525kv_1600mm2_bipole", model = "fem"),
        length_unit = :base, clip = false),
    joinpath(output_directory, "fem.json"))
export_data(:json, problem; file_name = joinpath(output_directory, "problem.json"))

coordinates = details(analytical_result).data.coordinates
coordinates == details(fem_result).data.coordinates ||
    error("Analytical and FEM conductor row/column orders differ")
open(joinpath(output_directory, "matrices.csv"), "w") do io
    println(io, "# Case: cable_525kv_1600mm2_bipole")
    println(io, "# Conductors in row/column order: ",
        join(("$i = $label" for (i, label) in enumerate(coordinates)), ", "))
    println(io, "# Units: frequency in Hz; Z in ohm/m and Y in S/m, per unit length")
    println(io, "# analytical: Unified earth, default shunt, 101 points; fem: Helmholtz, production defaults, 11 points")
    println(io, "source,frequency_hz,quantity,row,column,real,imaginary")
    for (source, result) in (("analytical", analytical_result), ("fem", fem_result))
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
println(comparison_io, "Case: cable_525kv_1600mm2_bipole")
println(comparison_io, "FEM: Helmholtz, Gamma = 0; production defaults; full conductor matrices")
println(comparison_io, "Metric (%): 100 * max|FEM - analytical| / max|analytical diagonal|; no acceptance filtering")
@printf(comparison_io, "Analytical wall seconds: %.1f\n", analytical_wall_seconds)
@printf(comparison_io, "FEM wall seconds: %.1f\n", fem_wall_seconds)
println(comparison_io, "FEM run directory: ", fem_run_directory)
println(comparison_io, "Primary reference: matrices.csv; fem.json omits the formulation record")
for (source, result) in (("analytical (default shunt)", analytical_result),)
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
