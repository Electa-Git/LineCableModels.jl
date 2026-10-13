# cable 380kv armoured ac flat
# Run from the repository root:
#     julia --project=test test/manual/cases/cable_380kv_armoured_ac_flat.jl
# Results: .linecablemodels/fem/cable_380kv_armoured_ac_flat/yyyymmdd-HHMMSS/
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
const strand_diameter = 0.00366
const core_lay_ratio = 15.0
const inner_semicon_thickness = 0.002
const insulation_thickness = 0.02
const outer_semicon_thickness = 0.0015
const water_blocking_thickness = 0.0003
const lead_sheath_thickness = 0.0033
const inner_pe_thickness = 0.003
const bedding_thickness = 0.009376520171274612
const armor_wires = 68
const armor_wire_diameter = 0.00582
const armor_lay_ratio = 15.0
const jacket_thickness = 0.01
const core_rho = 2.3853e-8
const core_mu_r = 1.01663
const xlpe_rho = 1.0e14
const xlpe_eps_r = 2.77781
const xlpe_mu_r = 1.59891
const lead_rho = 2.14e-7
const lead_mu_r = 0.999983
const pe_rho = 1.97e14
const pe_eps_r = 2.51862
const pp_rho = 1.0e15
const pp_eps_r = 2.8
const pp_mu_r = 1.1263
const steel_rho = 1.7475e-7
const steel_mu_r = 36.6326
const cable_x = (-1.0, 0.0, 1.0)
const cable_y = -1.0
const line_length = 20000.0
const voltage = 380.0
const temperature = 70.0
const earth_rho = 1000.0
const earth_eps_r = 10.0
const frequencies = 10.0 .^ range(-1, 7; length = 101)
const fem_frequency_indices = 1:10:101

# Construction only; verification stops before fem_problem.
problem = let
    core = LineCableModels.Material(
        :conductor, core_rho, 1.0, core_mu_r, 20.0, 0.00393
    )
    semicon_inner = LineCableModels.Material(:semicon, 1_000.0, 1_000.0)
    semicon_outer = LineCableModels.Material(:semicon, 500.0, 1_000.0)
    xlpe = LineCableModels.Material(
        :insulator, xlpe_rho, xlpe_eps_r, xlpe_mu_r
    )
    lead = LineCableModels.Material(
        :conductor, lead_rho, 1.0, lead_mu_r, 20.0, 0.004
    )
    pe = LineCableModels.Material(:insulator, pe_rho, pe_eps_r)
    pp = LineCableModels.Material(
        :insulator, pp_rho, pp_eps_r, pp_mu_r
    )
    steel = LineCableModels.Material(
        :conductor, steel_rho, 1.0, steel_mu_r, 20.0, 0.0045
    )
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    strand_radius = strand_diameter / 2
    strand_count = 1 + 3strand_layers * (strand_layers + 1)
    radius = sqrt(strand_count) * strand_radius
    packed_core = LineCableModels.stranded(
        core;
        shape = LineCableModels.Disk(strand_radius),
        lay = LineCableModels.LayRatio(core_lay_ratio),
        compact = true,
        boundary = LineCableModels.Disk(radius),
        fill = LineCableModels.Material(kind = :insulator, rho = Inf)
    )
    parts = LineCableModels.AbstractCablePart[packed_core]
    for (tag, layer_thickness, material) in (
        (:core_semicon_inner, inner_semicon_thickness, semicon_inner),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon_outer),
        (:core_water_blocking, water_blocking_thickness, semicon_outer)
    )
        push!(parts,
            LineCableModels.Region(tag, LineCableModels.Shell(layer_thickness), material))
        radius += layer_thickness
    end
    lead_outer = radius + lead_sheath_thickness
    push!(parts,
        LineCableModels.Group(
            :sheath,
            LineCableModels.Region(
                :sheath_lead,
                LineCableModels.Annulus(radius, lead_outer),
                lead
            )
        ))
    radius = lead_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_inner_pe,
            LineCableModels.Shell(inner_pe_thickness),
            pe
        ))
    radius += inner_pe_thickness

    armor_wire_radius = armor_wire_diameter / 2
    push!(parts,
        LineCableModels.Region(
            :sheath_bedding,
            LineCableModels.Shell(bedding_thickness),
            pp
        ))
    radius += bedding_thickness
    armor_outer = radius + 2armor_wire_radius
    armor = LineCableModels.Group(
        :armor,
        LineCableModels.Region(
            :armor_wires,
            LineCableModels.Disk(armor_wire_radius),
            steel
        );
        pattern = LineCableModels.Ring(
            armor_wires; r = radius + armor_wire_radius
        ),
        path = LineCableModels.Helix(LineCableModels.LayRatio(armor_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :armor_matrix,
            armor;
            primitive = LineCableModels.Annulus(radius, armor_outer),
            fill = matrix
        ))
    radius = armor_outer
    push!(parts,
        LineCableModels.Region(
            :armor_jacket,
            LineCableModels.Shell(jacket_thickness),
            pp
        ))

    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "cable_380kv_armoured",
        LineCableModels.Stack(parts);
        nominal_data = (
            designation_code = "cable_380kv_armoured",
            U0 = voltage / sqrt(3),
            U = voltage,
            conductor_cross_section = 1_600.0
        )
    )
    earth = LineCableModels.homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    connections = [Dict(:core => 3index - 2, :sheath => 3index - 1, :armor => 3index)
                   for index in 1:3]
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        fill(design, 3),
        [LineCableModels.Pose2(cable_x[index], cable_y) for index in 1:3];
        connections,
        environment = earth,
        system_id = "cable_380kv_armoured_ac_flat",
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
    pwd(), ".linecablemodels", "fem", "cable_380kv_armoured_ac_flat", Dates.format(now(), "yyyymmdd-HHMMSS")
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
        gridpoint = (id = nothing, case = "cable_380kv_armoured_ac_flat", model = "fem"),
        length_unit = :base, clip = false),
    joinpath(output_directory, "fem.json"))
export_data(:json, problem; file_name = joinpath(output_directory, "problem.json"))

coordinates = details(analytical_result).data.coordinates
coordinates == details(fem_result).data.coordinates ||
    error("Analytical and FEM conductor row/column orders differ")
open(joinpath(output_directory, "matrices.csv"), "w") do io
    println(io, "# Case: cable_380kv_armoured_ac_flat")
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
println(comparison_io, "Case: cable_380kv_armoured_ac_flat")
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
