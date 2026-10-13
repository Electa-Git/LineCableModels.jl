# cable 18kv 1000mm2 trefoil homogenized
# Run from the repository root:
#     julia --project=test test/manual/cases/cable_18kv_1000mm2_trefoil_homogenized.jl
# Results: .linecablemodels/fem/cable_18kv_1000mm2_trefoil_homogenized/yyyymmdd-HHMMSS/
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

const core_strand_diameter = 0.00419
const core_ring_lay_ratio_1 = 15.0
const core_ring_lay_ratio_2 = 13.5
const core_ring_lay_ratio_3 = 12.5
const core_ring_lay_ratio_4 = 11.0
const core_ring_wire_counts = (1, 6, 12, 18, 24)
const screen_wires = 49
const screen_wire_diameter = 0.00094
const screen_wire_lay_ratio = 10.0
const semicon_tape_thickness = 0.0003
const inner_semicon_thickness = 0.000768
const insulation_thickness = 0.0083
const outer_semicon_thickness = 0.000472
const copper_tape_thickness = 0.0001
const copper_tape_width = 0.01
const copper_tape_lay_ratio = 10.0
const water_blocking_thickness = 0.00094
const aluminum_tape_thickness = 0.00015
const pe_face_thickness = 5.0e-5
const jacket_thickness = 0.0034
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
    problem = let
    materials = LineCableModels.MaterialsLibrary(add_defaults = true)
    aluminum = LineCableModels.Material(materials, :aluminum)
    copper = LineCableModels.Material(materials, :copper)
    pe = LineCableModels.Material(materials, :pe)
    xlpe = LineCableModels.Material(materials, :xlpe)
    semicon1 = LineCableModels.Material(materials, :semicon1)
    semicon2 = LineCableModels.Material(materials, :semicon2)
    polyacrylate = LineCableModels.Material(materials, :polyacrylate)
    matrix = LineCableModels.Material(
        kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0
    )

    counts = core_ring_wire_counts
    lay_ratios = LineCableModels.LayRatio(
        core_ring_lay_ratio_1,
        core_ring_lay_ratio_2,
        core_ring_lay_ratio_3,
        core_ring_lay_ratio_4
    )
    wire_radius = core_strand_diameter / 2
    radius = (2length(counts) - 1) * wire_radius
    stranded_core = LineCableModels.stranded(
        aluminum;
        shape = LineCableModels.Disk(wire_radius),
        lay = lay_ratios,
        boundary = LineCableModels.Disk(radius),
        fill = matrix
    )
    core = stranded_core
    parts = LineCableModels.AbstractCablePart[core]
    for (tag, thickness, material) in (
        (:core_semicon_tape_inner, semicon_tape_thickness, polyacrylate),
        (:core_semicon_inner, inner_semicon_thickness, semicon1),
        (:core_insulation, insulation_thickness, xlpe),
        (:core_semicon_outer, outer_semicon_thickness, semicon2),
        (:core_semicon_tape_outer, semicon_tape_thickness, polyacrylate)
    )
        push!(parts, LineCableModels.Region(
            tag, LineCableModels.Shell(thickness), material
        ))
        radius += thickness
    end

    screen_radius = screen_wire_diameter / 2
    screen_inner = radius
    screen_centre = radius + screen_radius
    screen = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_wires, LineCableModels.Disk(screen_radius), copper
        );
        pattern = LineCableModels.Ring(screen_wires; r = screen_centre),
        path = LineCableModels.Helix(LineCableModels.LayRatio(screen_wire_lay_ratio))
    )
    radius += 2screen_radius
    tape_outer = radius + copper_tape_thickness
    tape = LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_copper_tape,
            LineCableModels.Rectangle(
                copper_tape_width,
                copper_tape_thickness
            ),
            copper
        );
        pattern = LineCableModels.Ring(1; r = (radius + tape_outer) / 2),
        path = LineCableModels.Helix(LineCableModels.LayRatio(copper_tape_lay_ratio))
    )
    push!(parts,
        LineCableModels.Enclosure(
            :screen_matrix,
            LineCableModels.Stack(screen, tape);
            primitive = LineCableModels.Annulus(screen_inner, tape_outer),
            fill = matrix
        ))
    radius = tape_outer
    push!(parts,
        LineCableModels.Region(
            :sheath_water_blocking,
            LineCableModels.Shell(water_blocking_thickness),
            polyacrylate
        ))
    radius += water_blocking_thickness

    aluminum_outer = radius + aluminum_tape_thickness
    push!(parts,
        LineCableModels.Group(
            :jacket,
            LineCableModels.Region(
                :jacket_aluminum_tape,
                LineCableModels.Annulus(radius, aluminum_outer),
                aluminum
            )
        ))
    push!(parts,
        LineCableModels.Region(
            :jacket_pe_face, LineCableModels.Shell(pe_face_thickness), pe
        ))
    push!(parts,
        LineCableModels.Region(
            :jacket_insulation, LineCableModels.Shell(jacket_thickness), pe
        ))
    design = LineCableModels.build(
        LineCableModels.CableDesign,
        "18kV_1000mm2",
        LineCableModels.Stack(parts)
    )
    formation_spacing = (2 + formation_clearance_ratio) *
                        LineCableModels.outer_radius(design)
    earth = LineCableModels.homogeneous(
        rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
    )
    placements = LineCableModels.trefoil(
        design;
        center = LineCableModels.at(formation_x, formation_y),
        spacing = formation_spacing,
        connections = (
            core = (1, 4, 7),
            sheath = (2, 5, 8),
            jacket = (3, 6, 9)
        )
    )
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        placements;
        system_id = "cable_18kv_1000mm2_trefoil",
        line_length = line_length
    )
    LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = temperature,
        earth_props = earth,
        frequencies = frequencies
    )
    end
    source = problem.system
    designs = map(LineCableModels.homogenize, source.designs)
    system = LineCableModels.build(
        LineCableModels.LineCableSystem,
        designs,
        source.positions;
        connections = source.connections,
        environment = source.environment,
        system_id = "cable_18kv_1000mm2_trefoil_homogenized",
        line_length = source.line_length
    )
    LineCableModels.Engine.LineParametersProblem(
        system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies
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
    pwd(), ".linecablemodels", "fem", "cable_18kv_1000mm2_trefoil_homogenized", Dates.format(now(), "yyyymmdd-HHMMSS")
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
        gridpoint = (id = nothing, case = "cable_18kv_1000mm2_trefoil_homogenized", model = "fem"),
        length_unit = :base, clip = false),
    joinpath(output_directory, "fem.json"))
export_data(:json, problem; file_name = joinpath(output_directory, "problem.json"))

coordinates = details(analytical_result).data.coordinates
coordinates == details(fem_result).data.coordinates ||
    error("Analytical and FEM conductor row/column orders differ")
coordinates == details(boundary_result).data.coordinates ||
    error("Boundary and default analytical conductor orders differ")
open(joinpath(output_directory, "matrices.csv"), "w") do io
    println(io, "# Case: cable_18kv_1000mm2_trefoil_homogenized")
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
println(comparison_io, "Case: cable_18kv_1000mm2_trefoil_homogenized")
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
