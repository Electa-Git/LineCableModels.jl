# Coaxial cable with a core of bent rectangular strands, sheath and steel armour
# Run from the repository root:
#     julia --project=test test/manual/cases/coaxial_rectangular_strands.jl
# Results: .linecablemodels/fem/coaxial_rectangular_strands/yyyymmdd-HHMMSS/
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

# Dimensions in m. The core is a straight centre wire and courses of 6, 12 and 18 members,
# given explicitly; the layers around it are the same in both coaxial cases.
const wire_diameter = 3.0e-3
const strand_width = 3.0e-3
const strand_height = 3.0e-3
const course_members = (6, 12, 18)
const core_lay_ratio = 12.0
const semicon_thickness = 1.0e-3
const insulation_thickness = 8.0e-3
const sheath_thickness = 1.0e-3
const bedding_thickness = 3.0e-3
const armour_wire_diameter = 4.0e-3
const armour_wires = 40
const armour_lay_ratio = 10.0
const jacket_thickness = 4.0e-3
const copper_rho = 1.7241e-8
const copper_alpha = 0.00393
const semicon_rho = 1000.0
const semicon_eps_r = 1000.0
const polymer_rho = 1.0e14
const polymer_eps_r = 2.3
const steel_rho = 1.7475e-7
const steel_mu_r = 36.6326
const steel_alpha = 0.0045
const cable_x = 0.0
const cable_y = -1.0
const line_length = 1000.0
const temperature = 20.0
const earth_rho = 100.0
const earth_eps_r = 10.0
const frequencies = 10.0 .^ range(-1, 7; length = 101)
# 11 points, one per 0.8 decade, 0.1 Hz to 10 MHz; all shared with the analytical grid.
const fem_frequency_indices = 1:10:101

copper = LineCableModels.Material(:conductor, copper_rho, 1.0, 1.0, 20.0, copper_alpha)
semicon = LineCableModels.Material(:semicon, semicon_rho, semicon_eps_r)
xlpe = LineCableModels.Material(:insulator, polymer_rho, polymer_eps_r)
pe = LineCableModels.Material(:insulator, polymer_rho, polymer_eps_r)
steel = LineCableModels.Material(:conductor, steel_rho, 1.0, steel_mu_r, 20.0, steel_alpha)
fill = LineCableModels.Material(kind = :insulator, rho = Inf, eps_r = 1.0, mu_r = 1.0)

wire_radius = wire_diameter / 2
# Each rectangle is bent into an annular sector of the same area over its course.
strand = LineCableModels.Rectangle(strand_width, strand_height)
courses = [LineCableModels.Group(
               Symbol(:core_course_, k),
               LineCableModels.Region(Symbol(:core_strands_, k), strand, copper);
               pattern = LineCableModels.Ring(members; r = wire_radius +
                                                           (k - 1 / 2) * strand_height),
               path = LineCableModels.Helix(LineCableModels.LayRatio(core_lay_ratio))
           ) for (k, members) in enumerate(course_members)]
radius = wire_radius + strand_height * length(course_members)
parts = LineCableModels.AbstractCablePart[LineCableModels.Enclosure(
    :core_matrix,
    LineCableModels.terminal(
        :core,
        LineCableModels.Region(:core_centre, LineCableModels.Disk(wire_radius), copper),
        courses...
    );
    primitive = LineCableModels.Disk(radius),
    fill = fill
)]
for (tag, layer_thickness, material) in (
    (:core_semicon_inner, semicon_thickness, semicon),
    (:core_insulation, insulation_thickness, xlpe),
    (:core_semicon_outer, semicon_thickness, semicon)
)
    push!(parts,
        LineCableModels.Region(tag, LineCableModels.Shell(layer_thickness), material))
    global radius += layer_thickness
end
push!(parts,
    LineCableModels.Group(
        :sheath,
        LineCableModels.Region(
            :sheath_copper,
            LineCableModels.Annulus(radius, radius + sheath_thickness),
            copper
        )
    ))
radius += sheath_thickness
push!(parts,
    LineCableModels.Region(:sheath_bedding, LineCableModels.Shell(bedding_thickness), pe))
radius += bedding_thickness
armour = LineCableModels.Group(
    :armour,
    LineCableModels.Region(
        :armour_wires, LineCableModels.Disk(armour_wire_diameter / 2), steel);
    pattern = LineCableModels.Ring(armour_wires; r = radius + armour_wire_diameter / 2),
    path = LineCableModels.Helix(LineCableModels.LayRatio(armour_lay_ratio))
)
push!(parts,
    LineCableModels.Enclosure(
        :armour_matrix,
        armour;
        primitive = LineCableModels.Annulus(radius, radius + armour_wire_diameter),
        fill = fill
    ))
radius += armour_wire_diameter
push!(parts,
    LineCableModels.Region(:armour_jacket, LineCableModels.Shell(jacket_thickness), pe))
design = LineCableModels.build(
    LineCableModels.CableDesign, "coaxial_rectangular_strands", LineCableModels.Stack(parts))
core_members = count(
    region -> region.source.material === copper &&
              !(region.primitive isa LineCableModels.Annulus),
    design.geometry.regions)
core_members == 1 + sum(course_members) ||
    error("The core has $core_members members, not $(1 + sum(course_members))")
cable_radius = LineCableModels.outer_radius(design)

earth = LineCableModels.homogeneous(
    rho = earth_rho, eps_r = earth_eps_r, mu_r = 1.0
)
system = LineCableModels.build(
    LineCableModels.LineCableSystem,
    design,
    LineCableModels.Pose2(cable_x, cable_y);
    connections = Dict(:core => 1, :sheath => 2, :armour => 3),
    system_id = "coaxial_rectangular_strands",
    line_length = line_length
)
problem = LineCableModels.Engine.LineParametersProblem(
    system;
    temperature = temperature,
    earth_props = earth,
    frequencies = frequencies
)
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
    pwd(), ".linecablemodels", "fem", "coaxial_rectangular_strands",
    Dates.format(now(), "yyyymmdd-HHMMSS")
)
mkpath(output_directory)
save(ObservedResult(analytical_result, (Z, Y); length_unit = :base, clip = false),
    joinpath(output_directory, "analytical.json"))
# The FEM formulation type cannot be archived, so this archive keeps the matrices,
# frequencies and conductor labels without the formulation record.
fem_matrices = LineParameters(
    PhaseDomain(), observe(fem_result, Z), observe(fem_result, Y),
    observe(fem_result, LineCableModels.frequencies); basis = :pul,
    details = ComputationDetails(coordinates = details(fem_result).data.coordinates))
save(
    ObservedResult(fem_matrices, (Z, Y);
        gridpoint = (id = nothing, case = "coaxial_rectangular_strands", model = "fem"),
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
    println(io,
        "# Case: Coaxial cable with a core of bent rectangular strands, sheath and steel armour (coaxial_rectangular_strands)")
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
println(comparison_io,
    "Case: Coaxial cable with a core of bent rectangular strands, sheath and steel armour (coaxial_rectangular_strands)")
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
    kappa_r_p = abs(kappa * cable_radius)
    @printf(comparison_io, "%10.4g %9.3f %12s %9.3f %12s %13d %11.3g\n",
        frequencies[analytical_index], e_Z, @sprintf("(%d,%d)", Tuple(Z_worst)...),
        e_Y, @sprintf("(%d,%d)", Tuple(Y_worst)...), G_sign_canary, kappa_r_p)
end
comparison = String(take!(comparison_io))
print(comparison)
write(joinpath(output_directory, "comparison.txt"), comparison)

println("Output directory: ", output_directory)
println("FEM run directory: ", fem_run_directory)
