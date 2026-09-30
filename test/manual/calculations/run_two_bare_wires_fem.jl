# Manual FEM sweep for two bare wires; installation_z selects the placement.
#
# Run from the repository root with:
#     julia -i --project=test test/manual/calculations/run_two_bare_wires_fem.jl
#
# The script leaves the selected FEM and analytical results in the corresponding
# `case_*_results` arrays, and their admittance tensors in the corresponding
# `case_*_admittance` arrays. FEM run directories are kept by the backend.
# Interactive G/B plots are retained in plot_1 and plot_2. Set save_tables=true
# below to also export CSV/XLSX tables.

# GetDP uses x horizontal, y vertical and z axial. Quasi-fw uses axial current [A] and normalized
# scalar plus vector-path voltage. Raw P is [ohm m], and Y=inv(P) is [S/m].
# Each terminal uses one path to its lowest CAD vertex. Overhead references
# lie on the local air/earth interface; buried references lie at the outer bottom boundary.
# Native extraction checks: test/extensions/fem_native_measurements.jl.
# Saved historical reference comparisons: fem_voltage_reference_notes.md.

using LineCableModels
using Printf
using Gmsh
native_review = any(in(ARGS), ("--compare-native", "--review-native"))
if !native_review
    using Revise
    using GLMakie
end

# Including a script from another checkout does not switch the loaded package.
# Fail before an expensive sweep if this session still owns the old backend.
script_checkout = normpath(joinpath(@__DIR__, "../../.."))
realpath(pkgdir(LineCableModels)) == realpath(script_checkout) || error(
    "This script belongs to $script_checkout, but LineCableModels was loaded from " *
    "$(pkgdir(LineCableModels)). Restart Julia with --project=$script_checkout/test.")
println("LineCableModels source: ", pkgdir(LineCableModels))

# All dimensions below are SI units. Frequencies span 0.1 Hz to 1 MHz.
frequency_grid = collect(10.0 .^ range(-1.0, 6.0, length = 10))
soil_relative_permittivity = 1.0
soil_relative_permeability = 1.0
line_length = 1.0
temperature = 20.0

# Public observation cutoffs in S/m; zero retains every nonzero G/B value.
plot_atol = (G = 0.0, B = 0.0)
save_tables = false
if save_tables
    using XLSX
end

# Case 1: fixed 4.25 cm radius, sweep soil resistivity.
fixed_external_radius = 4.25e-2
soil_resistivity_grid = Grid([0.1, 1.0, 100.0, 1000.0])
# soil_resistivity_grid = Grid([0.1])

# Case 2: fixed 0.1 Ohm.m soil, sweep external radius.  The 4.25 cm point is
# intentionally not rerun here: it is case_1_results[1], at rho = 0.1 Ohm.m.
fixed_soil_resistivity = 0.1
external_radius_grid = Grid([0.1e-2, 1.0e-2, 8.5e-2])
# external_radius_grid = Grid([8.5e-2])

installation_z = 1.0

function build_two_bare_wires_problem(
        external_radius, earth_rho, frequencies; vert = installation_z
)
    # Shared physical setup for the manual voltage-reference study.
    materials = MaterialsLibrary(add_defaults = true)
    copper = Material(materials, :copper)
    design = build(
        CableDesign,
        "two_bare_wires",
        Stack(
            Group(
            :core,
            Region(:core_metal, Disk(external_radius), copper)
        )
        )
    )
    system = build(
        LineCableSystem,
        [design, design],
        [Pose2(0.0, vert), Pose2(1.0, vert)];
        connections = [Dict(:core => 1), Dict(:core => 2)],
        system_id = "two_bare_wires",
        line_length = line_length
    )
    earth = homogeneous(
        rho = earth_rho,
        eps_r = soil_relative_permittivity,
        mu_r = soil_relative_permeability
    )
    return LineParametersProblem(
        system;
        temperature = temperature,
        earth_props = earth,
        frequencies = frequencies
    )
end

function make_problem_space(radius_values, earth_rho_values, frequencies; vert)
    build_point = (radius, earth_rho) -> build_two_bare_wires_problem(
        radius, earth_rho, frequencies; vert = vert)
    return Gridspace{LineParametersProblem}(
        build_point,
        (Grid(radius_values), Grid(earth_rho_values))
    )
end

fem_formulation = Formulation(
    :LineCableModelsFEM;
    options = (
        physics = :quasi_fw,
        reduce_bundle = false,
        kron_reduction = false,
        ideal_transposition = false
    )
)
fem_options = (
    mesh_policy = :remesh,
    resume_run_directory = nothing,
    keep_run_directory = true,
    trace = true,
    output_basis = :pul,
    verbosity = (default = 1,),
    gmsh_verbosity = 0,
    getdp_verbosity = 0,
    # Move the aerial PML farther away to resolve the small conductance signs.
    # Coarser physical elements keep this practical; percentage accuracy is
    # a separate convergence check, including the conductor-contour mesh.
    domain_skin_depths = 24.0,
    # Qualified signs and cost: see fem_pml_physical_mesh.md.
    # Counts follow propagation, clearance and stretch-coefficient variation.
    pml_resolution = (interpolation_cells = 72, coefficient_change = 0.12),
    mesh_size_factor = 5.0,
    # Coarsen the remote buffer and PML tangential mesh; retain local targets.
    exterior_mesh_size_factor = 10.0,
    # Increase to retain wider local interface refinement (minimum/default: 1).
    interface_refinement_factor = 1.0,
    # Fastest measured discretization in fem_native_performance/feature_review.jl.
    # Local runner selection; package defaults remain unchanged.
    pml_element_family = :quadrangle, # or :triangle
    volume_quadrature = 12,         # triangle rule, including triangular PML
    physical_volume_quadrature = 3,  # or nothing (inherit), 4, 7, 12, 13
    pml_quadrature = 9,             # quadrangle Gauss-Legendre: 4, 9, 16
    # Prescribed conductor resolution; independent of exterior coarsening.
    conductor_geometry_tolerance = 1e-3, # 96 segments per full circle
    conductor_skin_depth_elements = 6.0,
    conductor_mesh_growth = sqrt(1.25),
    conductor_skin_depths = 5.0,
    conductor_thickness_elements = 4,
    # Qualified on the unchanged mesh: distribute assembly across frequencies.
    frequency_workers = 4,
    solver_threads = 1,
    plot_field_maps = false
)

if native_review
    # One representative case, same ten-frequency scan. --review-native only
    # reads completed comparisons. Default plots are interactive GLMakie;
    # --save explicitly requests headless SVG/PNG output.
    include("fem_native_performance/feature_review.jl")
    "--compare-native" in ARGS &&
        run_native_review(frequency_grid, fem_formulation, fem_options)
    include("fem_native_performance/feature_review_plots.jl")
    native_review_plots = plot_native_review(frequency_grid)
else
    # One detached case: rho = 0.1 Ohm.m, radius = 4.25 cm, at frequency_grid.
    onelab_entry = export_data(
        :onelab,
        build_two_bare_wires_problem(fixed_external_radius, fixed_soil_resistivity, frequency_grid),
        fem_formulation;
        file_name = joinpath(
            script_checkout, ".linecablemodels", "fem", "onelab-two-bare-wires", "study.pro"),
        mesh_options = fem_options[(
            :domain_skin_depths, :pml_resolution, :mesh_size_factor, :exterior_mesh_size_factor,
            :interface_refinement_factor, :pml_element_family, :volume_quadrature,
            :physical_volume_quadrature, :pml_quadrature,
            :conductor_geometry_tolerance, :conductor_skin_depth_elements,
            :conductor_mesh_growth, :conductor_skin_depths, :conductor_thickness_elements)],
        overwrite = true
    )

    function run_fem_space(space, label, formulation, options)
        println("\nFEM ", label, " | ", length(space), " parameter points | ",
            length(frequency_grid), " frequencies")
        parametric_problem = ParametricProblem(space, ComputationOptions(options))
        return compute(parametric_problem, Combinatorial(formulation))
    end

    function write_admittance_csv(path, parameter_name, parameter_values, results)
        length(parameter_values) == length(results) ||
            throw(DimensionMismatch("parameter values and FEM results must have equal length"))
        mkpath(dirname(path))
        open(path, "w") do io
            println(io, parameter_name,
                ",frequency_hz,response_terminal,basis_terminal,real_s_per_m,imaginary_s_per_m")
            for (parameter, result) in zip(parameter_values, results)
                frequencies = observe(result, LineCableModels.frequencies)
                admittance = observe(result, Y)
                for frequency_index in eachindex(frequencies),
                    response in axes(admittance, 1), basis in axes(admittance, 2)
                    value = admittance[response, basis, frequency_index]
                    @printf(io,
                        "%.17g,%.17g,%d,%d,%.17g,%.17g\n",
                        parameter,
                        frequencies[frequency_index],
                        response,
                        basis,
                        real(value),
                        imag(value))
                end
            end
        end
        return path
    end

    function write_impedance_csv(path, parameter_name, parameter_values, results)
        length(parameter_values) == length(results) ||
            throw(DimensionMismatch("parameter values and FEM results must have equal length"))
        mkpath(dirname(path))
        open(path, "w") do io
            println(io, parameter_name,
                ",frequency_hz,response_terminal,basis_terminal,real_ohm_per_m,imaginary_ohm_per_m")
            for (parameter, result) in zip(parameter_values, results)
                frequencies = observe(result, LineCableModels.frequencies)
                impedance = observe(result, Z)
                for frequency_index in eachindex(frequencies),
                    response in axes(impedance, 1), basis in axes(impedance, 2)
                    value = impedance[response, basis, frequency_index]
                    @printf(io,
                        "%.17g,%.17g,%d,%d,%.17g,%.17g\n",
                        parameter,
                        frequencies[frequency_index],
                        response,
                        basis,
                        real(value),
                        imag(value))
                end
            end
        end
        return path
    end

    function write_line_parameters_xlsx(directory, file_prefix, parameter_values, results)
        length(parameter_values) == length(results) ||
            throw(DimensionMismatch("parameter values and XLSX results must have equal length"))
        mkpath(directory)
        paths = String[]
        for (parameter, result) in zip(parameter_values, results)
            token = replace(@sprintf("%.6g", float(parameter)), "." => "p", "-" => "m")
            path = joinpath(directory, "$(file_prefix)_$(token).xlsx")
            observed = ObservedResult(result, (R, X, G, B); length_unit = :base, clip = false)
            append!(paths, export_data(:xlsx, observed; file_name = path))
        end
        return paths
    end

    case_1_space = make_problem_space(
        Grid([fixed_external_radius]), soil_resistivity_grid, frequency_grid;
        vert = installation_z
    )
    case_1_result = run_fem_space(
        case_1_space, "case 1 (soil resistivity)", fem_formulation, fem_options
    )
    case_1_results = collect(case_1_result)
    case_1_admittance = [observe(result, Y) for result in case_1_results]
    case_1_run_directories = [details(result).data.fem.run.run_directory
                              for result in case_1_results]

    case_2_space = make_problem_space(
        external_radius_grid, Grid([fixed_soil_resistivity]), frequency_grid;
        vert = installation_z
    )
    case_2_result = run_fem_space(
        case_2_space, "case 2 (external radius)", fem_formulation, fem_options
    )
    case_2_results = collect(case_2_result)
    case_2_admittance = [observe(result, Y) for result in case_2_results]
    case_2_run_directories = [details(result).data.fem.run.run_directory
                              for result in case_2_results]

    # Default analytical/coaxial formulation, evaluated over the same two problem
    # spaces so that each result has the same parameter ordering as its FEM peer.
    analytical_formulation = Formulation(options = (
        reduce_bundle = false,
        kron_reduction = false,
        ideal_transposition = false
    ))

    function run_analytical_space(space, label, formulation)
        println("\nAnalytical ", label, " | ", length(space), " parameter points | ",
            length(frequency_grid), " frequencies")
        return compute(ParametricProblem(space), Combinatorial(formulation))
    end

    case_1_analytical_result = run_analytical_space(
        case_1_space, "case 1 (soil resistivity)", analytical_formulation
    )
    case_1_analytical_results = collect(case_1_analytical_result)
    case_1_analytical_admittance = [observe(result, Y)
                                    for result in case_1_analytical_results]

    case_2_analytical_result = run_analytical_space(
        case_2_space, "case 2 (external radius)", analytical_formulation
    )
    case_2_analytical_results = collect(case_2_analytical_result)
    case_2_analytical_admittance = [observe(result, Y)
                                    for result in case_2_analytical_results]

    # Each Y request opens separate G and B figures with the standard plot controls.
    # Keep the trailing semicolons on plot assignments: VS Code otherwise prints
    # the returned handles and their Makie internals when evaluating a statement.
    case_1_plot_results = vcat(case_1_results, case_1_analytical_results)
    case_1_plot_labels = vcat(
        ["FEM: ρ = $(ρ) Ω·m" for ρ in soil_resistivity_grid],
        ["Analytical: ρ = $(ρ) Ω·m" for ρ in soil_resistivity_grid]
    )
    case_1_plot_styles = vcat(
        fill((linestyle = :solid,), length(case_1_results)),
        fill((linestyle = :dash,), length(case_1_analytical_results))
    )
    plot_1 = Makie.plot(
        case_1_plot_results...;
        ydata = (Y,),
        atol = plot_atol,
        length_unit = :base,
        series_labels = case_1_plot_labels,
        series_attributes = case_1_plot_styles,
        layout = (2, 2),
        xscale = :log10,
        yscale = :log10,
        backend = :gl,
        display_plot = true,
        figure_title = "Two bare wires — case 1: soil resistivity",
        legend_title = "Model and physical assumptions",
        legend_attributes = (; valign = :center),
        legend_overflow = :show_all
    )

    case_2_plot_results = vcat(case_2_results, case_2_analytical_results)
    case_2_plot_labels = vcat(
        ["FEM: r = $(100r) cm" for r in external_radius_grid],
        ["Analytical: r = $(100r) cm" for r in external_radius_grid]
    )
    case_2_plot_styles = vcat(
        fill((linestyle = :solid,), length(case_2_results)),
        fill((linestyle = :dash,), length(case_2_analytical_results))
    )
    plot_2 = Makie.plot(
        case_2_plot_results...;
        ydata = (Y,),
        atol = plot_atol,
        length_unit = :base,
        series_labels = case_2_plot_labels,
        series_attributes = case_2_plot_styles,
        layout = (2, 2),
        xscale = :log10,
        yscale = :log10,
        backend = :gl,
        display_plot = true,
        figure_title = "Two bare wires — case 2: external radius",
        legend_title = "Model and physical assumptions",
        legend_attributes = (; valign = :center),
        legend_overflow = :show_all
    )

    if save_tables
        output_directory = joinpath(
            get(ENV, "LINECABLEMODELS_MANUAL_OUTPUT",
                joinpath(tempdir(), "linecablemodels-manual")),
            "two-bare-wires-fem")
        case_1_csv = write_admittance_csv(
            joinpath(output_directory, "case_1_soil_resistivity.csv"),
            "soil_resistivity_ohm_m",
            collect(soil_resistivity_grid),
            case_1_results
        )
        case_2_csv = write_admittance_csv(
            joinpath(output_directory, "case_2_external_radius.csv"),
            "external_radius_m",
            collect(external_radius_grid),
            case_2_results
        )
        case_1_impedance_csv = write_impedance_csv(
            joinpath(output_directory, "case_1_soil_resistivity_impedance.csv"),
            "soil_resistivity_ohm_m",
            collect(soil_resistivity_grid),
            case_1_results
        )
        case_2_impedance_csv = write_impedance_csv(
            joinpath(output_directory, "case_2_external_radius_impedance.csv"),
            "external_radius_m",
            collect(external_radius_grid),
            case_2_results
        )

        xlsx_directory = joinpath(output_directory, "xlsx")
        case_1_fem_xlsx = write_line_parameters_xlsx(
            xlsx_directory,
            "case_1_fem_rho_ohm_m",
            collect(soil_resistivity_grid),
            case_1_results
        )
        case_1_analytical_xlsx = write_line_parameters_xlsx(
            xlsx_directory,
            "case_1_analytical_rho_ohm_m",
            collect(soil_resistivity_grid),
            case_1_analytical_results
        )
        case_2_fem_xlsx = write_line_parameters_xlsx(
            xlsx_directory,
            "case_2_fem_radius_m",
            collect(external_radius_grid),
            case_2_results
        )
        case_2_analytical_xlsx = write_line_parameters_xlsx(
            xlsx_directory,
            "case_2_analytical_radius_m",
            collect(external_radius_grid),
            case_2_analytical_results
        )

        println("Case 1 CSV: ", case_1_csv)
        println("Case 2 CSV: ", case_2_csv)
        println("Case 1 impedance CSV: ", case_1_impedance_csv)
        println("Case 2 impedance CSV: ", case_2_impedance_csv)
        println("Case 1 FEM XLSX: ", case_1_fem_xlsx)
        println("Case 1 analytical XLSX: ", case_1_analytical_xlsx)
        println("Case 2 FEM XLSX: ", case_2_fem_xlsx)
        println("Case 2 analytical XLSX: ", case_2_analytical_xlsx)
    end

    println("\nCompleted.")
    println("Case 1 admittance tensors: ", [size(value) for value in case_1_admittance])
    println("Case 2 admittance tensors: ", [size(value) for value in case_2_admittance])
    println("The 4.25 cm, 0.1 Ohm.m point for case 2 is case_1_results[1].")
end # ordinary parametric sweep

using LineCableModels, Gmsh, GLMakie

mesh = import_data(:msh, ".linecablemodels/fem/onelab-two-bare-wires/study.msh")
p = LineCableModels.plot(mesh; color_by = :physical, inspect = :element, backend = :gl)

nothing
