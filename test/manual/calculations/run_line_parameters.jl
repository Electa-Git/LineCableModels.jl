# Disposable manual runner. Nothing here is part of the package or Gauntlet.
#
# Open this file in the IDE and run/include it. This script does not activate
# or replace the IDE's active project environment.
#
# One-time setup in the user's Julia 1.12 default environment (`@v1.12`):
#     import Pkg
#     Pkg.develop(path = "/home/amartins/Documents/KUL/LineCableModels")
# Then run this file with that consumer environment active. `develop` is used
# only because LineCableModels is not registered yet; dependency resolution is
# otherwise the same as for a registry install.
#
# Edit `run_fem`, `case_id`, and `frequency_grid` below, then re-include this
# file in a fresh session when changing the backend or case.

repository = normpath(joinpath(@__DIR__, "..", "..", ".."))
gauntlet_project = joinpath(repository, "gauntlet")
gauntlet_project in LOAD_PATH || push!(LOAD_PATH, gauntlet_project)

# The one backend switch: false runs the analytical/coaxial solver, true runs FEM.
run_fem = false

# Before: Revise was mandatory even for a one-shot include. It is an optional
# IDE convenience; load it in your session if you want live package revision.
using LineCableModels
run_fem && (@eval using Gmsh)
isdefined(@__MODULE__, :Gauntlet) || include(joinpath(gauntlet_project, "Gauntlet.jl"))

manual_output = joinpath(get(ENV, "LINECABLEMODELS_MANUAL_OUTPUT",
    joinpath(tempdir(), "linecablemodels-manual")), "line-parameters")
mkpath(manual_output)
fullfile(filename) = joinpath(manual_output, filename); #hide

# Edit these inputs for a different catalogue case or frequency sweep.
case_id = :cable_220kv_milliken_1x2500_252_trefoil
frequency_grid = 10.0 .^ range(-1, 7; length = 101)  # Hz; use [50.0] for one frequency.
loaded_case = Gauntlet.load_case(case_id;
    variation = Gauntlet.ExactOverrides(; frequencies = frequency_grid))
problem = loaded_case.problem
system = problem.system
earth = problem.earth_props

println("\n", run_fem ? "FEM" : "Analytical", ": ", case_id,
    " | ", length(frequency_grid), " frequencies")
println("Terminal order: ", loaded_case.port_order)

output_file = fullfile("pscad_export.pscx")
export_file = export_data(:pscad, system, earth, file_name = output_file);

if run_fem
    # FEM field-model choices. Execution controls belong to compute(...; options).
    fem_formulation = Formulation(:LineCableModelsFEM;
        options = (
            physics = :quasi_tem,
            reduce_bundle = false,
            kron_reduction = false,
            ideal_transposition = false
        ))
    fem_options = (
        mesh_policy = :remesh,
        resume_run_directory = nothing,
        keep_run_directory = true,
        trace = true,
        output_basis = :pul,
        verbosity = (default = 1,),
        gmsh_verbosity = 2,
        getdp_verbosity = 4,
        frequency_workers = 2,
        solver_threads = 1,
        plot_field_maps = false
    )

    fem_result = @time compute(problem, fem_formulation; options = fem_options)
    Zfem = observe(fem_result, Z)  # ohm/m
    Yfem = observe(fem_result, Y)  # S/m
    println("Run directory: ", details(fem_result).data.fem.run.run_directory)
    println("Zfem / Yfem: ", size(Zfem), " | result: fem_result")
end

# Analytical solver defaults. These are the same selections as Formulation().
analytical_formulation = Formulation(
    internal_impedance = :default,
    insulation_impedance = :default,
    earth_impedance = :default,
    shunt_model = :default,  # Coaxial annuli; :boundary opts into the local field solve.
    insulation_admittance = :default,
    semicon_admittance = :default,
    earth_admittance = :default,
    earth_properties = :default,
    pipe_impedance = :default,
    temperature_dependence = :default,
    options = (
        reduce_bundle = false,
        kron_reduction = false,
        ideal_transposition = false
    )
)
analytical_result = @time compute(problem, analytical_formulation)
Zan = observe(analytical_result, Z)  # ohm/m
Yan = observe(analytical_result, Y)  # S/m
println("Zresult / Yresult: ", size(Zan), " | result: analytical_result")

# `Results remain available in the REPL.
nothing
