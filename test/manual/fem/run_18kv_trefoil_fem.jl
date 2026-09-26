# Manual run: include("test/manual/fem/run_18kv_trefoil_fem.jl")
# Uses the current catalogue case, not saved campaign results. No gauntlet run.
# Before: this include activated Gauntlet. Now preserve the IDE's active project.
gauntlet_project = normpath(joinpath(@__DIR__, "..", "..", "..", "gauntlet"))
gauntlet_project in LOAD_PATH || push!(LOAD_PATH, gauntlet_project)
using LineCableModels, Gmsh
isdefined(@__MODULE__, :Gauntlet) ||
    include(joinpath(gauntlet_project, "Gauntlet.jl"))

# Edit these inputs and re-include. Each inclusion starts a fresh FEM run.
case_id = :cable_18kv_1000mm2_trefoil
frequency_grid = 10.0 .^ range(-1, 7; length=101)  # Hz; use [1e7] for one frequency.
loaded_case = Gauntlet.load_case(case_id;
    variation=Gauntlet.ExactOverrides(; frequencies=frequency_grid))
problem = loaded_case.problem
system = problem.system

# Keep core, wire screen and aluminium foil as separate terminals on each cable.
# Explicitly selects quasi-tem.pro, NOT the coupled quasi-full.pro.
fem_formulation = Formulation(:LineCableModelsFEM;
    options=(
        physics=:quasi_tem,
        reduce_bundle=false,
        kron_reduction=false,
        ideal_transposition=false,
    ))
fem_options = (
    mesh_policy=:remesh,
    resume_run_directory=nothing,
    keep_run_directory=true,
    trace=true,
    output_basis=:pul,
    verbosity=(default=1,),
    gmsh_verbosity=2,
    getdp_verbosity=4,
    frequency_workers=2,
    solver_threads=1,
    plot_field_maps=false,  # Set true to retain field maps for every excitation.
)

# Run through the owned FEM API. Campaign results are not read or overwritten.
println("\nFEM quasi-TEM: ", case_id, " | ", length(frequency_grid), " frequencies")
println("Terminal order: ", loaded_case.port_order)
fem_result = compute(problem, fem_formulation; options=fem_options)

# Raw per-unit-length matrices, indexed [response, source, frequency], for the IDE.
Zfem = observe(fem_result, Z)  # ohm/m
Yfem = observe(fem_result, Y)  # S/m
Pfem = details(fem_result).data.fem.primitive.P_primitive  # ohm m; Y = inv(P) per frequency.
run_directory = details(fem_result).data.fem.run.run_directory
println("\nRun directory: ", run_directory)
println("Zfem / Yfem / Pfem: ", size(Zfem), " | result: fem_result")
# Meshes, solver inputs, raw columns and logs remain under run_directory.
# Optional matrix plots after the solve:
# using GLMakie
# Before: blocks selected matrix pages. Now layout is the only panel capacity.
# plots = LineCableModels.plot(fem_result; ydata=(R, X, G, B), layout=(3, 3))
nothing
