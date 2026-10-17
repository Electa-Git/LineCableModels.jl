# Disposable local-shunt inspection. No FEM or PSCAD launches.
# Run from the repository root:
# `julia --project=test test/manual/performance/benchmark_local_shunt.jl`
# It times the boundary shunt stage of one known case, cold and warm, and shows the numerical
# record of each solve and the capacitance couplings of the first block. The known-case items
# compare the results of the `:boundary` shunt model with FEM. This script only times and
# inspects.
using LineCableModels, LinearAlgebra, DataFrames
include(joinpath(@__DIR__, "..", "..", "support", "scenarios.jl"))
using .CurrentScenarios: known_case_problem

shunt_case = :cable_18kv_1000mm2_trefoil
shunt_warm_repeats = 1 # Increase for stable coefficient-construction timings. Each is a fresh solve.
shunt_blas_threads = 2
BLAS.set_num_threads(shunt_blas_threads)
println("Julia threads: ", Threads.nthreads(), "; BLAS threads: ", BLAS.get_num_threads())

shunt_problem = known_case_problem(shunt_case; frequencies = [50.0])
shunt_designs = shunt_problem.system.designs
shunt_blueprints = [LineCableModels.DataModel.flatten(design, Float64) for design in shunt_designs]
shunt_formulation = Formulation(shunt_model = :boundary)
shunt_response() = only(LineCableModels.Engine.internal_potential_coefficients(
    shunt_designs, shunt_blueprints, [shunt_formulation]))

shunt_cold = @timed shunt_response()
shunt_completed = shunt_cold.value
shunt_measurements = [(pass = "cold", seconds = shunt_cold.time,
    allocated_MiB = shunt_cold.bytes/1024^2, solves = shunt_completed.details.solves)]
for repetition in 1:shunt_warm_repeats
    sample = @timed shunt_response()
    push!(shunt_measurements, (pass = "warm $repetition", seconds = sample.time,
        allocated_MiB = sample.bytes/1024^2, solves = sample.value.details.solves))
end
shunt_C = first(shunt_completed.blocks).C
display(DataFrame(shunt_measurements))
display(DataFrame(shunt_completed.details.diagnostics))
display(DataFrame(coupling = ["inner-open", "inner-reference", "open-reference"],
    nF_per_m = [-shunt_C[1, 2], sum(shunt_C[1, :]), sum(shunt_C[2, :])] .* 1e9))
println("Allocated MiB is cumulative allocation, not peak or live memory.")

# The `:boundary` shunt model changes the shunt admittance only.
shunt_physical = (reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
shunt_sweep_problem = known_case_problem(shunt_case;
    frequencies = 10.0 .^ range(-1, 7; length = 11))
shunt_boundary = compute(shunt_sweep_problem,
    Formulation(shunt_model = :boundary; options = shunt_physical))
shunt_equivalent = compute(shunt_sweep_problem,
    Formulation(shunt_model = :equivalent; options = shunt_physical))
@assert observe(shunt_boundary, Z) == observe(shunt_equivalent, Z)
println("Boundary solves in the sweep: ", details(shunt_boundary).data.shunt_model.solves)
