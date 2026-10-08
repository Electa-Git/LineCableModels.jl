# Disposable public-API timings. Uses the active environment. Never activates,
# installs, updates, exports, or writes Gauntlet campaign artifacts.
# Study dependencies come from the local Gauntlet environment through LOAD_PATH.
gauntlet_project = normpath(joinpath(@__DIR__, "..", "..", "..", "gauntlet"))
gauntlet_project in LOAD_PATH || push!(LOAD_PATH, gauntlet_project)
using LineCableModels, LinearAlgebra
isdefined(@__MODULE__, :Gauntlet) ||
    include(joinpath(gauntlet_project, "Gauntlet.jl"))

shunt_benchmark_cases = (:cable_18kv_1000mm2_trefoil, :cable_132kv_630mm2_flathor)
shunt_benchmark_frequencies = 10.0 .^ range(-1, 7; length = 101)
shunt_benchmark_repeats = 3
shunt_benchmark_boundary = true
shunt_benchmark_rows = NamedTuple[]
println("Julia ", VERSION, "; Julia threads ", Threads.nthreads(),
    "; BLAS threads ", BLAS.get_num_threads())
for case_id in shunt_benchmark_cases
    # Before: nominal_problem silently discarded the edited frequency override.
    # problem is the materialized case with ExactOverrides applied.
    problem = Gauntlet.load_case(case_id;
        variation = Gauntlet.ExactOverrides(
            frequencies = shunt_benchmark_frequencies)).problem
    physical=(reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    equivalent=Formulation(; options = physical)
    println(
        "\n", case_id, " | default equivalent annular layer | ", length(problem.frequencies), " frequencies")
    @time compute(problem, equivalent) # Warm-up. Do not mix with the warm timings.
    for repetition in 1:shunt_benchmark_repeats
        sample=@timed @time compute(problem, equivalent)
        push!(shunt_benchmark_rows,
            (; case_id, model = :equivalent, repetition,
                seconds = sample.time, MiB = sample.bytes/1024^2))
        @assert details(sample.value).data.shunt_model.solves==0
    end
    if shunt_benchmark_boundary
        boundary=Formulation(shunt_model = :boundary; options = physical)
        println("Boundary compute, including blueprint construction; first call")
        @time compute(problem, boundary)
        for repetition in 1:shunt_benchmark_repeats
            sample=@timed @time compute(problem, boundary)
            push!(shunt_benchmark_rows,
                (; case_id, model = :boundary, repetition,
                    seconds = sample.time, MiB = sample.bytes/1024^2))
        end
    end
end
display(shunt_benchmark_rows)
nothing
