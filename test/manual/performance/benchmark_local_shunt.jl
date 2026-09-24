# Disposable local-shunt audit. No FEM/PSCAD launches and no campaign writes.
# Run from a fresh Julia:
# julia --project=gauntlet test/manual/performance/benchmark_local_shunt.jl
gauntlet_project=normpath(joinpath(@__DIR__,"..","..","..","gauntlet"))
gauntlet_project in LOAD_PATH || push!(LOAD_PATH,gauntlet_project)
using LineCableModels, LinearAlgebra, BenchmarkTools, DataFrames, TOML, JLD2, SHA
isdefined(@__MODULE__,:Gauntlet) || include(joinpath(gauntlet_project,"Gauntlet.jl"))

shunt_case = :cable_18kv_1000mm2_trefoil
shunt_warm_repeats = 1 # Increase for stable preparation timings; each is a fresh solve.
shunt_compute_sweep = true
shunt_compare_saved = true
shunt_saved_folder = joinpath(gauntlet_project,".work","all-references",
    "benchmark_18kv_1000mm2_trefoil_fem")
shunt_blas_threads = 2
BLAS.set_num_threads(shunt_blas_threads)
println("Julia threads: ",Threads.nthreads(),"; BLAS threads: ",BLAS.get_num_threads())

# Numerical inspection only; this does not validate the auxiliary solver files.
function shunt_read_arrays(path)
    bytes2hex(open(sha256,path)) == first(split(read(path*".sha256",String))) ||
        error("Changed numerical payload: $path")
    jldopen(path,"r") do file
        file["basis"] === :pul || error("Per-unit-length results required")
        LineParameters(PhaseDomain,file["Z"],file["Y"],file["frequencies"];basis=:pul)
    end
end

shunt_loaded = Gauntlet.load_case(shunt_case)
shunt_problem = shunt_loaded.nominal_problem
shunt_engine = LineCableModels.Engine
shunt_blueprints = shunt_engine.flatten.(Ref(LineCableModelsCoaxial()),shunt_problem.system.designs)
shunt_domains = shunt_engine.ShuntModel.internal_shunt_domains(shunt_problem.system.designs,shunt_blueprints)
isempty(shunt_domains) && error("No qualified local domains in $shunt_case")
shunt_selection = formula(:boundary;options=(audit=true,))
shunt_prepare() = only(shunt_engine.flatten(LineCableModelsCoaxial(),
    shunt_problem.system.designs, eltype(shunt_problem), [Formulation(shunt_model=shunt_selection)]))

shunt_cold = @timed shunt_prepare()
shunt_prepared = shunt_engine.LocalCableData(shunt_cold.value)
shunt_measurements = [(pass="cold",seconds=shunt_cold.time,allocated_MiB=shunt_cold.bytes/1024^2,
    solves=shunt_prepared.shunt_details.solves)]
for repetition in 1:shunt_warm_repeats
    sample = @timed shunt_prepare()
    push!(shunt_measurements,(pass="warm $repetition",seconds=sample.time,
        allocated_MiB=sample.bytes/1024^2,solves=shunt_engine.LocalCableData(sample.value).shunt_details.solves))
end
shunt_timing_df = DataFrame(shunt_measurements)
shunt_diagnostics_df = DataFrame(shunt_prepared.shunt_details.diagnostics)
shunt_C = first(shunt_prepared.shunt).C
shunt_couplings_df = DataFrame(coupling=["inner-open","inner-reference","open-reference"],
    nF_per_m=[-shunt_C[1,2],sum(shunt_C[1,:]),sum(shunt_C[2,:])].*1e9)
display(shunt_timing_df)
display(shunt_diagnostics_df)
display(shunt_couplings_df)
println("Allocated MiB is cumulative allocation, not peak/live memory. Dense matrix bytes are in diagnostics.")

if shunt_compute_sweep
    shunt_physical = (reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    # Before: prescribed Γ was under parameters. Now Unified owns it in options;
    # its value remains verbatim and this is not a comparison with withdrawn formulas.
    shunt_formulations = [Formulation(shunt_model=shunt_selection;options=shunt_physical),
        Formulation(shunt_model=shunt_selection,
            earth_impedance=formula(:unified;options=(Γ=1e-4im,));options=shunt_physical)]
    shunt_sweep = @timed compute(shunt_problem,shunt_formulations;options=(trace=true,))
    shunt_results = shunt_sweep.value
    shunt_default = first(shunt_results)
    @assert details(shunt_default).data.shunt_model.solves == shunt_prepared.shunt_details.solves
    shunt_coaxial = Formulation(shunt_model=:coaxial;options=shunt_physical)
    shunt_annular = compute(shunt_problem,shunt_coaxial;options=(trace=true,))
    @assert observe(shunt_default,Z) == observe(shunt_annular,Z)
    @assert details(shunt_default).data.trace.Pg == details(shunt_annular).data.trace.Pg
    shunt_loop_input = shunt_engine.lineinput(shunt_problem, shunt_cold.value)
    shunt_loop_options = computation_options(LineCableModelsCoaxial, ComputationOptions(trace=true))
    # Before: the private loop accepted five arguments. Now completion also takes
    # the captured physical inputs and original point identity. Capture once here,
    # just as the engine does, to keep this measurement's preparation scope honest.
    shunt_inputs = shunt_engine.completed_inputs(shunt_problem)
    shunt_id = LineCableModels.Grammar.gridpoint_id()
    shunt_loop_trial = @benchmark shunt_engine._compute(LineCableModelsCoaxial(),
        $shunt_problem,$(first(shunt_formulations)),$shunt_loop_options,$shunt_loop_input,
        $shunt_inputs,$shunt_id) samples=3 evals=1
    shunt_annular_trial = @benchmark compute($shunt_problem,$shunt_coaxial;
        options=(trace=true,)) samples=3 evals=1
    display(shunt_loop_trial)
    display(shunt_annular_trial)
    println("End-to-end formulation collection: ",shunt_sweep.time," s; ",shunt_sweep.bytes/1024^2," MiB allocated.")
    if shunt_compare_saved && isfile(joinpath(shunt_saved_folder,"state.toml"))
        shunt_saved_state = TOML.parsefile(joinpath(shunt_saved_folder,"state.toml"))
        shunt_attempt = joinpath(shunt_saved_folder,shunt_saved_state["current"])
        shunt_reference = shunt_read_arrays(joinpath(shunt_attempt,"reference","calculation.jld2"))
        shunt_saved = shunt_read_arrays(joinpath(shunt_attempt,"candidate","points","1","calculation.jld2"))
        shunt_comparison_rows = NamedTuple[]
        for (model,value) in (("saved annular",shunt_saved),("integrated finite strip",shunt_default))
            comparisons = map(band->shunt_engine.compare(shunt_reference,value,B;band),
                (:all,:dc,:harmonic,:narrow,:wide))
            for (coupling,i,j) in (("core-screen",1,2),("core-foil",1,3),("screen-foil",2,3))
                errors = NamedTuple{(:all,:dc,:harmonic,:narrow,:wide)}(
                    Tuple(100comparison.relative[i,j] for comparison in comparisons))
                push!(shunt_comparison_rows,merge((;model,coupling),errors))
            end
        end
        shunt_comparison_df = DataFrame(shunt_comparison_rows)
        println("B relative RMS [%] versus saved FEM; this is a model comparison, not a convergence bound.")
        show(stdout,MIME"text/plain"(),shunt_comparison_df;allcols=true);println()
    end
end
