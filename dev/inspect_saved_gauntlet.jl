# Disposable manual inspector. Include from the REPL or IDE; the active project is unchanged.
# All comparisons below use saved numerical results; nothing is sampled or solved.
gauntlet_environment = normpath(joinpath(@__DIR__, "..", "gauntlet"))
gauntlet_environment in LOAD_PATH || push!(LOAD_PATH, gauntlet_environment)
using LineCableModels, DataFrames, Statistics, Measurements
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
# A failed include can leave the module name bound without its reader.
if !isdefined(@__MODULE__, :Gauntlet) || !isdefined(Gauntlet, :read_benchmark)
    include(joinpath(gauntlet_environment, "Gauntlet.jl"))
end

# Edit these and re-include.
campaign_directory = normpath(joinpath(@__DIR__, "..", "gauntlet", ".work", "all-references"))
benchmark_id = :benchmark_30kv_na2xs2y_630mm2_trefoil_fem
analysis_snapshot = nothing       # Optional explicit snapshot.jld2 path.
ydata = (R, X, G, B)

bands = (:all, :dc, :harmonic, :narrow, :wide)
detail_bands = (:all,)             # Print worst pairs for these bands; IDE tables retain all bands.
show_native_timings = true
blocks = nothing
problem_index = nothing          # Required when the result has several outer points.
plot_band = nothing
rms_metric = :relative            # Or :absolute.
make_plots = true
plot_backend = :gl                # Use :cairo for headless/inline plots.
display_plot = true
fig_size = (1400, 1000)

println("\nLoading saved benchmark: ", benchmark_id, " (no solver or MC run)")
path = analysis_snapshot === nothing ? joinpath(campaign_directory, string(benchmark_id)) : analysis_snapshot
benchmark = Gauntlet.read_benchmark(path; load_results=true, evidence=:numerical)
println("Building tables from retained results...")
requests = ydata
definition = BenchmarkTableDefinition(requests; bands)
benchmark_report = report(definition, (reference=benchmark.reference, candidate=benchmark.candidate,
    context=(id=benchmark.id,case_id=Symbol(first(benchmark.analyses)["case_id"]),collection=:manual),
    measurements=benchmark.measurements))
inspection_tables = benchmark_report.table

# These variables are ordinary DataFrames/collections available in the REPL and IDE.
feature_tables = filter(feature -> problem_index === nothing ||
    feature.problem_index == problem_index, inspection_tables.features)
isempty(feature_tables) && throw(ArgumentError("Selected parameter point is not present in this report"))
feature_dataframes = [getproperty(feature, rms_metric)
    for feature in feature_tables]
terms_df = filter(row -> problem_index === nothing || row.problem_index == problem_index, inspection_tables.terms)
maxima_df = filter(row -> problem_index === nothing || row.problem_index == problem_index, inspection_tables.maxima)
formulations_df = inspection_tables.formulations
formula_details_df = inspection_tables.formula_details
statistics_df = inspection_tables.statistics
sampling_df = inspection_tables.sampling
mean_sampling_precision_df = inspection_tables.mean_sampling_precision
execution_df = inspection_tables.execution
source_timings_df = inspection_tables.source_timings
performance_df = inspection_tables.performance
performance_comparison_df = inspection_tables.performance_comparison
performance_samples_df = inspection_tables.performance_samples
performance_environment_df = inspection_tables.performance_environment
performance_policy_df = inspection_tables.performance_policy

# Presentation only: labels, moments, RMS, units and timing scopes come from ReportBuilder.
# Full owner tables above remain available in the IDE through inspection_tables.
timing_labels = Dict(group.role[1] => (
    group.role[1] === :candidate && nrow(group) > 1 ?
    "Candidate formulation batch" :
    join(group.label, " / ")) for group in groupby(formulations_df, :role))
performance_display_df = isempty(performance_df) ? DataFrame() :
    select(performance_df, :role => ByRow(role -> timing_labels[role]) => :method,
        :median_seconds, :samples => :timed_calls, :allocated_MiB,
        :allocation_statistic, :allocation_scope, :scope, :reused)
timing_ratio_df = isempty(performance_comparison_df) ? DataFrame() :
    select(performance_comparison_df, :reference_over_candidate, :comparable)
execution_display_df = isempty(execution_df) ? DataFrame() :
    select(filter(row -> ismissing(row.point), execution_df),
        :role => ByRow(role -> timing_labels[role]) => :method,
        :seconds, :scope, :reused)
# Transpose backend-owned fields instead of making a window-wide row. In particular,
# accumulated worker times retain their scope; they are not campaign elapsed time.
source_timing_values_df = isempty(source_timings_df) ? DataFrame() :
    stack(source_timings_df, Not([:role, :method, :point]), [:method, :point];
        variable_name = :measurement, value_name = :value)

# Keep separate winners: the largest absolute error need not be at the largest
# relative error's terminal pair. Never put the two maxima beside a shared pair.
worst_relative_df = select(maxima_df,
    :method => :candidate,
    :quantity, :statistic, :problem_index => :point, :band,
    :relative_term_response => :response, :relative_term_excitation => :excitation,
    :maximum_relative_rms_percent => :RMS_percent, :compared, :unavailable)
worst_absolute_df = select(maxima_df,
    :method => :candidate,
    :quantity, :statistic, :problem_index => :point, :band,
    :absolute_term_response => :response, :absolute_term_excitation => :excitation,
    :maximum_absolute_rms => :RMS, :absolute_unit => :unit)
band_coverage_df = unique(select(maxima_df,
    :problem_index => :point, :band, :samples => :frequency_count,
    :actual_lower_Hz => :first_Hz, :actual_upper_Hz => :last_Hz))

println("\n", benchmark.id, " — saved calculations, current report")
println(join(formulations_df.label[formulations_df.role .=== :reference], " / "))
println("RMS tables show the worst eligible matrix entry in each band; entries are not averaged together.")
for (feature, frame) in zip(feature_tables, feature_dataframes)
    title = feature.statistic === :std ? "standard deviation (propagated uncertainty)" :
        feature.statistic === :mean ? "mean" : string(feature.statistic)
    println("\n", feature.quantity, " · ", title, " · parameter point ", feature.problem_index,
        " — maximum per-term ", rms_metric, " RMS ",
        rms_metric === :relative ? "[%]" : "[" * feature.absolute_unit * "]")
    show(stdout, MIME"text/plain"(), frame;
        allrows = true, allcols = true, truncate = 0)
    println()
end
for (label, frame) in (
        "Worst relative terms and counts" => worst_relative_df,
        "Worst absolute terms (independent maxima)" => worst_absolute_df)
    selected = filter(row -> row.band in detail_bands, frame)
    isempty(selected) && continue
    println("\n", label, " — ", join(string.(detail_bands), ", "))
    show(stdout, MIME"text/plain"(), selected; allrows = true, allcols = true, truncate = 0)
    println()
end
println("\nFrequency coverage")
show(stdout, MIME"text/plain"(), band_coverage_df; allrows = true, allcols = true)
println()
for (label, frame) in (
        "Recorded calculation wall times (not controlled repetitions)" => execution_display_df,
        "Controlled calculation measurements" => performance_display_df,
        "Recorded reference / candidate time ratio" => timing_ratio_df)
    isempty(frame) && continue
    println("\n", label)
    show(stdout, MIME"text/plain"(), frame; allrows = true, allcols = true, truncate = 0)
    println()
end
if isempty(performance_display_df)
    println("\nNo controlled performance measurements were saved for this case.")
else
    println("\nAllocated MiB are cumulative Julia allocations per call, not peak memory.")
    any(==(1), performance_display_df.timed_calls) &&
        println("A method has only one timed call; its timing variability cannot be assessed.")
    println("The ratio is reference time / candidate time for the recorded workloads, not per-formula cost.")
end
if show_native_timings && !isempty(source_timing_values_df)
    println("\nBackend timing records — scopes and counts as recorded by the backend")
    show(stdout, MIME"text/plain"(), source_timing_values_df;
        allrows = true, allcols = true, truncate = 0)
    println()
end
println("\nMissing relative RMS: near-zero operand or unavailable comparison; see terms_df.reason.")
println("All term errors and coordinates: terms_df. Full timing details: inspection_tables.")
for operand in (benchmark.reference,benchmark.candidate)
    isempty(operand.metadata.evidence_issues) || @warn "Auxiliary solver evidence differs; numerical payload was verified" issues=operand.metadata.evidence_issues
end

benchmark_plots = if make_plots
    plotting_environment=joinpath(@__DIR__,"plotting")
    plotting_environment in LOAD_PATH || push!(LOAD_PATH,plotting_environment)
    if plot_backend===:cairo
        @eval import CairoMakie
    elseif plot_backend===:gl
        @eval import GLMakie
    else
        throw(ArgumentError("plot_backend must be :gl or :cairo"))
    end
    LineCableModels.plot(benchmark_report,ydata;problem=problem_index,
        band=plot_band,blocks,backend=plot_backend,display_plot,fig_size,
        xscale=:log10,legend_position=:bottom)
else
    nothing
end
nothing
