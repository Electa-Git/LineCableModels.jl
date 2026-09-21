# Disposable manual inspector. Include from the REPL or IDE; the active project is unchanged.
# Recompare saved numerical operands explicitly; nothing is sampled or solved.
# Before: report construction mixed reading and comparison. Now the raw-operand
# report convenience completes comparison, constructs ObservedResult, then tabulates.
# Existing current snapshots can instead be selected with analysis_snapshot.
gauntlet_environment = normpath(joinpath(@__DIR__, "..", "gauntlet"))
gauntlet_environment in LOAD_PATH || push!(LOAD_PATH, gauntlet_environment)
using LineCableModels, DataFrames, Statistics, Measurements
using LineCableModels.ReportBuilder: BenchmarkTableDefinition
# A failed include can leave the module name bound without its reader.
if !isdefined(@__MODULE__, :Gauntlet) || !isdefined(Gauntlet, :read_benchmark)
    include(joinpath(gauntlet_environment, "Gauntlet.jl"))
end

include(joinpath(@__DIR__, "inspect_saved_gauntlet_inputs.jl"))

# Edit these and re-include.
campaign_directory = normpath(joinpath(
    @__DIR__, "..", "gauntlet", ".work", "all-references"))
benchmark_id = :benchmark_220kv_milliken_1x2500_252_trefoil_pscad
analysis_snapshot = nothing       # Optional CURRENT observed snapshot.jld2.
use_previous_complete = false    # Explicitly opt into an older completed attempt.
ydata = (R, X, G, B)

bands = (:all,)
detail_bands = ()                  # Optional worst-pair detail, e.g. (:all,); all bands remain in the IDE.
show_native_timings = false
show_full_report = false          # Full-width scientific descriptions can be very long.
series_labels = nothing           # Optional candidate labels, then the reference label.
# Before: blocks controlled pagination. Now layout is the sole panel capacity;
# each requested quantity/statistic produces its own figure family.
layout = nothing
problem_index = nothing          # Optional original physical-point index.
plot_band = nothing
rms_metric = :relative            # Or :absolute.
make_plots = true
plot_backend = :gl                # Use :cairo for headless/inline plots.
enable_svg_export = false         # Also load CairoMakie for the SVG toolbar button.
errorbar_sampling = :staggered     # Or :all to inspect every retained interval.
display_plot = true
fig_size = (1400, 1000)

println("\nLoading saved benchmark: ", benchmark_id, " (no solver or MC run)")
benchmark = analysis_snapshot === nothing ?
    read_inspection_operands(campaign_directory,benchmark_id;previous=use_previous_complete) :
    Gauntlet.read_benchmark(analysis_snapshot;load_results=true,evidence=:numerical)
println(analysis_snapshot===nothing ? "Comparing saved operands and constructing current observations..." : "Tabulating retained observations...")
requests = ydata
definition = BenchmarkTableDefinition(requests; bands)
benchmark_report = if analysis_snapshot===nothing
    # Fresh comparisons/observations from checksummed numerical inputs; no saved
    # analysis is modified and no old comparison-generation adapter is installed.
    report(definition,(reference=benchmark.reference,candidate=benchmark.candidate,
        measurements=benchmark.measurements,context=(id=benchmark_id,));
        requests=Tuple(unique((ydata...,requests...))))
else
    report(definition,benchmark)
end
inspection_tables = benchmark_report.tables

# Before: historical metadata expanded every numerical setting into plot legends.
# Use explicit input identities for this manual view; the complete owner-provided
# descriptions and their label mapping remain in plot_series_df for inspection.
plot_points = benchmark_report.observed isa ObservedResult ? [benchmark_report.observed] : collect(benchmark_report.observed)
plot_labels = ["Candidate $index" for index in eachindex(plot_points)]
if benchmark_report.reference !== nothing
    push!(plot_points, benchmark_report.reference)
    push!(plot_labels, "Reference")
end
series_labels === nothing || (plot_labels = collect(series_labels))
plot_series_df = DataFrame(label=plot_labels,
    description=LineCableModels.Grammar.observation_labels(plot_points))
println("Plot labels and complete scientific descriptions: plot_series_df.")
show(stdout, MIME"text/plain"(), plot_series_df; allrows=true, truncate=100)
println()

# These variables are ordinary DataFrames/collections available in the REPL and IDE.
feature_tables = filter(
    feature -> problem_index === nothing ||
               feature.problem_index == problem_index,
    inspection_tables.features)
isempty(feature_tables) &&
    throw(ArgumentError("Selected parameter point is not present in this report"))
feature_dataframes = [getproperty(feature, rms_metric)
                      for feature in feature_tables]
terms_df = filter(row -> problem_index === nothing || row.problem_index == problem_index,
    inspection_tables.terms)
maxima_df = filter(row -> problem_index === nothing || row.problem_index == problem_index,
    inspection_tables.maxima)
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

# Compact display tables are owned by ReportBuilder, shared with documentation.
overview_tables = inspection_tables.overview
performance_display_df = overview_tables.performance
timing_ratio_df = overview_tables.timing_ratio
execution_display_df = overview_tables.execution
source_timing_values_df = overview_tables.source_timings

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
band_coverage_df = overview_tables.coverage

println("\n", benchmark.id, " — saved calculations, current report")
if show_full_report
    show(stdout, MIME"text/plain"(), benchmark_report;
        metric = rms_metric, problem = problem_index, native_timings = show_native_timings)
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
println("\nMissing relative RMS: near-zero operand or unavailable comparison; see terms_df.reason.")
println("All term errors and coordinates: terms_df. Full timing details: inspection_tables.")
for operand in (benchmark.reference, benchmark.candidate)
    isempty(operand.metadata.evidence_issues) ||
        @warn "Auxiliary solver evidence differs; numerical payload was verified" issues=operand.metadata.evidence_issues
end

benchmark_plots = if make_plots
    plotting_environment=joinpath(@__DIR__, "plotting")
    plotting_environment in LOAD_PATH || push!(LOAD_PATH, plotting_environment)
    if enable_svg_export
        @eval import CairoMakie
    end
    if plot_backend===:cairo
        @eval import CairoMakie
    elseif plot_backend===:gl
        @eval import GLMakie
    else
        throw(ArgumentError("plot_backend must be :gl or :cairo"))
    end
    LineCableModels.plot(benchmark_report, ydata; problem = problem_index,
        band = plot_band, layout, backend = plot_backend, display_plot, fig_size,
        xscale = :log10, legend_position = :bottom, series_labels=plot_labels, errorbar_sampling)
else
    nothing
end
nothing
