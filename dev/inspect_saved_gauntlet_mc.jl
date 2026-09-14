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
campaign_directory = normpath(joinpath(
    @__DIR__, "..", "gauntlet", ".work", "all-references"))
benchmark_id = :benchmark_30kv_na2xs2y_630mm2_trefoil_lep_montecarlo
analysis_snapshot = nothing       # Optional explicit snapshot.jld2 path.
ydata = (R, L, G, C)
plot_statistics = (std, mean)      # Statistics compared in the RMS tables.
make_statistic_plots = false      # Optional separate mean-only/std-only plots, in addition to mean ± std.
inspection_term = (1, 1)          # Row/column for a small table of retained statistics.
inspection_frequency_Hz = 50.0    # Preview uses the nearest saved frequency, never interpolation.
bands = (:all, :dc, :harmonic, :narrow, :wide)
detail_bands = ()                  # Optional worst-pair detail, e.g. (:all,); all bands remain in the IDE.
show_native_timings = false
blocks = nothing
problem_index = nothing          # Required when the result has several outer points.
plot_band = nothing
rms_metric = :relative            # Or :absolute.
make_plots = true
plot_backend = :gl                # Use :cairo for headless/inline plots.
display_plot = true
fig_size = (1400, 1000)

println("\nLoading saved benchmark: ", benchmark_id, " (no solver or MC run)")
path = analysis_snapshot === nothing ? joinpath(campaign_directory, string(benchmark_id)) :
       analysis_snapshot
benchmark = Gauntlet.read_benchmark(path; load_results = true, evidence = :numerical)
println("Building tables from retained results...")
requests = Tuple((statistics, q, statistic) for statistic in plot_statistics for q in ydata)
definition = BenchmarkTableDefinition(requests; bands)
benchmark_report = report(definition,
    (reference = benchmark.reference, candidate = benchmark.candidate,
        context = (
            id = benchmark.id, case_id = Symbol(first(benchmark.analyses)["case_id"]),
            collection = :manual),
        measurements = benchmark.measurements))
inspection_tables = benchmark_report.table

# These variables are ordinary DataFrames/collections available in the REPL and IDE.
feature_tables = filter(
    feature -> problem_index === nothing ||
               feature.problem_index == problem_index,
    inspection_tables.features)
isempty(feature_tables) &&
    throw(ArgumentError("Selected parameter point is not present in this report"))
configurations = unique(feature.problem_index for feature in feature_tables)
length(configurations) == 1 || throw(ArgumentError(
    "Set problem_index to select one parameter configuration for this manual inspector"))
configuration_index = only(configurations)
multiple_configurations = length(unique(inspection_tables.maxima.problem_index)) > 1
feature_dataframes = [getproperty(feature, rms_metric)
                      for feature in feature_tables]
std_error_dataframes = [frame
                        for (feature, frame) in zip(feature_tables, feature_dataframes)
                        if feature.statistic === :std]
mean_error_dataframes = [frame
                         for (feature, frame) in zip(feature_tables, feature_dataframes)
                         if feature.statistic === :mean]
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
# This manual overlay selects one method per role; use its owned label unchanged.
timing_labels = Dict(role => only(formulations_df.label[formulations_df.role .=== role])
for role in (:reference, :candidate))

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
# Preview existing scalar publications at one explicit matrix coordinate.
sampling_display_df = overview_tables.sampling
cdf_precision_df = overview_tables.cdf_precision
statistics_values_df = DataFrame()
statistics_preview_df = DataFrame()
mean_values_df = DataFrame()
std_values_df = DataFrame()
comparison_dataframes = Dict{Symbol, DataFrame}()
mean_precision_preview_df = DataFrame()
if !isempty(statistics_df)
    preview_point = configuration_index
    point_statistics = filter(
        row -> row.point == preview_point &&
               (row.row, row.column) == inspection_term,
        statistics_df)
    isempty(point_statistics) &&
        throw(ArgumentError("No retained statistics for the selected point/term"))
    saved_frequencies = unique(point_statistics.frequency)
    preview_frequency_Hz = nominal(saved_frequencies[argmin(abs.(saved_frequencies .-
                                                                 inspection_frequency_Hz))])
    quantity_columns = unique(feature.quantity for feature in feature_tables)
    # Publications can use mH or µF even with a per-metre length basis. Convert
    # with their owned unit metadata so the preview and mean standard errors
    # use the same units; do not relabel display-scaled values as native values.
    column_contracts = LineCableModels.ReportBuilder.observation_columns(statistics_df)
    native_units = Dict(quantity => LineCableModels.Units.native_unit(
                            column_contracts[quantity].quantity, metadata(statistics_df, "basis"))
    for quantity in quantity_columns)
    statistics_values_df = select(statistics_df,
        :role => ByRow(role -> timing_labels[role]) => :method,
        :point, :frequency =>
            ByRow(value -> iszero(uncertainty(value)) ? nominal(value) : value) =>
                :frequency_Hz,
        :row, :column, :statistic,
        (quantity =>
             (values -> values .* LineCableModels.Units.scale_factor(
                 column_contracts[quantity].unit, native_units[quantity])) =>
                 Symbol(string(quantity, " [", LineCableModels.Units.label(native_units[quantity]), "]"))
        for quantity in quantity_columns)...)
    statistics_preview_df = filter(
        row -> row.point == preview_point &&
               (row.row, row.column) == inspection_term &&
               row.frequency_Hz == preview_frequency_Hz,
        statistics_values_df)
    multiple_configurations || select!(statistics_preview_df, Not(:point))
    mean_values_df = select(filter(:statistic => ==(:mean), statistics_preview_df), Not(:statistic))
    std_values_df = select(filter(:statistic => ==(:std), statistics_preview_df), Not(:statistic))
    # One row per frequency, both methods side by side. No frequency aggregation
    # or statistical estimation is done here: this only pivots the owned table.
    selected_moments = filter(
        row -> row.point == preview_point &&
               (row.row, row.column) == inspection_term && row.statistic in (:mean, :std),
        statistics_values_df)
    for quantity in quantity_columns
        unit = LineCableModels.Units.label(native_units[quantity])
        column = Symbol(string(quantity, " [", unit, "]"))
        frame = select(selected_moments, :frequency_Hz, :row, :column,
            [:method, :statistic] =>
                ByRow((method, stat) -> string(method, " · ", stat, " [", unit, "]")) =>
                    :series,
            column => :value)
        comparison_dataframes[quantity] = sort(
            unstack(frame,
                [:frequency_Hz, :row, :column], :series, :value), :frequency_Hz)
    end
    if !isempty(mean_sampling_precision_df)
        mean_precision_preview_df = select(
            filter(
                row -> row.point == preview_point &&
                       (row.row, row.column) == inspection_term &&
                       row.frequency_Hz == preview_frequency_Hz,
                mean_sampling_precision_df),
            :method, :frequency_Hz, :row, :column, :quantity, :mean_standard_error, :unit)
    end
    println("\nRetained statistics",
        multiple_configurations ? " — configuration $preview_point" : "",
        " — matrix entry ", inspection_term, ", ", preview_frequency_Hz,
        " Hz (requested ", inspection_frequency_Hz, " Hz)")
    println("These are values at this frequency, not maxima or band averages.")
    for (label, frame) in ("Mean — mean_values_df" => mean_values_df, "Std — std_values_df" =>
        std_values_df)
        println("\n", label)
        show(
            stdout, MIME"text/plain"(), frame; allrows = true, allcols = true, truncate = 0)
        println()
    end
    println("\nFull-frequency comparisons for entry ", inspection_term,
        ": comparison_dataframes, keyed by ", join(string.(quantity_columns), ", "), ".")
end
# Direct IDE table bindings; unselected quantities are nothing.
R_values_df = get(comparison_dataframes, :R, nothing)
X_values_df = get(comparison_dataframes, :X, nothing)
L_values_df = get(comparison_dataframes, :L, nothing)
G_values_df = get(comparison_dataframes, :G, nothing)
B_values_df = get(comparison_dataframes, :B, nothing)
C_values_df = get(comparison_dataframes, :C, nothing)

show(stdout, MIME"text/plain"(), benchmark_report;
    metric = rms_metric, problem = problem_index, native_timings = show_native_timings)
println()
for (label, frame) in (
    "Worst relative terms and counts" => worst_relative_df,
    "Worst absolute terms (independent maxima)" => worst_absolute_df)
    selected = filter(row -> row.band in detail_bands, frame)
    multiple_configurations || select!(selected, Not(:point))
    isempty(selected) && continue
    println("\n", label, " — ", join(string.(detail_bands), ", "))
    show(stdout, MIME"text/plain"(), selected; allrows = true, allcols = true, truncate = 0)
    println()
end
println("\nMissing relative RMS: near-zero operand or unavailable comparison; see terms_df.reason.")
println("All term errors and coordinates: terms_df. Full timing details: inspection_tables.")
for (label, frame) in (
    "MC standard error of the mean — same point, entry and frequency" =>
    mean_precision_preview_df,)
    isempty(frame) && continue
    println("\n", label)
    visible = !multiple_configurations && :point in propertynames(frame) ?
              select(frame, Not(:point)) : frame
    show(stdout, MIME"text/plain"(), visible; allrows = true, allcols = true, truncate = 0)
    println()
end
println("\nStd compares propagated physical uncertainty. MC standard error measures sampling noise in the mean.")
println("The CDF bound is not an error bar on mean or std. MC std sampling precision is not retained.")
println("Full scalar statistics with units: statistics_values_df. MC mean standard errors: mean_sampling_precision_df.")
for operand in (benchmark.reference, benchmark.candidate)
    isempty(operand.metadata.evidence_issues) ||
        @warn "Auxiliary solver evidence differs; numerical payload was verified" issues=operand.metadata.evidence_issues
end

mean_std_plots = nothing
statistic_plots = nothing
benchmark_plots = if make_plots
    plotting_environment=joinpath(@__DIR__, "plotting")
    plotting_environment in LOAD_PATH || push!(LOAD_PATH, plotting_environment)
    if plot_backend===:cairo
        @eval import CairoMakie
    elseif plot_backend===:gl
        @eval import GLMakie
    else
        throw(ArgumentError("plot_backend must be :gl or :cairo"))
    end
    all(quantity -> quantity in (R, X, L, G, B, C), ydata) || throw(ArgumentError(
        "Use R/X/L and G/B/C for marginal mean ± std plots; complex-magnitude uncertainty requires joint statistics"))
    # The UQ owner reads one uncertainty-bearing core: MC reconstructs marginal
    # means/stds; LEP preserves its native values and Measurement dependencies.
    uncertainty_sources = map((benchmark.reference, benchmark.candidate)) do operand
        uncertain(operand.result, configuration_index)
    end
    if plot_band !== nothing
        comparisons = filter(row -> LineCableModels.details(row.error).band == plot_band,
            benchmark_report.published.comparisons)
        isempty(comparisons) && throw(ArgumentError("plot_band must be included in bands"))
        samples = only(unique(LineCableModels.details(row.error).indices
        for row in comparisons))
        uncertainty_sources = map(source -> source[samples], uncertainty_sources)
    end
    println("\nPlotting both methods: frequency on x, mean on y, error bars ±1 std (not standard error).")
    mean_std_plots = LineCableModels.plot(uncertainty_sources...; ydata,
        series_labels = (timing_labels[:reference], timing_labels[:candidate]),
        series_attributes = (
            (marker = :circle, markersize = 5), (marker = :utriangle, markersize = 5)),
        title = string(benchmark.id, " — mean ± 1σ"),
        blocks, backend = plot_backend, display_plot, fig_size,
        xscale = :log10, signed_ylog = true, legend_position = :bottom)
    if make_statistic_plots
        statistic_plots = LineCableModels.plot(
            benchmark_report, requests; problem = problem_index,
            band = plot_band, blocks, backend = plot_backend, display_plot, fig_size,
            xscale = :log10, legend_position = :bottom)
    end
    mean_std_plots
else
    nothing
end
nothing
