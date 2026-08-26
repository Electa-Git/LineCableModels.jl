using Test
using LineCableModels
using LineCableModels.PlotBuilder: axis!, plotwindow, register!
using GLMakie

const MONTE_CARLO_GALLERY_SMOKE_ONLY = lowercase(
    get(ENV, "LINECABLEMODELS_GL_GALLERY_SMOKE", "false")
) == "true"

set_backend!(:gl)

# The fixture is a completed Monte Carlo result. Every chart below observes
# retained products; none of the plotting calls reruns the calculation.
sample_values = collect(range(1.0, 5.0; length = 81))
summary = SampleSummary(sample_values)
histogram = HistogramDensity(
    collect(range(1.0, 5.0; length = 9)),
    fill(0.25, 8)
)
formulation = MonteCarlo(
    Formulation();
    trials = length(sample_values),
    seed = 41,
    return_samples = true,
    return_histograms = true
)
cable_result = MonteCarloResult(
    formulation,
    [CableConstants(3.0, 3.0, 3.0)],
    [(R = summary, L = summary, C = summary)],
    [(R = copy(sample_values), L = copy(sample_values), C = copy(sample_values))],
    [(R = histogram, L = histogram, C = histogram)],
    UInt64(41),
    UInt64[41],
    [length(sample_values)]
)

display_plot = !MONTE_CARLO_GALLERY_SMOKE_ONLY
gallery = Pair{String, UIPlot}[]

# Native Makie function identity chooses the primitive. The LineCableModels
# overload first publishes the requested marginal, then creates the standard
# shell and calls hist! on a managed axis.
push!(gallery, "Sample histogram" => Makie.hist(
    cable_result,
    R;
    bins = 12,
    normalization = :pdf,
    quantity_units = :milli,
    backend = :gl,
    display_plot
))

# Model-only charts use the retained HistogramDensity. If histograms were not
# retained, stairs and lines derive the same model from retained samples.
push!(gallery, "Model probability density" => Makie.stairs(
    cable_result,
    R;
    quantity_units = :milli,
    color = :red,
    backend = :gl,
    display_plot
))
push!(gallery, "Empirical cumulative distribution" => Makie.ecdfplot(
    cable_result,
    R;
    quantity_units = :milli,
    color = :blue,
    backend = :gl,
    display_plot
))
push!(gallery, "Model cumulative distribution" => Makie.lines(
    cable_result,
    R;
    quantity_units = :milli,
    color = :red,
    backend = :gl,
    display_plot
))
push!(gallery, "Sample/model Q-Q" => Makie.qqplot(
    cable_result,
    R;
    quantity_units = :milli,
    qqline = :identity,
    backend = :gl,
    display_plot
))

# Composition remains inside the same shell. This callback publishes samples
# and the model once, creates one managed axis, invokes Makie's mutating
# primitives directly, and registers both groups for the common legend,
# visibility, reset, formatter, status, and SVG-replay machinery.
combined = plotwindow(
    title = "Samples and model density",
    backend = :gl,
    display_plot = display_plot
) do context
    published = observables(
        cable_result,
        (
            sample = (samples, R, 1, Colon()),
            model = (histograms, R, 1)
        )
    )
    sample = published.sample
    model = published.model
    axis = axis!(
        context,
        context.canvas[1, 1],
        sample,
        nothing;
        title = "R samples and density"
    )
    bars = hist!(axis, sample.values; bins = 12, normalization = :pdf)
    curve = stairs!(
        axis,
        model.values.edges,
        [model.values.density; last(model.values.density)];
        step = :post,
        color = :red
    )
    register!(
        context,
        axis;
        groups = (samples = (bars,), model = (curve,)),
        labels = (samples = "samples", model = "model PDF"),
        data = (
            (xdata = sample.values, ydata = nothing,
                group = :samples, label = "samples"),
            (xdata = model.values.edges,
                ydata = [model.values.density; last(model.values.density)],
                group = :model, label = "model PDF")
        )
    )
end
push!(gallery, "Samples and model PDF" => combined)

# Matrix-valued results use the observable request itself for the conductor and
# frequency coordinate. There is no plotting-specific ijk keyword.
frequencies = [50.0, 500.0, 5_000.0]
frequency_scales = [0.8, 1.0, 1.4]
line_samples = reshape(
    [scale * sample_values[trial]
     for scale in frequency_scales, trial in eachindex(sample_values)],
    1,
    1,
    length(frequencies),
    length(sample_values)
)
summaries = reshape(
    [SampleSummary(vec(line_samples[1, 1, index, :]))
     for index in eachindex(frequencies)],
    1,
    1,
    length(frequencies)
)
histograms = reshape(
    [HistogramDensity(histogram.edges .* scale, histogram.density ./ scale)
     for scale in frequency_scales],
    1,
    1,
    length(frequencies)
)
line_result = MonteCarloResult(
    formulation,
    [LineParameters(
        zeros(ComplexF64, 1, 1, length(frequencies)),
        zeros(ComplexF64, 1, 1, length(frequencies)),
        frequencies
    )],
    [(R = summaries, L = summaries, C = summaries, G = summaries)],
    [(R = line_samples, L = line_samples, C = line_samples, G = line_samples)],
    [(R = histograms, L = histograms, C = histograms, G = histograms)],
    UInt64(41),
    UInt64[41],
    [length(sample_values)]
)
line_request = @observe R[1, 1, 2]
push!(gallery, "Line parameters R[1,1,2]" => Makie.hist(
    line_result,
    line_request;
    bins = 12,
    backend = :gl,
    display_plot
))

@testset "manual GL Monte Carlo gallery" begin
    @test length(gallery) == 7
    @test all(pair -> pair.second isa UIPlot, gallery)
    @test all(pair -> pair.second.context.backend === :gl, gallery)
    @test all(pair -> !haskey(pair.second.controls, :xlog), gallery)
    @test all(pair -> !haskey(pair.second.controls, :ylog), gallery)
end

println("Built seven shell-managed Monte Carlo inspection pages.")

if MONTE_CARLO_GALLERY_SMOKE_ONLY
    @test all(pair -> pair.second.context.window === nothing, gallery)
    println("GL gallery smoke-only mode complete; native windows were not opened.")
    exit()
end

println("Close all windows to finish, or press Ctrl+C here.")
try
    while any(pair -> isopen(pair.second.context.window), gallery)
        sleep(0.1)
    end
finally
    GLMakie.closeall()
end
