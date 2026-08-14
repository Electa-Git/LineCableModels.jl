struct MCDistributionPlotSpec <: PlotBuilder.AbstractPlotSpec end

function _plot_exponent(series, field::Symbol)
    maximum_value = 0.0
    for item in series
        samples = field === :x ? item.xdata : item.ydata
        samples === nothing && continue
        for sample in samples
            nominal = abs(value(sample))
            nominal isa Real && isfinite(nominal) &&
                (maximum_value = max(maximum_value, Float64(nominal)))
        end
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

function _mc_selection(result::CableConstantsMC, quantity::Symbol, ijk)
    ijk === nothing || throw(ArgumentError("CableConstantsMC does not use matrix indices"))
    sample_values = has_samples(result) ? samples(result, quantity) : nothing
    distribution_value = has_distributions(result) ? distribution(result, quantity) :
                         nothing
    return sample_values, distribution_value, nothing
end

function _mc_selection(result::LineParametersMC, quantity::Symbol, ijk)
    selection = isnothing(ijk) ? (1, 1, 1) : ijk
    selection isa NTuple{3, Int} || throw(ArgumentError("ijk must be a tuple (i, j, k)"))
    i, j, k = selection
    sample_values = has_samples(result) ? samples(result, quantity, i, j, k) : nothing
    distribution_value = has_distributions(result) ?
                         distribution(result, quantity, i, j, k) : nothing
    return sample_values, distribution_value, selection
end

function _mc_target_unit(result, quantity::Symbol, length_unit, quantity_units)
    quantity in (:R, :L, :C, :G) ||
        throw(ArgumentError("quantity must be :R, :L, :C, or :G"))
    result isa CableConstantsMC && quantity === :G &&
        throw(
            ArgumentError("CableConstantsMC does not contain conductance"),
        )
    result_basis = result isa CableConstantsMC ? :per_length : basis(result)
    resolved = UnitHandler.line_component_unit(
        quantity,
        result_basis;
        length_unit,
        quantity_units
    )
    return resolved.quantity, resolved.units, resolved.scale
end

function _scaled_distribution(distribution_value::HistogramPDF, conversion)
    return HistogramPDF(
        distribution_value.edges .* conversion,
        distribution_value.density ./ conversion
    )
end

function _mc_plot_series(mode, data, values, model, bins, normalization)
    series = PlotBuilder.SeriesSpec[]
    if mode in (:hist, :pdf)
        if mode === :hist && data in (:samples, :both)
            values === nothing && throw(ArgumentError("samples were not retained"))
            push!(
                series,
                PlotBuilder.SeriesSpec(
                    :histogram,
                    values,
                    nothing,
                    nothing,
                    "samples";
                    attributes = (; bins, normalization)
                )
            )
        end
        if mode === :pdf || data in (:pdf, :both)
            model === nothing && throw(ArgumentError("distributions were not retained"))
            y = [model.density; last(model.density)]
            push!(
                series,
                PlotBuilder.SeriesSpec(
                    :stairs,
                    model.edges,
                    y,
                    nothing,
                    "model PDF";
                    attributes = (; step = :post, color = :red, linewidth = 2)
                )
            )
        end
    elseif mode === :ecdf
        values === nothing &&
            throw(ArgumentError("ECDF plotting requires retained samples"))
        model === nothing &&
            throw(ArgumentError("ECDF plotting requires retained distributions"))
        lower, upper = extrema(model.edges)
        padding = iszero(upper - lower) ? one(lower) : 0.05 * (upper - lower)
        x = collect(range(lower - padding, upper + padding; length = 500))
        empirical = StatsBase.ecdf(values)
        push!(
            series,
            PlotBuilder.SeriesSpec(
                :line,
                x,
                Distributions.cdf.(Ref(model), x),
                nothing,
                "model CDF";
                attributes = (; color = :red, linewidth = 2)
            )
        )
        push!(
            series,
            PlotBuilder.SeriesSpec(
                :line,
                x,
                empirical.(x),
                nothing,
                "empirical";
                attributes = (; color = :blue, linestyle = :dash, linewidth = 2)
            )
        )
    elseif mode === :qq
        values === nothing && throw(ArgumentError("Q-Q plotting requires retained samples"))
        model === nothing &&
            throw(ArgumentError("Q-Q plotting requires retained distributions"))
        sample_quantiles = sort(values)
        probabilities = ((1:length(values)) .- 0.5) ./ length(values)
        model_quantiles = Distributions.quantile.(Ref(model), probabilities)
        limits = extrema(vcat(sample_quantiles, model_quantiles))
        push!(
            series,
            PlotBuilder.SeriesSpec(
                :scatter,
                sample_quantiles,
                model_quantiles,
                nothing,
                "quantiles";
                attributes = (; color = :steelblue, markersize = 6)
            )
        )
        push!(
            series,
            PlotBuilder.SeriesSpec(
                :line,
                collect(limits),
                collect(limits),
                nothing,
                "perfect fit";
                attributes = (; color = :black, linestyle = :dash, linewidth = 2)
            )
        )
    else
        throw(ArgumentError("mode must be :hist, :pdf, :ecdf, or :qq"))
    end
    return series
end

function PlotBuilder.make_render(
        ::Type{MCDistributionPlotSpec},
        result::Union{CableConstantsMC, LineParametersMC};
        quantity::Symbol = :R,
        ijk = nothing,
        mode::Symbol = :hist,
        data::Symbol = :samples,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        nbins::Union{Nothing, Int} = nothing,
        normalization::Symbol = :none,
        fig_size::Tuple{Int, Int} = (800, 400),
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    PlotBuilder._validate_export_theme(export_theme)
    data in (:samples, :pdf, :both) || throw(
        ArgumentError("data must be :samples, :pdf, or :both"),
    )
    nbins === nothing || nbins > 0 || throw(ArgumentError("nbins must be positive"))
    sample_values, distribution_value, selection = _mc_selection(result, quantity, ijk)
    tag, target, conversion = _mc_target_unit(
        result,
        quantity,
        length_unit,
        quantity_units
    )
    scaled_values = sample_values === nothing ? nothing :
                    collect(sample_values) .* conversion
    scaled_distribution = distribution_value === nothing ? nothing :
                          _scaled_distribution(distribution_value, conversion)
    bin_count = isnothing(nbins) && scaled_values !== nothing ?
                _auto_nbins(scaled_values) : nbins
    bins = if scaled_distribution !== nothing
        scaled_distribution.edges
    elseif scaled_values !== nothing
        collect(fit(Histogram, scaled_values; nbins = bin_count, closed = :left).edges[1])
    else
        Float64[]
    end
    effective_normalization = data in (:pdf, :both) ? :pdf : normalization
    series = _mc_plot_series(
        mode,
        data,
        scaled_values,
        scaled_distribution,
        bins,
        effective_normalization
    )
    symbol = UnitHandler.get_symbol(tag)
    suffix = selection === nothing ? "" : "[$(join(selection, ','))]"
    title = mode === :hist ? "$symbol$suffix histogram" :
            mode === :pdf ? "$symbol$suffix probability density" :
            mode === :ecdf ? "$symbol$suffix cumulative distribution" :
            "$symbol$suffix Q-Q plot"
    x_label = mode === :qq ? "sample quantiles [$(UnitHandler.get_label(target))]" :
              "$(UnitHandler.get_label(tag)) [$(UnitHandler.get_label(target))]"
    y_label = mode === :hist ? String(effective_normalization) :
              mode === :pdf ? "density" :
              mode === :ecdf ? "cumulative probability" :
              "model quantiles [$(UnitHandler.get_label(target))]"
    xaxis = PlotBuilder.AxisSpec(:x, tag, target, x_label, :linear)
    yaxis = PlotBuilder.AxisSpec(
        :y,
        mode === :qq ? tag : UnitHandler.QuantityTag{:dimensionless}(),
        mode === :qq ? target : UnitHandler.Units(),
        y_label,
        :linear
    )
    view_spec = PlotBuilder.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        title,
        series,
        (; quantity, selection, mode)
    )
    page = PlotBuilder.PageSpec(
        title,
        fig_size,
        :single,
        PlotBuilder.ViewSpec[view_spec],
        (;
            quantity,
            selection,
            mode,
            data,
            x_exponent = _plot_exponent(series, :x),
            y_exponent = _plot_exponent(series, :y),
            export_theme,
            open_export,
            controls = PlotBuilder.control_definitions(xlog = false, ylog = false),
            configuration = (;
                quantity,
                selection,
                mode,
                data,
                length_unit,
                quantity_units,
                nbins,
                normalization
            )
        )
    )
    return PlotBuilder.RenderSpec(MCDistributionPlotSpec, PlotBuilder.PageSpec[page])
end
