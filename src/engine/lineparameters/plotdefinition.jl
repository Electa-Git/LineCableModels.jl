"""
$(TYPEDEF)

Prepare line-parameter observations for drawing by a loaded plotting extension.
"""
struct LineParameterPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Store aligned detached observations and panel state for one line dashboard.

$(TYPEDFIELDS)
"""
struct LineDashboardPayload{F, R, O, C, P, T, XS, YS, A, L, K, S}
    "Published frequency observation."
    frequency::F
    "Original scientific requests in display order."
    requests::R
    "Detached scientific observations."
    observations::O
    "Resolved row, column, and frequency indices."
    coordinates::C
    "Canvas positions in display order."
    positions::P
    "Panel titles in display order."
    titles::T
    "Admitted abscissa scales for each panel."
    xscales::XS
    "Admitted ordinate scales for each panel."
    yscales::YS
    "Axis visibility attributes for each panel."
    attributes::A
    "Presentation-only legend labels."
    legend_labels::L
    "Comparison colors, or an empty tuple for one source."
    colors::K
    "Captured runtime state used for current-state SVG replay."
    runtime::S
end

function _indices(selector, count::Int)
    selector === nothing && return collect(1:count)
    selector isa Colon && return collect(1:count)
    selector isa Integer && return [Int(selector)]
    selector isa AbstractRange && return collect(Int, selector)
    selector isa AbstractVector && return collect(Int, selector)
    throw(ArgumentError("conductor selections must be integers, ranges, vectors, `:`, or nothing"))
end

function _conductor_pairs(object, selector)
    row_count, column_count, _ = size(object)
    selector === nothing ||
        (selector isa Tuple && length(selector) == 2) ||
        throw(
            ArgumentError("conductor selection must be a tuple (rows, columns) or nothing"),
        )
    row_selector, column_selector = selector === nothing ? (nothing, nothing) : selector
    rows = _indices(row_selector, row_count)
    columns = _indices(column_selector, column_count)
    all(index -> index in 1:row_count, rows) || throw(BoundsError(object, rows))
    all(index -> index in 1:column_count, columns) || throw(BoundsError(object, columns))
    return [(i, j) for i in rows for j in columns]
end

_line_request_family(request) = Units.family(request_quantity(request))
_family_parent(::Val{:series}) = Z
_family_parent(::Val{:shunt}) = Y

function _line_input_defaults(frequencies)
    return (;
        frequencies,
        requests = (),
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip = true,
        title = nothing,
        labels = nothing,
        legend = nothing,
        xscale = :linear,
        yscale = :linear
    )
end

function PlotBuilder.entitle(
        ::Type{LineParameterPlotDefinition},
        source::Union{LineParameters, SeriesImpedance, ShuntAdmittance}
)
    return source
end
function PlotBuilder.input_defaults(::Type{LineParameterPlotDefinition}, ::LineParameters)
    _line_input_defaults(nothing)
end
function PlotBuilder.input_defaults(
        ::Type{LineParameterPlotDefinition},
        ::Union{SeriesImpedance, ShuntAdmittance}
)
    _line_input_defaults(nothing)
end
function PlotBuilder.renderer_defaults(
        ::Type{LineParameterPlotDefinition},
        ::Union{LineParameters, SeriesImpedance, ShuntAdmittance}
)
    (; fig_size = (800, 400))
end

function PlotBuilder.resolve(
        ::Type{LineParameterPlotDefinition},
        object,
        request::NamedTuple
)
    input = request.input
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    input.clip isa Bool || throw(ArgumentError("clip must be true or false"))
    object isa Union{SeriesImpedance, ShuntAdmittance} &&
        input.frequencies === nothing &&
        throw(
            ArgumentError("frequencies are required for SeriesImpedance and ShuntAdmittance"),
        )
    supplied_frequencies = input.frequencies === nothing ? nothing :
                           collect(input.frequencies)
    object isa Union{SeriesImpedance, ShuntAdmittance} &&
        length(supplied_frequencies) != size(object, 3) &&
        throw(DimensionMismatch(
            "frequency vector length does not match the parameter depth",
        ))
    requests = input.requests
    requests isa Tuple || throw(ArgumentError("requests must be a tuple"))
    isempty(requests) && throw(ArgumentError(
        "at least one explicit observable request is required",
    ))
    all(request -> request isa Tuple, requests) || throw(ArgumentError(
        "line plots accept explicit observable request tuples",
    ))
    identities = map(request_identity, requests)
    length(unique(identities)) == length(identities) || throw(ArgumentError(
        "line plots do not accept duplicate scientific quantities",
    ))
    validate_observables(object, requests)
    all(request -> length(request_indices(request)) == 3, requests) || throw(ArgumentError(
        "line plots require row, column, and frequency indices",
    ))
    if supplied_frequencies !== nothing
        all(isfinite, supplied_frequencies) ||
            throw(ArgumentError("frequencies must be finite"))
        input.xscale === :log10 && any(<=(0), supplied_frequencies) &&
            throw(
                DomainError(
                supplied_frequencies,
                "logarithmic frequency axes require positive frequencies"
            ),
            )
        any(request -> request_identity(request) in (L, C), requests) &&
            any(iszero, supplied_frequencies) &&
            throw(DomainError(
                supplied_frequencies,
                "inductance and capacitance are undefined at zero frequency"
            ))
    end
    request.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    all(>(0), request.renderer.fig_size) || throw(
        ArgumentError("fig_size dimensions must be positive"),
    )
    input.title === nothing || input.title isa AbstractString ||
        throw(
            ArgumentError("title must be a string or nothing"),
        )
    labels = input.labels
    labels === nothing || labels isa Tuple ||
        throw(
            ArgumentError("labels must be a tuple or nothing"),
        )
    labels === nothing || length(labels) == length(requests) ||
        throw(
            DimensionMismatch("labels must contain one entry per normalized request"),
        )
    labels === nothing || all(label -> label isa AbstractString, labels) ||
        throw(
            ArgumentError("labels must contain strings"),
        )
    input.legend === nothing || input.legend isa Tuple ||
        throw(
            ArgumentError("legend must be a tuple or nothing"),
        )
    return merge(
        request,
        (; input = merge(input, (; frequencies = supplied_frequencies)))
    )
end

function _frequency_payload(values, target)
    native = Units.units(:base, :hertz)
    factor = Units.scale_factor(native, target)
    return (
        values = map(value -> value * factor, values),
        quantity = Units.Quantity{:frequency}(),
        unit = target
    )
end

function _published_frequency(object, input, selector)
    target = Units.units(input.freq_unit, :hertz)
    if !(object isa LineParameters)
        return _frequency_payload(input.frequencies[selector], target)
    end
    published = only(observables(
        object,
        ((frequencies, selector),);
        units = (target,),
        clip = input.clip
    ))
    if input.frequencies !== nothing
        supplied = _frequency_payload(input.frequencies[selector], target)
        supplied.values == published.values || throw(
            ArgumentError("supplied frequencies do not match the LineParameters frequency axis"),
        )
    end
    return published
end

function _publish_request(object, request, target, clip::Bool)
    return only(observables(
        object,
        (request,);
        units = (target,),
        clip
    ))
end

function _request_coordinates(object, request)
    row, column, frequency = request_indices(request)
    dimensions = object isa LineParameters ? size(Z(object)) : size(object)
    rows = observation_indices(row, dimensions[1])
    columns = observation_indices(column, dimensions[2])
    frequency_count = object isa LineParameters ? nfrequencies(object) : size(object, 3)
    samples = observation_indices(frequency, frequency_count)
    return rows, columns, samples
end

function _materialized_line_request(object, input, request)
    rows, columns, samples = _request_coordinates(object, request)
    identity = observation_request(object, request).identity
    prefix = identity isa Tuple ? identity : (identity,)
    if object isa SeriesImpedance && identity === L
        prefix = (L, input.frequencies)
    elseif object isa ShuntAdmittance && identity === C
        prefix = (C, input.frequencies)
    end
    return (prefix..., rows, columns, samples)
end

function _publish_line_source(object, input, requests)
    coordinates = map(request -> _request_coordinates(object, request), requests)
    sample_indices = last.(coordinates)
    all(==(first(sample_indices)), sample_indices) || throw(DimensionMismatch(
        "all requests on one line dashboard must select the same frequency indices",
    ))
    frequency = _published_frequency(object, input, first(sample_indices))
    targets = unit_targets(
        requests,
        basis(object);
        length_prefix = input.length_unit,
        overrides = input.quantity_units
    )
    observations = map(requests, targets) do request, target
        _publish_request(
            object,
            _materialized_line_request(object, input, request),
            target,
            input.clip
        )
    end
    all(observation -> size(observation.values, 3) == length(frequency.values), observations) ||
        throw(
            DimensionMismatch("frequency count does not match line-parameter samples"),
        )
    legend = input.legend
    if legend !== nothing
        matrix_pairs = map(coordinates) do (rows, columns, _)
            Tuple((row, column) for row in rows for column in columns)
        end
        all(==(first(matrix_pairs)), matrix_pairs) || throw(DimensionMismatch(
            "legend overrides require every request to select the same matrix indices",
        ))
        length(legend) == length(first(matrix_pairs)) || throw(DimensionMismatch(
            "legend must contain one label per selected matrix element",
        ))
        all(label -> label isa AbstractString, legend) || throw(ArgumentError(
            "legend must contain strings",
        ))
    end
    length(frequency.values) <= 1 &&
        @warn "Frequency vector has $(length(frequency.values)) sample(s); nothing to plot."
    return (; frequency, observations, coordinates)
end

function _supports_log_values(samples)
    found = false
    samples === nothing && return false
    for sample in samples
        found = true
        value = nominal(sample)
        uncertainty_value = abs(uncertainty(sample))
        value isa Real && isfinite(value) && isfinite(uncertainty_value) &&
        value - uncertainty_value > 0 || return false
    end
    return found
end

_axis_scales(values) = _supports_log_values(values) ? (:linear, :log10) : (:linear,)

function _panel_positions(count::Int)
    columns = max(1, ceil(Int, sqrt(count)))
    return Tuple(((index - 1) ÷ columns + 1, (index - 1) % columns + 1)
    for index in 1:count)
end

function _line_page(configuration, published, family, page_index::Int)
    selected = Tuple(
        (index, scientific_request, observation, coordinates)
    for (index, (scientific_request, observation, coordinates)) in enumerate(zip(
            configuration.input.requests,
            published.observations,
            published.coordinates
        ))
    if _line_request_family(scientific_request) === family
    )
    positions = _panel_positions(length(selected))
    titles = map(selected) do selected_request
        request_index, scientific_request, observation, coordinates = selected_request
        return configuration.input.labels === nothing ? Units.label(observation.quantity) :
               String(configuration.input.labels[request_index])
    end
    xscales = map(selected) do _
        admitted = _axis_scales(published.frequency.values)
        configuration.input.xscale in admitted || throw(DomainError(
            published.frequency.values,
            "logarithmic frequency axes require positive finite data and uncertainty bounds"
        ))
        return admitted
    end
    yscales = map(selected) do selected_request
        observation = selected_request[3]
        admitted = _axis_scales(observation.values)
        configuration.input.yscale in admitted || throw(DomainError(
            observation.values,
            "logarithmic ordinate axes require positive finite data and uncertainty bounds"
        ))
        return admitted
    end
    parent = _family_parent(family)
    title = configuration.input.title === nothing ?
            Units.label(Units.quantity(parent)) : String(configuration.input.title)
    key = (;
        page = page_index,
        family,
        requests = getindex.(selected, 2)
    )
    payload = LineDashboardPayload(
        published.frequency,
        getindex.(selected, 2),
        getindex.(selected, 3),
        getindex.(selected, 4),
        positions,
        titles,
        xscales,
        yscales,
        ntuple(_ -> (;), length(selected)),
        configuration.input.legend,
        (),
        (;
            xscale = configuration.input.xscale,
            yscale = configuration.input.yscale,
            panels = nothing
        )
    )
    return PlotBuilder.PlotPage(
        title,
        configuration.renderer.fig_size,
        key,
        payload;
        legend = PlotBuilder.LegendDefinition(),
        export_definition = PlotBuilder.ExportDefinition(
            theme = configuration.renderer.export_theme,
            name = title,
            open_file = configuration.renderer.open_export
        )
    )
end

_line_families(::SeriesImpedance) = (Val(:series),)
_line_families(::ShuntAdmittance) = (Val(:shunt),)
_line_families(::LineParameters) = (Val(:series), Val(:shunt))

function PlotBuilder.fetch(
        ::Type{LineParameterPlotDefinition},
        object,
        request::NamedTuple
)
    published = _publish_line_source(object, request.input, request.input.requests)
    length(published.frequency.values) <= 1 && return PlotBuilder.PlotPage[]
    families = Tuple(family
    for family in _line_families(object)
    if any(item -> _line_request_family(item) === family, request.input.requests))
    return PlotBuilder.PlotPage[_line_page(request, published, family, page_index)
                                for (page_index, family) in enumerate(families)]
end
