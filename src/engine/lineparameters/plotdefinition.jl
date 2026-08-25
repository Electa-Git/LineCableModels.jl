"""
$(TYPEDEF)

Prepare line-parameter observations for drawing by a loaded plotting extension.
"""
struct LineParameterPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Store one line in a detached line-parameter panel.

$(TYPEDFIELDS)
"""
struct LineCurve{Y, S}
    "Displayed ordinate values in the unit recorded by the panel observation."
    values::Y
    "Legend text."
    label::String
    "Identity shared by lines controlled by one legend entry."
    group::Symbol
    "Line attributes owned by the scientific plot definition."
    style::S
end

"""
$(TYPEDEF)

Store the observations, curves, and axis behavior for one line-parameter panel.

$(TYPEDFIELDS)
"""
struct LinePanelPayload{R, X, Y, C, XS <: Tuple, YS <: Tuple, A}
    "Function-valued scientific request represented by the panel."
    request::R
    "One-based row and column in the page canvas."
    position::Tuple{Int, Int}
    "Displayed panel title."
    title::String
    "Published abscissa observation."
    x_observation::X
    "Published ordinate observation restricted to the displayed curves."
    y_observation::Y
    "Detached curves in draw order."
    curves::C
    "Initial abscissa scale."
    xscale::Symbol
    "Initial ordinate scale."
    yscale::Symbol
    "Scales admitted for the abscissa."
    xscales::XS
    "Scales admitted for the ordinate."
    yscales::YS
    "Axis visibility attributes."
    attributes::A
end

"""
$(TYPEDEF)

Store one detached line-parameter page for the standard plotting shell.

$(TYPEDFIELDS)
"""
struct LinePagePayload{K, P, E}
    "Displayed page title."
    title::String
    "Scientific and placement identity of the page."
    key::K
    "Detached panels in draw order."
    panels::P
    "Legend behavior supplied to the standard shell."
    legend::PlotBuilder.LegendDefinition
    "SVG export behavior supplied to the standard shell."
    export_definition::E
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

const _DISPLAY_ZERO_RTOL = sqrt(eps(Float64))

_request_parent(::typeof(R)) = Z
_request_parent(::typeof(X)) = Z
_request_parent(::typeof(L)) = Z
_request_parent(::typeof(G)) = Y
_request_parent(::typeof(B)) = Y
_request_parent(::typeof(C)) = Y
_request_parent(request::Tuple) = first(request)

_default_line_requests(::LineParameters) = (R, X, G, B)
_default_line_requests(::SeriesImpedance) = (R, X)
_default_line_requests(::ShuntAdmittance) = (G, B)

function _line_requests(::SeriesImpedance, accessor)
    accessor === R && return (R,)
    accessor === X && return (X,)
    accessor === L && return (L,)
    accessor === real && return (R,)
    accessor === imag && return (X,)
    accessor === abs && return ((Z, abs),)
    accessor === angle && return ((Z, angle),)
    accessor === Z && return (R, X)
    throw(ArgumentError("accessor $(accessor) is not defined for SeriesImpedance presentation"))
end

function _line_requests(::ShuntAdmittance, accessor)
    accessor === G && return (G,)
    accessor === B && return (B,)
    accessor === C && return (C,)
    accessor === real && return (G,)
    accessor === imag && return (B,)
    accessor === abs && return ((Y, abs),)
    accessor === angle && return ((Y, angle),)
    accessor === Y && return (G, B)
    throw(ArgumentError("accessor $(accessor) is not defined for ShuntAdmittance presentation"))
end

function _line_requests(::LineParameters, accessor)
    accessor === Z && return (R, X)
    accessor === Y && return (G, B)
    accessor === real && return (R, G)
    accessor === imag && return (X, B)
    accessor === abs && return ((Z, abs), (Y, abs))
    accessor === angle && return ((Z, angle), (Y, angle))
    accessor === R && return (R,)
    accessor === X && return (X,)
    accessor === L && return (L,)
    accessor === G && return (G,)
    accessor === B && return (B,)
    accessor === C && return (C,)
    throw(ArgumentError("accessor $(accessor) is not defined for LineParameters presentation"))
end

function _resolve_line_requests(object, quantities)
    quantities isa Tuple || throw(ArgumentError("quantities must be a tuple of accessors"))
    selected = isempty(quantities) ? _default_line_requests(object) :
               Tuple(
        request for accessor in quantities
    for request in _line_requests(object, accessor)
    )
    isempty(selected) && throw(ArgumentError("at least one line-parameter accessor is required"))
    length(unique(selected)) == length(selected) || throw(
        ArgumentError("line-parameter accessors select duplicate quantities"),
    )
    return selected
end

function _request_quantity(request)
    request isa Function && return Units.quantity(request)
    return Units.quantity(request...)
end

function _quantity_prefix(quantity_units, request, fallback::Symbol)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    if quantity_units isa NamedTuple
        request isa Function || return fallback
        key = nameof(request)
        return haskey(quantity_units, key) ? getproperty(quantity_units, key) : fallback
    end
    quantity_units isa AbstractDict || throw(
        ArgumentError("quantity_units must be a prefix, selector-keyed collection, or nothing"),
    )
    haskey(quantity_units, request) && return quantity_units[request]
    request isa Function || return fallback
    key = nameof(request)
    return haskey(quantity_units, key) ? quantity_units[key] : fallback
end

function _request_target(request, parameter_basis, length_unit, quantity_units)
    scientific_quantity = _request_quantity(request)
    default = Units.display_unit(
        scientific_quantity,
        parameter_basis;
        length_prefix = length_unit
    )
    isempty(default.numerator) && return default
    fallback = first(default.numerator).prefix
    selected = _quantity_prefix(quantity_units, request, fallback)
    selected isa Units.UnitExpr && return selected
    selected isa Symbol || throw(
        ArgumentError("quantity-unit overrides must be prefixes or UnitExpr values"),
    )
    return Units.display_unit(
        scientific_quantity,
        parameter_basis;
        length_prefix = length_unit,
        prefix = selected
    )
end

function _line_input_defaults(frequencies)
    return (;
        frequencies,
        quantities = (),
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        con = nothing,
        xscale = :linear,
        yscale = :linear
    )
end

function PlotBuilder.dispatch_on(::Type{LineParameterPlotDefinition})
    Union{LineParameters, SeriesImpedance, ShuntAdmittance}
end
function PlotBuilder.input_kwargs(::Type{LineParameterPlotDefinition})
    (
        :frequencies,
        :quantities,
        :freq_unit,
        :length_unit,
        :quantity_units,
        :con,
        :xscale,
        :yscale
    )
end
PlotBuilder.renderer_kwargs(::Type{LineParameterPlotDefinition}) = (:fig_size,)
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

function PlotBuilder.resolve(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    input = recipe.input
    requests = _resolve_line_requests(recipe.object, input.quantities)
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    recipe.object isa Union{SeriesImpedance, ShuntAdmittance} &&
        input.frequencies === nothing && throw(
        ArgumentError("frequencies are required for SeriesImpedance and ShuntAdmittance"),
    )
    supplied_frequencies = input.frequencies === nothing ? nothing : collect(input.frequencies)
    if supplied_frequencies !== nothing
        all(isfinite, supplied_frequencies) || throw(ArgumentError("frequencies must be finite"))
        input.xscale === :log10 && any(<=(0), supplied_frequencies) && throw(
            DomainError(
                supplied_frequencies,
                "logarithmic frequency axes require positive frequencies"
            ),
        )
        any(request -> request === L || request === C, requests) &&
            any(iszero, supplied_frequencies) && throw(DomainError(
            supplied_frequencies,
            "inductance and capacitance are undefined at zero frequency"
        ))
    end
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    all(>(0), recipe.renderer.fig_size) || throw(
        ArgumentError("fig_size dimensions must be positive"),
    )
    return PlotBuilder.PlotRecipe(
        LineParameterPlotDefinition,
        recipe.object,
        merge(input, (; frequencies = supplied_frequencies, requests)),
        recipe.renderer
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

function _published_frequency(object, input)
    target = Units.units(input.freq_unit, :hertz)
    object isa LineParameters || return _frequency_payload(input.frequencies, target)
    published = observables(
        object,
        (frequency = (frequencies, Colon()),);
        units = (frequency = target,)
    ).frequency
    if input.frequencies !== nothing
        supplied = _frequency_payload(input.frequencies, target)
        supplied.values == published.values || throw(
            ArgumentError("supplied frequencies do not match the LineParameters frequency axis"),
        )
    end
    return published
end

function _observable_request(object, request, supplied_frequencies)
    object isa LineParameters && return request
    request === L && return (L, supplied_frequencies)
    request === C && return (C, supplied_frequencies)
    return request
end

function _publish_request(object, request, target, supplied_frequencies)
    selected = _observable_request(object, request, supplied_frequencies)
    return observables(
        object,
        (value = selected,);
        units = (value = target,)
    ).value
end

function _reference_payload(object, request, target, supplied_frequencies)
    return _publish_request(
        object,
        _request_parent(request),
        target,
        supplied_frequencies
    )
end

_is_cartesian_request(request) = request === R || request === X || request === G || request === B

function _maximum_nominal_magnitude(values)
    return mapreduce(
        value -> begin
            magnitude = abs(nominal(value))
            isfinite(magnitude) ? Float64(magnitude) : 0.0
        end,
        max,
        values;
        init = 0.0
    )
end

function _suppress_display_residue(payload, reference, request)
    _is_cartesian_request(request) || return payload
    component_scale = _maximum_nominal_magnitude(payload.values)
    iszero(component_scale) && return payload
    reference_scale = _maximum_nominal_magnitude(reference.values)
    iszero(reference_scale) && return payload
    component_scale <= _DISPLAY_ZERO_RTOL * reference_scale || return payload
    return (;
        values = zero.(payload.values), quantity = payload.quantity, unit = payload.unit)
end

function _publish_line_source(object, input, requests)
    frequency = _published_frequency(object, input)
    observations = map(requests) do request
        target = _request_target(
            request,
            basis(object),
            input.length_unit,
            input.quantity_units
        )
        payload = _publish_request(object, request, target, input.frequencies)
        _is_cartesian_request(request) || return payload
        reference = _reference_payload(object, request, target, input.frequencies)
        return _suppress_display_residue(payload, reference, request)
    end
    sample = first(observations)
    size(sample.values, 3) == length(frequency.values) || throw(
        DimensionMismatch("frequency count does not match line-parameter samples"),
    )
    _conductor_pairs(sample.values, input.con)
    length(frequency.values) <= 1 &&
        @warn "Frequency vector has $(length(frequency.values)) sample(s); nothing to plot."
    return (; frequency, observations)
end

function _supports_log_values(samples)
    found = false
    samples === nothing && return false
    for sample in samples
        found = true
        value = nominal(sample)
        uncertainty = abs(standard_uncertainty(sample))
        value isa Real && isfinite(value) && isfinite(uncertainty) &&
        value - uncertainty > 0 || return false
    end
    return found
end

_axis_scales(values) = _supports_log_values(values) ? (:linear, :log10) : (:linear,)

function _panel_positions(count::Int)
    columns = max(1, ceil(Int, sqrt(count)))
    return Tuple(((index - 1) ÷ columns + 1, (index - 1) % columns + 1)
    for index in 1:count)
end

function _axis_observation(observation, curves)
    values = collect(Iterators.flatten(curve.values for curve in curves))
    return (; values, quantity = observation.quantity, unit = observation.unit)
end

function _line_curves(observation, parent, selector)
    pairs = _conductor_pairs(observation.values, selector)
    family = parent === Z ? "series" : "shunt"
    scientific_symbol = Units.symbol(Units.quantity(parent))
    return map(pairs) do (row, column)
        LineCurve(
            collect(view(observation.values, row, column, :)),
            "$scientific_symbol[$row,$column]",
            Symbol("$(family)_$(row)_$(column)"),
            (; linewidth = 2)
        )
    end
end

function _line_page(recipe, published, parent)
    selected = Tuple(
        (request, observation)
    for (request, observation) in zip(recipe.input.requests, published.observations)
    if _request_parent(request) === parent
    )
    positions = _panel_positions(length(selected))
    panels = map(selected, positions) do selected_request, position
        request, observation = selected_request
        curves = _line_curves(
            observation,
            parent,
            recipe.input.con
        )
        xscales = _axis_scales(published.frequency.values)
        y_observation = _axis_observation(observation, curves)
        yscales = _axis_scales(y_observation.values)
        recipe.input.xscale in xscales || throw(DomainError(
            published.frequency.values,
            "logarithmic frequency axes require positive finite data and uncertainty bounds"
        ))
        recipe.input.yscale in yscales || throw(DomainError(
            y_observation.values,
            "logarithmic ordinate axes require positive finite data and uncertainty bounds"
        ))
        LinePanelPayload(
            request,
            position,
            Units.label(observation.quantity),
            published.frequency,
            y_observation,
            curves,
            recipe.input.xscale,
            recipe.input.yscale,
            xscales,
            yscales,
            (;)
        )
    end
    title = Units.label(Units.quantity(parent))
    return LinePagePayload(
        title,
        (; family = parent, requests = first.(selected)),
        panels,
        PlotBuilder.LegendDefinition(),
        PlotBuilder.ExportDefinition(
            theme = recipe.renderer.export_theme,
            name = title,
            open_file = recipe.renderer.open_export
        )
    )
end

_line_parents(::SeriesImpedance) = (Z,)
_line_parents(::ShuntAdmittance) = (Y,)
_line_parents(::LineParameters) = (Z, Y)

function PlotBuilder.fetch(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    published = _publish_line_source(recipe.object, recipe.input, recipe.input.requests)
    pages = length(published.frequency.values) <= 1 ? () : Tuple(
        _line_page(recipe, published, parent)
    for parent in _line_parents(recipe.object)
    if any(request -> _request_parent(request) === parent, recipe.input.requests)
    )
    return PlotBuilder.PlotRecipe(
        LineParameterPlotDefinition,
        recipe.object,
        merge(recipe.input, (; pages)),
        recipe.renderer
    )
end

PlotBuilder._recipe_variant(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe) =
    Val(:direct)
function PlotBuilder._composition(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        ::PlotBuilder.PlotRecipe
)
    return Val(:empty)
end
function PlotBuilder._page_keys(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        ::Val{:empty},
        recipe::PlotBuilder.PlotRecipe
)
    return eachindex(recipe.input.pages)
end

function _line_page_payload(recipe, page_index::Integer)
    return recipe.input.pages[page_index]
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        recipe::PlotBuilder.PlotRecipe,
        page_index::Integer,
        ::Nothing
)
    return _line_page_payload(recipe, page_index).title
end

function PlotBuilder.default_figsize(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        recipe::PlotBuilder.PlotRecipe,
        ::Integer
)
    return recipe.renderer.fig_size
end

function PlotBuilder.layout_spec(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        ::PlotBuilder.PlotRecipe,
        ::Integer
)
    return :grid
end

function PlotBuilder.page_identity(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        recipe::PlotBuilder.PlotRecipe,
        page_index::Integer
)
    return merge((; page = page_index), _line_page_payload(recipe, page_index).key)
end

function PlotBuilder.legend_spec(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        recipe::PlotBuilder.PlotRecipe,
        page_index::Integer
)
    return _line_page_payload(recipe, page_index).legend
end

function PlotBuilder.export_spec(
        ::Type{LineParameterPlotDefinition},
        ::Val{:direct},
        recipe::PlotBuilder.PlotRecipe,
        page_index::Integer,
        ::AbstractString
)
    return _line_page_payload(recipe, page_index).export_definition
end
