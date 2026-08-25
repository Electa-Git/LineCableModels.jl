struct LineParameterPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

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

function _finite_exponent(curves)
    maximum_value = 0.0
    for curve in curves, sample in curve

        value = abs(nominal(sample))
        isfinite(value) && (maximum_value = max(maximum_value, value))
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

struct LinePageKey{K, C} end

_line_parent(::LinePageKey{K, C}) where {K, C} = K
_line_component(::LinePageKey{K, C}) where {K, C} = C
_line_parent(::Val{K}) where {K} = K

const _SERIES_COMPONENTS = (:R, :X, :L, :Z_re, :Z_im, :Z_abs, :Z_angle)
const _SHUNT_COMPONENTS = (:G, :B, :C, :Y_re, :Y_im, :Y_abs, :Y_angle)
const _CARTESIAN_COMPONENTS = (:R, :X, :G, :B, :Z_re, :Z_im, :Y_re, :Y_im)
const _DISPLAY_ZERO_RTOL = sqrt(eps(Float64))

_component_request(::Val{:R}, frequencies) = R
_component_request(::Val{:X}, frequencies) = X
_component_request(::Val{:L}, frequencies) = frequencies === nothing ? L : (L, frequencies)
_component_request(::Val{:G}, frequencies) = G
_component_request(::Val{:B}, frequencies) = B
_component_request(::Val{:C}, frequencies) = frequencies === nothing ? C : (C, frequencies)
_component_request(::Val{:Z_re}, frequencies) = R
_component_request(::Val{:Z_im}, frequencies) = X
_component_request(::Val{:Z_abs}, frequencies) = (Z, abs)
_component_request(::Val{:Z_angle}, frequencies) = (Z, angle)
_component_request(::Val{:Y_re}, frequencies) = G
_component_request(::Val{:Y_im}, frequencies) = B
_component_request(::Val{:Y_abs}, frequencies) = (Y, abs)
_component_request(::Val{:Y_angle}, frequencies) = (Y, angle)

function _request_quantity(request)
    request isa Function && return Units.quantity(request)
    selector = first(request)
    length(request) >= 2 && applicable(Units.quantity, selector, request[2]) &&
        return Units.quantity(selector, request[2])
    return Units.quantity(selector)
end

function _quantity_prefix(quantity_units, component::Symbol, fallback::Symbol)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    quantity_units isa NamedTuple || quantity_units isa AbstractDict ||
        throw(
            ArgumentError("quantity_units must be a prefix, keyed collection, or nothing"),
        )
    return haskey(quantity_units, component) ? quantity_units[component] : fallback
end

function _component_target(component, request, parameter_basis, length_unit, quantity_units)
    scientific_quantity = _request_quantity(request)
    default = Units.display_unit(scientific_quantity, parameter_basis; length_prefix = length_unit)
    isempty(default.numerator) && return default
    fallback = first(default.numerator).prefix
    selected = _quantity_prefix(quantity_units, component, fallback)
    selected isa Units.UnitExpr && return selected
    selected isa Symbol ||
        throw(ArgumentError("quantity-unit overrides must be prefixes or UnitExpr values"))
    return Units.display_unit(
        scientific_quantity,
        parameter_basis;
        length_prefix = length_unit,
        prefix = selected
    )
end

function _line_page_keys(::Val{K}, ::Val{Components}) where {K, Components}
    allowed = K === :series ? _SERIES_COMPONENTS : _SHUNT_COMPONENTS
    return Tuple(LinePageKey{K, component}()
    for component in Components if component in allowed)
end

_default_line_components(::LineParameters) = (:Z_re, :Z_im, :Y_re, :Y_im)
_default_line_components(::SeriesImpedance) = (:Z_re, :Z_im)
_default_line_components(::ShuntAdmittance) = (:Y_re, :Y_im)

function _line_components(::SeriesImpedance, accessor)
    accessor === R && return (:R,)
    accessor === X && return (:X,)
    accessor === L && return (:L,)
    accessor === real && return (:Z_re,)
    accessor === imag && return (:Z_im,)
    accessor === abs && return (:Z_abs,)
    accessor === angle && return (:Z_angle,)
    accessor === Z && return (:Z_re, :Z_im)
    throw(ArgumentError("accessor $(accessor) is not defined for SeriesImpedance presentation"))
end

function _line_components(::ShuntAdmittance, accessor)
    accessor === G && return (:G,)
    accessor === B && return (:B,)
    accessor === C && return (:C,)
    accessor === real && return (:Y_re,)
    accessor === imag && return (:Y_im,)
    accessor === abs && return (:Y_abs,)
    accessor === angle && return (:Y_angle,)
    accessor === Y && return (:Y_re, :Y_im)
    throw(ArgumentError("accessor $(accessor) is not defined for ShuntAdmittance presentation"))
end

function _line_components(parameters::LineParameters, accessor)
    accessor === Z && return (:Z_re, :Z_im)
    accessor === Y && return (:Y_re, :Y_im)
    accessor === real && return (:Z_re, :Y_re)
    accessor === imag && return (:Z_im, :Y_im)
    accessor === abs && return (:Z_abs, :Y_abs)
    accessor === angle && return (:Z_angle, :Y_angle)
    accessor === R && return (:R,)
    accessor === X && return (:X,)
    accessor === L && return (:L,)
    accessor === G && return (:G,)
    accessor === B && return (:B,)
    accessor === C && return (:C,)
    throw(ArgumentError("accessor $(accessor) is not defined for LineParameters presentation"))
end

function _resolve_line_components(object, quantities)
    quantities isa Tuple || throw(ArgumentError("quantities must be a tuple of accessors"))
    selected = isempty(quantities) ? _default_line_components(object) :
               Tuple(
        component for accessor in quantities
    for component in _line_components(object, accessor)
    )
    isempty(selected) &&
        throw(ArgumentError("at least one line-parameter accessor is required"))
    length(unique(selected)) == length(selected) || throw(
        ArgumentError("line-parameter accessors select duplicate quantities"),
    )
    return selected
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
function PlotBuilder.input_defaults(::Type{LineParameterPlotDefinition}, parameters::LineParameters)
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
    components = _resolve_line_components(recipe.object, input.quantities)
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    recipe.object isa Union{SeriesImpedance, ShuntAdmittance} &&
        input.frequencies === nothing &&
        throw(
            ArgumentError("frequencies are required for SeriesImpedance and ShuntAdmittance"),
        )
    frequencies = input.frequencies === nothing ? nothing : collect(input.frequencies)
    if frequencies !== nothing
        all(isfinite, frequencies) || throw(ArgumentError("frequencies must be finite"))
        input.xscale === :log10 && any(<=(0), frequencies) &&
            throw(
                DomainError(frequencies, "logarithmic frequency axes require positive frequencies"),
            )
        any(component -> component in (:L, :C), components) && any(iszero, frequencies) &&
            throw(DomainError(
                frequencies,
                "inductance and capacitance are undefined at zero frequency"
            ))
    end
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    return PlotBuilder.PlotRecipe(
        LineParameterPlotDefinition,
        recipe.object,
        merge(input, (; frequencies, components)),
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
    object isa LineParameters || return _frequency_payload(
        input.frequencies,
        target
    )
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

function _reference_payload(object, component, target)
    parent = component in _SERIES_COMPONENTS ? Z : Y
    return observables(
        object,
        (reference = parent,);
        units = (reference = target,)
    ).reference
end

function _suppress_display_residue(payload, reference, component)
    component in _CARTESIAN_COMPONENTS || return payload
    component_scale = _maximum_nominal_magnitude(payload.values)
    iszero(component_scale) && return payload
    reference_scale = _maximum_nominal_magnitude(reference.values)
    iszero(reference_scale) && return payload
    component_scale <= _DISPLAY_ZERO_RTOL * reference_scale || return payload
    return (;
        values = zero.(payload.values), quantity = payload.quantity, unit = payload.unit)
end

function _publish_line_source(object, input, components)
    frequency = _published_frequency(object, input)
    native_frequencies = object isa LineParameters ? nothing : input.frequencies
    requests = NamedTuple{components}(map(
        component -> _component_request(Val(component), native_frequencies),
        components
    ))
    targets = NamedTuple{components}(map(components) do component
        request = _component_request(Val(component), native_frequencies)
        _component_target(
            component,
            request,
            basis(object),
            input.length_unit,
            input.quantity_units
        )
    end)
    published = observables(object, requests; units = targets)
    payloads = map(
        components,
        Base.values(published),
        Base.values(targets)
    ) do component, payload, target
        component in _CARTESIAN_COMPONENTS || return payload
        reference = _reference_payload(object, component, target)
        _suppress_display_residue(payload, reference, component)
    end
    component_payloads = NamedTuple{components}(payloads)
    sample = first(Base.values(component_payloads))
    size(sample.values, 3) == length(frequency.values) || throw(
        DimensionMismatch("frequency count does not match line-parameter samples"),
    )
    _conductor_pairs(sample.values, input.con)
    length(frequency.values) <= 1 &&
        @warn "Frequency vector has $(length(frequency.values)) sample(s); nothing to plot."
    return (; frequency, component_payloads)
end

function PlotBuilder.fetch(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    published = _publish_line_source(recipe.object, recipe.input, recipe.input.components)
    input = merge(recipe.input, published)
    return PlotBuilder.PlotRecipe(
        LineParameterPlotDefinition,
        recipe.object,
        input,
        recipe.renderer
    )
end

function PlotBuilder._recipe_variant(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    return Val(recipe.input.components)
end

function PlotBuilder._composition(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:panels)
end

_line_family_facets(::SeriesImpedance) = (Val(:series),)
_line_family_facets(::ShuntAdmittance) = (Val(:shunt),)
function _line_family_facets(::LineParameters)
    return (Val(:series), Val(:shunt))
end

function _line_component_facets(
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        ::Val{K}
) where {K}
    return _line_page_keys(Val(K), mode)
end

function PlotBuilder._page_keys(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe
)
    length(recipe.input.frequency.values) <= 1 && return ()
    return Tuple(
        family
    for family in _line_family_facets(recipe.object)
    if !isempty(_line_component_facets(mode, recipe, family))
    )
end

function PlotBuilder._view_keys(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val
)
    return _line_component_facets(mode, recipe, page_key)
end

function PlotBuilder._series_keys(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:panels},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    payload = getproperty(recipe.input.component_payloads, _line_component(view_key))
    return _conductor_pairs(payload.values, recipe.input.con)
end

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

function PlotBuilder.axis_payload(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.frequency
end

function PlotBuilder.axis_payload(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    return getproperty(recipe.input.component_payloads, _line_component(page_key))
end

function PlotBuilder.axis_payload(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    return PlotBuilder.axis_payload(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        view_key,
        nothing
    )
end

function PlotBuilder.axis_scale(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.xscale
end

function PlotBuilder.axis_scale(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.yscale
end

function PlotBuilder.series_values(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series_key
)
    return recipe.input.frequency.values
end

function PlotBuilder.series_values(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key,
        series_key::Tuple{Int, Int}
)
    payload = getproperty(recipe.input.component_payloads, _line_component(page_key))
    return collect(view(payload.values, series_key[1], series_key[2], :))
end

function PlotBuilder.series_values(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey,
        series_key::Tuple{Int, Int}
)
    return PlotBuilder.series_values(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        view_key,
        nothing,
        series_key
    )
end

function PlotBuilder.legend_label(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key,
        series_key::Tuple{Int, Int}
)
    quantity = PlotBuilder.axis_payload(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    ).quantity
    return "$(Units.symbol(quantity))[$(series_key[1]),$(series_key[2])]"
end

function PlotBuilder.legend_label(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{K},
        view_key::LinePageKey,
        series_key::Tuple{Int, Int}
) where {K}
    symbol = K === :series ? "Z" : "Y"
    return "$symbol[$(series_key[1]),$(series_key[2])]"
end

function PlotBuilder.series_group(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{K},
        view_key::LinePageKey,
        series_key::Tuple{Int, Int}
) where {K}
    return Symbol("$(K)_$(series_key[1])_$(series_key[2])")
end

function PlotBuilder.series_attributes(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series_key
)
    return (; linewidth = 2)
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    quantity = PlotBuilder.axis_payload(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    ).quantity
    return Units.label(quantity)
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{:series},
        ::Nothing
)
    return Units.label(Units.quantity(Z))
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{:shunt},
        ::Nothing
)
    return Units.label(Units.quantity(Y))
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    return PlotBuilder.default_title(
        LineParameterPlotDefinition,
        mode,
        recipe,
        view_key,
        nothing
    )
end

function PlotBuilder.view_key(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    return (; component = _line_component(page_key))
end

function PlotBuilder.view_key(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    return (; component = _line_component(view_key))
end

function PlotBuilder.layout_spec(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val
)
    return :grid
end

function PlotBuilder.default_figsize(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.renderer.fig_size
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

function _supports_log(series, dim::Symbol)
    found = false
    for item in series
        samples = dim === :x ? item.xdata : item.ydata
        samples === nothing && continue
        for sample in samples
            found = true
            value = nominal(sample)
            uncertainty = abs(standard_uncertainty(sample))
            value isa Real && isfinite(value) && isfinite(uncertainty) &&
            value - uncertainty > 0 || return false
        end
    end
    return found
end

function PlotBuilder.axis_scales(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{dim},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series::Vector{PlotBuilder.SeriesSpec}
) where {dim}
    supports_log = dim === :x ?
                   _supports_log_values(recipe.input.frequency.values) :
                   _supports_log(series, dim)
    return supports_log ? (:linear, :log10) : (:linear,)
end

function PlotBuilder.axis_exponent(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{dim},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key,
        series::Vector{PlotBuilder.SeriesSpec}
) where {dim}
    return _finite_exponent(
        (dim === :x ? item.xdata : item.ydata for item in series)
    )
end

function PlotBuilder.axis_exponent(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        dim::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey,
        series::Vector{PlotBuilder.SeriesSpec}
)
    return PlotBuilder.axis_exponent(
        LineParameterPlotDefinition,
        mode,
        dim,
        recipe,
        view_key,
        nothing,
        series
    )
end

function PlotBuilder.page_identity(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return (;
        family = _line_parent(page_key),
        component = _line_component(page_key),
        components = recipe.input.components,
        conductors = recipe.input.con
    )
end

function PlotBuilder.page_identity(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val
)
    return (;
        family = _line_parent(page_key),
        components = Tuple(
            _line_component(component)
        for component in _line_component_facets(mode, recipe, page_key)
        ),
        conductors = recipe.input.con
    )
end
