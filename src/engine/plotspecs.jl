struct LineParameterPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

function _component_unit(
        component::Symbol,
        parameter_basis::Symbol,
        length_unit::Symbol,
        quantity_units
)
    resolved = UnitHandler.line_component_unit(
        component,
        parameter_basis;
        length_unit,
        quantity_units
    )
    return resolved.quantity, resolved.units, resolved.scale
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

_line_component_values(::Val{:R}, object, frequencies) = R(object)
_line_component_values(::Val{:X}, object, frequencies) = X(object)
_line_component_values(::Val{:G}, object, frequencies) = G(object)
_line_component_values(::Val{:B}, object, frequencies) = B(object)
_line_component_values(::Val{:Z_re}, object, frequencies) = R(object)
_line_component_values(::Val{:Z_im}, object, frequencies) = X(object)
_line_component_values(::Val{:Z_abs}, object, frequencies) = abs.(Z(object))
function _line_component_values(::Val{:Z_angle}, object, frequencies)
    angle.(Z(object)) .* (180 / π)
end
_line_component_values(::Val{:Y_re}, object, frequencies) = G(object)
_line_component_values(::Val{:Y_im}, object, frequencies) = B(object)
_line_component_values(::Val{:Y_abs}, object, frequencies) = abs.(Y(object))
function _line_component_values(::Val{:Y_angle}, object, frequencies)
    angle.(Y(object)) .* (180 / π)
end

_line_component_values(::Val{:L}, parameters::LineParameters, frequencies) = L(parameters)
_line_component_values(::Val{:C}, parameters::LineParameters, frequencies) = C(parameters)

function _frequency_component_values(component::Symbol, reactive, frequencies)
    any(iszero, frequencies) && throw(
        DomainError(frequencies, "$component is undefined at zero frequency"),
    )
    return reactive ./ reshape(2π .* frequencies, 1, 1, :)
end

function _line_component_values(::Val{:L}, object::SeriesImpedance, frequencies)
    _frequency_component_values(:L, X(object), frequencies)
end
function _line_component_values(::Val{:C}, object::ShuntAdmittance, frequencies)
    _frequency_component_values(:C, B(object), frequencies)
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
    accessor in (R, X, L) && return _line_components(Z(parameters), accessor)
    accessor in (G, B, C) && return _line_components(Y(parameters), accessor)
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

function _line_sources(parameters::LineParameters)
    values = observables(parameters)
    return (values.series_impedance, values.shunt_admittance)
end
_line_sources(object::Union{SeriesImpedance, ShuntAdmittance}) = (object,)

_line_source(object::SeriesImpedance, ::LinePageKey{:series, C}) where {C} = object
_line_source(object::ShuntAdmittance, ::LinePageKey{:shunt, C}) where {C} = object
function _line_source(parameters::LineParameters, ::LinePageKey{:series, C}) where {C}
    observables(parameters).series_impedance
end
function _line_source(parameters::LineParameters, ::LinePageKey{:shunt, C}) where {C}
    observables(parameters).shunt_admittance
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
    _line_input_defaults(frequencies(parameters))
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

function PlotBuilder.resolve_input(::Type{LineParameterPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    input = recipe.input
    components = _resolve_line_components(recipe.object, input.quantities)
    input.xscale in (:linear, :log10) || throw(
        ArgumentError("xscale must be :linear or :log10"),
    )
    input.yscale in (:linear, :log10) || throw(
        ArgumentError("yscale must be :linear or :log10"),
    )
    input.frequencies === nothing && throw(
        ArgumentError("frequencies are required for SeriesImpedance and ShuntAdmittance"),
    )
    frequencies = collect(input.frequencies)
    all(isfinite, frequencies) || throw(ArgumentError("frequencies must be finite"))
    input.xscale === :log10 && any(<=(0), frequencies) &&
        throw(
            DomainError(frequencies, "logarithmic frequency axes require positive frequencies"),
        )
    any(component -> component in (:L, :C), components) && any(iszero, frequencies) &&
        throw(
            DomainError(frequencies, "inductance and capacitance are undefined at zero frequency"),
        )
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    for source in _line_sources(recipe.object)
        size(source, 3) == length(frequencies) || throw(
            DimensionMismatch("frequency count does not match line-parameter samples"),
        )
        _conductor_pairs(source, input.con)
    end
    length(frequencies) <= 1 &&
        @warn "Frequency vector has $(length(frequencies)) sample(s); nothing to plot."
    return PlotBuilder.PlotRecipe(
        LineParameterPlotDefinition,
        recipe.object,
        merge(input, (; frequencies, components)),
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
    length(recipe.input.frequencies) <= 1 && return ()
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
    source = _line_source(recipe.object, view_key)
    return _conductor_pairs(source, recipe.input.con)
end

function _line_values(recipe::PlotBuilder.PlotRecipe, page_key::LinePageKey)
    source = _line_source(recipe.object, page_key)
    component = _line_component(page_key)
    _, _,
    conversion = _component_unit(
        component,
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    values = _line_component_values(
        Val(component),
        recipe.object,
        recipe.input.frequencies
    )
    return _display_values(values, source, component), conversion
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

function _display_values(values, source, component::Symbol)
    component in _CARTESIAN_COMPONENTS || return values
    component_scale = _maximum_nominal_magnitude(values)
    iszero(component_scale) && return values
    reference_scale = _maximum_nominal_magnitude(source)
    iszero(reference_scale) && return values
    component_scale <= _DISPLAY_ZERO_RTOL * reference_scale || return values
    return zero.(values)
end

function PlotBuilder.axis_quantity(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return UnitHandler.QuantityTag{:frequency}()
end

function PlotBuilder.axis_quantity(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    source = _line_source(recipe.object, page_key)
    quantity, _,
    _ = _component_unit(
        _line_component(page_key),
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    return quantity
end

function PlotBuilder.axis_quantity(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    return PlotBuilder.axis_quantity(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        view_key,
        nothing
    )
end

function PlotBuilder.axis_unit(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        quantity::UnitHandler.QuantityTag,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return UnitHandler.units(recipe.input.freq_unit, :hertz)
end

function PlotBuilder.axis_unit(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    source = _line_source(recipe.object, page_key)
    _, target,
    _ = _component_unit(
        _line_component(page_key),
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    return target
end

function PlotBuilder.axis_unit(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey
)
    return PlotBuilder.axis_unit(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        quantity,
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

function PlotBuilder.series_data(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series_key
)
    quantity = UnitHandler.QuantityTag{:frequency}()
    target = UnitHandler.units(recipe.input.freq_unit, :hertz)
    conversion = UnitHandler.scale_factor(UnitHandler.default_unit(quantity), target)
    return recipe.input.frequencies .* conversion
end

function PlotBuilder.series_data(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key,
        series_key::Tuple{Int, Int}
)
    values, conversion = _line_values(recipe, page_key)
    return collect(view(values, series_key[1], series_key[2], :)) .* conversion
end

function PlotBuilder.series_data(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val,
        view_key::LinePageKey,
        series_key::Tuple{Int, Int}
)
    return PlotBuilder.series_data(
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
    quantity = PlotBuilder.axis_quantity(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    )
    return "$(UnitHandler.get_symbol(quantity))[$(series_key[1]),$(series_key[2])]"
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
    quantity = PlotBuilder.axis_quantity(
        LineParameterPlotDefinition,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    )
    return UnitHandler.get_label(quantity)
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{:series},
        ::Nothing
)
    return "Series impedance"
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::Val{:shunt},
        ::Nothing
)
    return "Shunt admittance"
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
                   _supports_log_values(recipe.input.frequencies) :
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
