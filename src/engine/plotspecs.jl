struct LineParameterPlotSpec <: PlotBuilder.AbstractPlotSpec end

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

_line_parameter_kind(::SeriesImpedance) = :series
_line_parameter_kind(::ShuntAdmittance) = :shunt

function _finite_exponent(curves)
    maximum_value = 0.0
    for curve in curves, sample in curve

        nominal = abs(Measurements.value(sample))
        isfinite(nominal) && (maximum_value = max(maximum_value, nominal))
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

struct LinePageKey{K, C} end

_line_parent(::LinePageKey{K, C}) where {K, C} = K
_line_component(::LinePageKey{K, C}) where {K, C} = C

_line_components(::Val{:series}, ::Val{(:RLCG, :cart)}) = (:R, :L)
_line_components(::Val{:series}, ::Val{(:RLCG, :polar)}) = (:R, :L)
_line_components(::Val{:shunt}, ::Val{(:RLCG, :cart)}) = (:G, :C)
_line_components(::Val{:shunt}, ::Val{(:RLCG, :polar)}) = (:G, :C)
_line_components(::Val{:series}, ::Val{(:ZY, :cart)}) = (:Z_re, :Z_im)
_line_components(::Val{:series}, ::Val{(:ZY, :polar)}) = (:Z_abs, :Z_angle)
_line_components(::Val{:shunt}, ::Val{(:ZY, :cart)}) = (:Y_re, :Y_im)
_line_components(::Val{:shunt}, ::Val{(:ZY, :polar)}) = (:Y_abs, :Y_angle)

function _line_page_keys(::Val{K}, mode::Val) where {K}
    return Tuple(LinePageKey{K, component}()
    for component in _line_components(Val(K), mode))
end

_line_sources(parameters::LineParameters) = (parameters.Z, parameters.Y)
_line_sources(object::Union{SeriesImpedance, ShuntAdmittance}) = (object,)

_line_source(object::SeriesImpedance, ::LinePageKey{:series, C}) where {C} = object
_line_source(object::ShuntAdmittance, ::LinePageKey{:shunt, C}) where {C} = object
_line_source(parameters::LineParameters, ::LinePageKey{:series, C}) where {C} = parameters.Z
_line_source(parameters::LineParameters, ::LinePageKey{:shunt, C}) where {C} = parameters.Y

function _line_input_defaults(frequencies)
    return (;
        frequencies,
        mode = :ZY,
        coord = :cart,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        con = nothing,
        xscale = :linear,
        yscale = :linear
    )
end

function PlotBuilder.dispatch_on(::Type{LineParameterPlotSpec})
    Union{LineParameters, SeriesImpedance, ShuntAdmittance}
end
function PlotBuilder.input_kwargs(::Type{LineParameterPlotSpec})
    (
        :frequencies,
        :mode,
        :coord,
        :freq_unit,
        :length_unit,
        :quantity_units,
        :con,
        :xscale,
        :yscale
    )
end
PlotBuilder.renderer_kwargs(::Type{LineParameterPlotSpec}) = (:fig_size,)
function PlotBuilder.input_defaults(::Type{LineParameterPlotSpec}, parameters::LineParameters)
    _line_input_defaults(parameters.f)
end
function PlotBuilder.input_defaults(
        ::Type{LineParameterPlotSpec},
        ::Union{SeriesImpedance, ShuntAdmittance}
)
    _line_input_defaults(nothing)
end
function PlotBuilder.renderer_defaults(
        ::Type{LineParameterPlotSpec},
        ::Union{LineParameters, SeriesImpedance, ShuntAdmittance}
)
    (; fig_size = (800, 400))
end

function PlotBuilder.resolve_input(::Type{LineParameterPlotSpec}, recipe::PlotBuilder.PlotRecipe)
    input = recipe.input
    input.mode in (:RLCG, :ZY) || throw(ArgumentError("mode must be :RLCG or :ZY"))
    input.coord in (:cart, :polar) || throw(ArgumentError("coord must be :cart or :polar"))
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
    input.mode === :RLCG && any(iszero, frequencies) &&
        throw(
            DomainError(frequencies, "RLCG plotting is undefined at zero frequency"),
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
        recipe.object,
        merge(input, (; frequencies)),
        recipe.renderer
    )
end

function PlotBuilder.recipe_mode(::Type{LineParameterPlotSpec}, recipe::PlotBuilder.PlotRecipe)
    return Val((recipe.input.mode, recipe.input.coord))
end

function PlotBuilder.grouping_mode(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:faceted_pages)
end

function PlotBuilder.page_facets(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    length(recipe.input.frequencies) <= 1 && return ()
    return _line_page_facets(mode, recipe.object)
end

_line_page_facets(mode::Val, ::SeriesImpedance) = _line_page_keys(Val(:series), mode)
_line_page_facets(mode::Val, ::ShuntAdmittance) = _line_page_keys(Val(:shunt), mode)
function _line_page_facets(mode::Val, ::LineParameters)
    return (_line_page_keys(Val(:series), mode)..., _line_page_keys(Val(:shunt), mode)...)
end

function _line_values(recipe::PlotBuilder.PlotRecipe, page_key::LinePageKey)
    source = _line_source(recipe.object, page_key)
    component = _line_component(page_key)
    _, _, conversion = _component_unit(
        component,
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    values = UnitHandler.line_component_values(
        component,
        source.values,
        recipe.input.frequencies
    )
    return values, conversion
end

function PlotBuilder.group_facets(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    source = _line_source(recipe.object, page_key)
    pairs = _conductor_pairs(source, recipe.input.con)
    values, conversion = _line_values(recipe, page_key)
    return [pair
            for pair in pairs
            if any(
        sample -> abs(Measurements.value(sample)) > eps(Float64),
        view(values, pair[1], pair[2], :) .* conversion
    )]
end

function PlotBuilder.axis_quantity(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return UnitHandler.QuantityTag{:freq}()
end

function PlotBuilder.axis_quantity(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    source = _line_source(recipe.object, page_key)
    quantity, _, _ = _component_unit(
        _line_component(page_key),
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    return quantity
end

function PlotBuilder.axis_unit(
        ::Type{LineParameterPlotSpec},
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
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    source = _line_source(recipe.object, page_key)
    _, target, _ = _component_unit(
        _line_component(page_key),
        basis(source),
        recipe.input.length_unit,
        recipe.input.quantity_units
    )
    return target
end

function PlotBuilder.axis_scale(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.xscale
end

function PlotBuilder.axis_scale(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:y},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.yscale
end

function PlotBuilder.series_data(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{:x},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series_key
)
    quantity = UnitHandler.QuantityTag{:freq}()
    target = UnitHandler.units(recipe.input.freq_unit, :hertz)
    conversion = UnitHandler.scale_factor(UnitHandler.default_unit(quantity), target)
    return recipe.input.frequencies .* conversion
end

function PlotBuilder.series_data(
        ::Type{LineParameterPlotSpec},
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

function PlotBuilder.legend_label(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key,
        series_key::Tuple{Int, Int}
)
    quantity = PlotBuilder.axis_quantity(
        LineParameterPlotSpec,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    )
    return "$(UnitHandler.get_symbol(quantity))[$(series_key[1]),$(series_key[2])]"
end

function PlotBuilder.series_attributes(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series_key
)
    return (; linewidth = 2)
end

function PlotBuilder.default_title(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    quantity = PlotBuilder.axis_quantity(
        LineParameterPlotSpec,
        mode,
        Val(:y),
        recipe,
        page_key,
        view_key
    )
    return UnitHandler.get_label(quantity)
end

function PlotBuilder.view_key(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey,
        view_key
)
    return (; component = _line_component(page_key))
end

function PlotBuilder.default_figsize(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.renderer.fig_size
end

function _supports_log(series, dim::Symbol)
    found = false
    for item in series
        samples = dim === :x ? item.xdata : item.ydata
        samples === nothing && continue
        for sample in samples
            found = true
            nominal = Measurements.value(sample)
            uncertainty = abs(Measurements.uncertainty(sample))
            nominal isa Real && isfinite(nominal) && isfinite(uncertainty) &&
            nominal - uncertainty > 0 || return false
        end
    end
    return found
end

function PlotBuilder.axis_scales(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        ::Val{dim},
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key,
        series::Vector{PlotBuilder.SeriesSpec}
) where {dim}
    return _supports_log(series, dim) ? (:linear, :log10) : (:linear,)
end

function PlotBuilder.axis_exponent(
        ::Type{LineParameterPlotSpec},
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

function PlotBuilder.page_identity(
        ::Type{LineParameterPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key::LinePageKey
)
    return (;
        family = _line_parent(page_key),
        component = _line_component(page_key),
        mode = recipe.input.mode,
        coordinates = recipe.input.coord,
        conductors = recipe.input.con
    )
end
