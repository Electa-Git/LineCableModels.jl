const COMMON_RENDERER_KWARGS = (:export_theme, :open_export, :layout)
const EXPORT_THEMES = (:default, :publication)

function _validate_export_theme(value::Symbol)
    value in EXPORT_THEMES || throw(
        ArgumentError("export_theme must be :default or :publication"),
    )
    return value
end

"""
$(TYPEDSIGNATURES)

Export the current typed state of a `UIPlot` through an explicitly loaded
CairoMakie extension.
"""
function export_svg end

"Return the domain object accepted by a plot definition."
function entitle(::Type{S}, object) where {S <: AbstractPlotDefinition}
    expected = dispatch_on(S)
    object isa expected || throw(
        ArgumentError("$S accepts $expected, not $(typeof(object))"),
    )
    return object
end

"""
    dispatch_on(::Type{S})

Return the domain type accepted by plot definition `S`.
"""
dispatch_on(::Type{S}) where {S <: AbstractPlotDefinition} = Any

"""
    input_kwargs(::Type{S})

Return semantic keyword names accepted by plot definition `S`.
"""
input_kwargs(::Type{S}) where {S <: AbstractPlotDefinition} = ()

"""
    renderer_kwargs(::Type{S})

Return recipe-specific renderer keyword names accepted by plot definition `S`.
"""
renderer_kwargs(::Type{S}) where {S <: AbstractPlotDefinition} = ()

"""
    input_defaults(::Type{S}, object)

Return defaults for every name declared by `input_kwargs(S)`.
"""
input_defaults(::Type{S}, object) where {S <: AbstractPlotDefinition} = (;)

"""
    renderer_defaults(::Type{S}, object)

Return defaults for every name declared by `renderer_kwargs(S)`.
"""
renderer_defaults(::Type{S}, object) where {S <: AbstractPlotDefinition} = (;)

function _symbol_tuple(::Type{S}, values, accessor::Symbol) where {S <:
                                                                   AbstractPlotDefinition}
    values isa Tuple || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    all(value -> value isa Symbol, values) || throw(
        ArgumentError("$accessor($S) must return a tuple of symbols"),
    )
    length(unique(values)) == length(values) || throw(
        ArgumentError("$accessor($S) cannot declare duplicate keywords"),
    )
    return values
end

function _select_kwargs(kwargs::NamedTuple, names::Tuple)
    selected = Tuple(pair for pair in pairs(kwargs) if first(pair) in names)
    return NamedTuple(selected)
end

function _validate_defaults(::Type{S}, defaults::NamedTuple, names::Tuple,
        accessor::Symbol) where {
        S <: AbstractPlotDefinition,
}
    actual = Tuple(keys(defaults))
    Set(actual) == Set(names) || throw(
        ArgumentError(
        "$accessor($S) must define exactly the declared keywords; " *
        "declared $(collect(names)), received $(collect(actual))"
    ),
    )
    return defaults
end

"""
    parse(::Type{S}, object, kwargs)

Validate and split caller input into the semantic and renderer options declared
by a recipe. Unsupported keywords are errors.
"""
function parse(
        ::Type{S},
        object,
        kwargs::NamedTuple
) where {S <: AbstractPlotDefinition}
    input_names = _symbol_tuple(S, input_kwargs(S), :input_kwargs)
    renderer_names = _symbol_tuple(S, renderer_kwargs(S), :renderer_kwargs)
    collisions = intersect(input_names, renderer_names)
    isempty(collisions) || throw(
        ArgumentError("$S declares keywords in both input and renderer options: $(join(collisions, ", "))"),
    )
    common_collisions = intersect((input_names..., renderer_names...), COMMON_RENDERER_KWARGS)
    isempty(common_collisions) || throw(
        ArgumentError("$S redeclares common renderer keywords: $(join(common_collisions, ", "))"),
    )
    allowed_renderer = (COMMON_RENDERER_KWARGS..., renderer_names...)
    allowed = (input_names..., allowed_renderer...)
    unsupported = Tuple(name for name in keys(kwargs) if name ∉ allowed)
    isempty(unsupported) || throw(
        ArgumentError("unsupported plot keyword(s) for $S: $(join(unsupported, ", "))"),
    )

    declared_input = input_defaults(S, object)
    declared_input isa NamedTuple || throw(
        ArgumentError("input_defaults($S) must return a NamedTuple"),
    )
    declared_renderer = renderer_defaults(S, object)
    declared_renderer isa NamedTuple || throw(
        ArgumentError("renderer_defaults($S) must return a NamedTuple"),
    )
    _validate_defaults(S, declared_input, input_names, :input_defaults)
    _validate_defaults(S, declared_renderer, renderer_names, :renderer_defaults)

    input = merge(declared_input, _select_kwargs(kwargs, input_names))
    renderer = merge(
        (; export_theme = :default, open_export = true, layout = nothing),
        declared_renderer,
        _select_kwargs(kwargs, allowed_renderer)
    )
    _validate_export_theme(renderer.export_theme)
    renderer.open_export isa Bool || throw(ArgumentError("open_export must be Bool"))
    renderer.layout isa Union{Nothing, Symbol, LayoutSpec} || throw(
        ArgumentError("layout must be nothing, a preset symbol, or LayoutSpec"),
    )
    return PlotRecipe(S, object, input, renderer)
end

function parse(::Type{S}, object; kwargs...) where {S <: AbstractPlotDefinition}
    parse(S, object, (; kwargs...))
end
