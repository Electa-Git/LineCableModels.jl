const COMMON_RENDERER_KWARGS = (:export_theme, :open_export)
const EXPORT_THEMES = (:default, :publication)

"""
$(TYPEDSIGNATURES)

Validate an SVG export theme and return it unchanged.
"""
function validate_export_theme(value::Symbol)
    value in EXPORT_THEMES || throw(
        ArgumentError("export_theme must be :default or :publication"),
    )
    return value
end

"""
$(TYPEDSIGNATURES)

Export the current typed state of a [`UIPlot`](@ref) through an explicitly
loaded CairoMakie extension.
"""
function export_svg end

"""
$(SIGNATURES)

Build one standard plotting window and let `callback` draw directly into its
native Makie canvas. A Makie extension supplies the implementation.
"""
function plotwindow end

"""
$(SIGNATURES)

Create a Makie axis on a plotting window, derive scientific labels from the
published observations, and register the axis with the standard controls.
"""
function axis! end

"""
$(SIGNATURES)

Register an ordinary Makie axis with the standard formatter, limits, legend,
interaction, and export machinery.
"""
function register! end

"Return defaults for every semantic keyword accepted by a plot definition."
input_defaults(::Type{D}, source) where {D <: AbstractPlotDefinition} = (;)

"Return defaults for every definition-specific renderer keyword accepted."
renderer_defaults(::Type{D}, source) where {D <: AbstractPlotDefinition} = (;)

"Accept one definition-owned source before any later stage runs."
function entitle end

function _select_kwargs(kwargs::NamedTuple, names::Tuple)
    selected = Tuple(pair for pair in pairs(kwargs) if first(pair) in names)
    return NamedTuple(selected)
end

"""
$(SIGNATURES)

Normalize caller keywords into separate semantic and renderer named tuples.
Unsupported keywords are errors.
"""
function parse(
        ::Type{D},
        source,
        kwargs::NamedTuple
) where {D <: AbstractPlotDefinition}
    declared_input = input_defaults(D, source)
    declared_input isa NamedTuple || throw(
        ArgumentError("input_defaults($D) must return a NamedTuple"),
    )
    declared_renderer = renderer_defaults(D, source)
    declared_renderer isa NamedTuple || throw(
        ArgumentError("renderer_defaults($D) must return a NamedTuple"),
    )
    input_names = Tuple(keys(declared_input))
    renderer_names = Tuple(keys(declared_renderer))
    collisions = intersect(input_names, renderer_names)
    isempty(collisions) || throw(
        ArgumentError("$D declares keywords in both input and renderer options: $(join(collisions, ", "))"),
    )
    common_collisions = intersect((input_names..., renderer_names...), COMMON_RENDERER_KWARGS)
    isempty(common_collisions) || throw(
        ArgumentError("$D redeclares common renderer keywords: $(join(common_collisions, ", "))"),
    )
    allowed_renderer = (COMMON_RENDERER_KWARGS..., renderer_names...)
    allowed = (input_names..., allowed_renderer...)
    unsupported = Tuple(name for name in keys(kwargs) if name ∉ allowed)
    isempty(unsupported) || throw(
        ArgumentError("unsupported plot keyword(s) for $D: $(join(unsupported, ", "))"),
    )

    input = merge(declared_input, _select_kwargs(kwargs, input_names))
    renderer = merge(
        (; export_theme = :default, open_export = true),
        declared_renderer,
        _select_kwargs(kwargs, allowed_renderer)
    )
    validate_export_theme(renderer.export_theme)
    renderer.open_export isa Bool || throw(ArgumentError("open_export must be Bool"))
    return (; input, renderer)
end

function parse(::Type{D}, source; kwargs...) where {D <: AbstractPlotDefinition}
    return parse(D, source, (; kwargs...))
end

"Validate semantic choices and resolve definition-owned selections."
function resolve end

"Construct completed detached pages from one resolved plot request."
function fetch end

"Validate page cardinality and construct one final detached plot recipe."
function finish(
        ::Type{D},
        pages::P
) where {D <: AbstractPlotDefinition, P <: Union{Tuple, AbstractVector}}
    return PlotRecipe(D, pages)
end
