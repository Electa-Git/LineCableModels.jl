"Apply page-level decorations after layout resolution."
function decorate(
        ::Type{S}, recipe::PlotRecipe, pages::Vector{PageSpec}
) where {S <: AbstractPlotDefinition}
    return pages
end

"Finish a completed renderer-independent plot recipe."
function finish(
        ::Type{S}, recipe::PlotRecipe, pages::Vector{PageSpec}
) where {S <: AbstractPlotDefinition}
    return PlotRecipe(S, recipe.object, recipe.input, recipe.renderer, pages)
end

"""
    make_render(::Type{S}, object; kwargs...)

Run the fixed PlotBuilder sequence for a domain object. Plot definitions
specialise its methods without replacing the sequence.
"""
function make_render(::Type{S}, object; kwargs...) where {S <: AbstractPlotDefinition}
    entitled = entitle(S, object)
    recipe = parse(S, entitled; kwargs...)
    recipe = resolve(S, recipe)
    recipe isa PlotRecipe || throw(
        ArgumentError("resolve($S) must return PlotRecipe"),
    )
    recipe = fetch(S, recipe)
    recipe isa PlotRecipe || throw(
        ArgumentError("fetch($S) must return PlotRecipe"),
    )
    mode = _recipe_variant(S, recipe)
    mode isa Val || throw(ArgumentError("plot definition variant must be a Val"))
    grouping = _composition(S, mode, recipe)
    grouping isa Val || throw(
        ArgumentError("plot definition composition must be a Val"),
    )
    axes = make_axes(S, mode, grouping, recipe)
    series = make_series(S, mode, grouping, recipe, axes)
    views = make_views(S, mode, grouping, recipe, series)
    pages = make_pages(S, mode, grouping, recipe, views)
    pages = decorate(S, recipe, pages)
    return finish(S, recipe, pages)
end
