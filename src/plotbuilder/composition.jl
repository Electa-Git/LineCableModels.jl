"""
    resolve_input(::Type{S}, recipe)

Validate and enrich parsed recipe input before materialisation.
"""
resolve_input(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = recipe

"Resolve the semantic observations consumed by subsequent plot stages."
observe(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition} = recipe

"""
    _recipe_variant(::Type{S}, recipe)

Return the value-dispatched plotting mode for a resolved recipe.
"""
function _recipe_variant(::Type{S}, recipe::PlotRecipe) where {S <: AbstractPlotDefinition}
    Val(:default)
end

"""
    _composition(::Type{S}, variant, recipe)

Return the value-dispatched grouping mode for a resolved recipe.
"""
function _composition(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    Val(:overlay)
end

"""
    _page_items(::Type{S}, variant, recipe)

Return semantic page facets for recipes using `:faceted_pages`.
"""
function _page_items(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    (nothing,)
end

"""
    _series_items(::Type{S}, variant, recipe, page_key)

Return semantic series facets for the selected recipe page.
"""
function _series_items(
        ::Type{S},
        mode::Val,
        recipe::PlotRecipe,
        page_key
) where {S <: AbstractPlotDefinition}
    (nothing,)
end

function _page_keys(::Type{S}, mode::Val, ::Val{:overlay},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _page_keys(::Type{S}, mode::Val, ::Val{:panels},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _page_keys(::Type{S}, mode::Val, ::Val{:pages},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    _series_items(S, mode, recipe, nothing)
end
function _page_keys(::Type{S}, mode::Val, ::Val{:faceted_pages},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    _page_items(S, mode, recipe)
end
function _page_keys(::Type{S}, mode::Val, ::Val{:empty},
        recipe::PlotRecipe) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _page_keys(
        ::Type{S},
        mode::Val,
        grouping::Val,
        recipe::PlotRecipe
) where {S <: AbstractPlotDefinition}
    grouping_name = only(typeof(grouping).parameters)
    throw(
        ArgumentError(
        "unsupported grouping mode :$grouping_name for $S; " *
        "specialise PlotBuilder._page_keys for this local definition variant"
    ),
    )
end

function _view_keys(::Type{S}, mode::Val, ::Val{:overlay}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _view_keys(::Type{S}, mode::Val, ::Val{:panels}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    _series_items(S, mode, recipe, page_key)
end
function _view_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _view_keys(::Type{S}, mode::Val, ::Val{:faceted_pages},
        recipe::PlotRecipe, page_key) where {
        S <: AbstractPlotDefinition,
}
    (nothing,)
end
function _view_keys(::Type{S}, mode::Val, ::Val{:empty}, recipe::PlotRecipe,
        page_key) where {
        S <: AbstractPlotDefinition,
}
    ()
end

function _series_keys(::Type{S}, mode::Val, ::Val{:overlay}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    _series_items(S, mode, recipe, page_key)
end
function _series_keys(::Type{S}, mode::Val, ::Val{:panels}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    (view_key,)
end
function _series_keys(::Type{S}, mode::Val, ::Val{:pages}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    (page_key,)
end
function _series_keys(::Type{S}, mode::Val, ::Val{:faceted_pages}, recipe::PlotRecipe,
        page_key, view_key) where {S <: AbstractPlotDefinition}
    _series_items(S, mode, recipe, page_key)
end
