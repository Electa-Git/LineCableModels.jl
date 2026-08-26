"""
$(TYPEDEF)

Supertype for renderer-independent plot-definition identities.
"""
abstract type AbstractPlotDefinition end

const LEGEND_OVERFLOW_MODES = (:ellipsis, :show_all)

"""
$(TYPEDEF)

Control the standard legend supplied by the plotting shell.

$(TYPEDFIELDS)
"""
struct LegendDefinition
    "Whether the shell renders a legend."
    enabled::Bool
    "Whether legend entries control the visibility of registered plots."
    interactive::Bool
    "How entries exceeding the available height are handled."
    overflow::Symbol
    function LegendDefinition(
            enabled::Bool,
            interactive::Bool,
            overflow::Symbol
    )
        overflow in LEGEND_OVERFLOW_MODES || throw(
            ArgumentError("legend overflow must be :ellipsis or :show_all"),
        )
        return new(enabled, interactive, overflow)
    end
end

function LegendDefinition(;
        enabled::Bool = true,
        interactive::Bool = true,
        overflow::Symbol = :ellipsis
)
    return LegendDefinition(enabled, interactive, overflow)
end

"""
$(TYPEDEF)

Define one color scale displayed in the standard side dock.

$(TYPEDFIELDS)
"""
struct ColorbarDefinition{C, T}
    "Displayed scale label."
    label::String
    "Ordered colors used by the scale."
    colormap::C
    "Finite, strictly increasing color limits."
    limits::Tuple{Float64, Float64}
    "Tick positions and their displayed labels."
    ticks::T
    function ColorbarDefinition(
            label::String,
            colormap::C,
            limits::Tuple{Float64, Float64},
            ticks::T
    ) where {C, T}
        lower, upper = limits
        isfinite(lower) && isfinite(upper) && lower < upper || throw(
            ArgumentError("colorbar limits must be finite and strictly increasing"),
        )
        ticks isa Tuple && length(ticks) == 2 || throw(
            ArgumentError("colorbar ticks must be a tuple of positions and labels"),
        )
        positions, labels = ticks
        applicable(length, positions) && applicable(length, labels) || throw(
            ArgumentError("colorbar tick positions and labels must be collections"),
        )
        length(positions) == length(labels) || throw(
            DimensionMismatch("colorbar tick positions and labels must have equal lengths"),
        )
        all(position -> position isa Real && isfinite(position), positions) || throw(
            ArgumentError("colorbar tick positions must be finite real values"),
        )
        all(position -> lower <= position <= upper, positions) || throw(
            ArgumentError("colorbar tick positions must lie within the color limits"),
        )
        return new{C, T}(label, colormap, limits, ticks)
    end
end

function ColorbarDefinition(
        label::AbstractString,
        colormap,
        limits::Tuple{<:Real, <:Real},
        ticks
)
    lower, upper = Float64.(limits)
    return ColorbarDefinition(String(label), colormap, (lower, upper), ticks)
end

"""
$(TYPEDEF)

Define the filename and theme used when the standard shell exports SVG.

$(TYPEDFIELDS)
"""
struct ExportDefinition
    "Export theme, either `:default` or `:publication`."
    theme::Symbol
    "Unsanitized base filename."
    name::String
    "Whether to open the completed file with the system application."
    open_file::Bool
    function ExportDefinition(theme::Symbol, name::String, open_file::Bool)
        validate_export_theme(theme)
        isempty(strip(name)) && throw(ArgumentError("export name cannot be empty"))
        return new(theme, name, open_file)
    end
end

function ExportDefinition(;
        theme::Symbol = :default,
        name::AbstractString = "linecablemodels_plot",
        open_file::Bool = true
)
    return ExportDefinition(theme, String(name), open_file)
end

"Supertype for definition-owned controls added to the standard toolbar."
abstract type AbstractWidgetDefinition end

"""
$(TYPEDEF)

Store one completed detached page produced by a domain plot definition.

$(TYPEDFIELDS)
"""
struct PlotPage{K, P}
    "Displayed page title."
    title::String
    "Figure width and height in pixels."
    size::Tuple{Int, Int}
    "Definition-owned semantic page identity."
    key::K
    "Detached definition-owned drawing payload."
    payload::P
    function PlotPage(
            title::String,
            size::Tuple{Int, Int},
            key::K,
            payload::P
    ) where {K, P}
        all(>(0), size) || throw(ArgumentError("page dimensions must be positive"))
        return new{K, P}(title, size, key, payload)
    end
end

function PlotPage(
        title::AbstractString,
        size::Tuple{<:Integer, <:Integer},
        key,
        payload
)
    return PlotPage(String(title), Tuple(Int.(size)), key, payload)
end

"""
$(TYPEDEF)

Store completed detached pages for one plot definition.

`PlotRecipe` is a final rendering product. It does not retain the plotted
source or parsed and resolved request state.

$(TYPEDFIELDS)
"""
struct PlotRecipe{D <: AbstractPlotDefinition, P}
    "Plot definition type."
    definition::Type{D}
    "Concrete tuple or homogeneous collection of completed pages."
    pages::P
end

function _check_recipe_pages(pages)
    all(page -> page isa PlotPage, pages) || throw(
        ArgumentError("PlotRecipe pages must contain only PlotPage values"),
    )
    for first_index in eachindex(pages), second_index in eachindex(pages)

        first_index < second_index || continue
        first_page = pages[first_index]
        second_page = pages[second_index]
        applicable(isempty, first_page.key) && isempty(first_page.key) && continue
        isequal(first_page.key, second_page.key) && throw(
            ArgumentError("PlotRecipe pages must have unique nonempty keys"),
        )
    end
    return pages
end

function PlotRecipe(
        definition::Type{D},
        pages::P
) where {D <: AbstractPlotDefinition, P <: Union{Tuple, AbstractVector}}
    _check_recipe_pages(pages)
    return PlotRecipe{D, P}(definition, pages)
end
