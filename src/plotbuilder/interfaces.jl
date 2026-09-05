"""
    plot(source[, selection]; kwargs...)

Create a compact native Makie plot for a supported result or observation
publication. The optional Makie extension infers distinct physical quantities
from the requested observables; `layout` controls only where those inferred
axes are placed.
"""
function plot end

function plot(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling plot.",
    ))
end

"""
    preview(source; kwargs...)

Preview a cable design, a collection of cable designs, or a cable system with
a loaded Makie backend.
"""
function preview end

function preview(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling preview.",
    ))
end

"""
    show_material_scale(; kwargs...)

Display the three independently defined material colour schemes as a compact
reference figure. Use [`materialscale!`](@ref) to place any one scheme in a
caller-owned Makie layout.
"""
function show_material_scale end

function show_material_scale(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling show_material_scale.",
    ))
end

"""
    export_svg(plot::UIPlot; path=nothing, theme=nothing, open_file=nothing)

Save the current live Makie figure in `plot` as SVG through CairoMakie and
return the absolute output path. `theme` may be `:default` or `:publication`.
"""
function export_svg end

"""
    figurelegend!(plot::UIPlot; position=:right, anchor=:rt, kwargs...)

Create or replace the figure-scoped native Makie legend from the semantic
groups retained by the shell. All remaining keywords are forwarded to Makie's
`Legend`.
"""
function figurelegend! end

"""
    panellegend!(plot::UIPlot, panel; kwargs...)

Create or replace a native Makie legend scoped to one logical plot panel.
`panel` may be the stable panel identity returned by a recipe or its compatible
grid position.
"""
function panellegend! end

"""
    figuretitle!(plot::UIPlot, title; kwargs...)

Create, replace, or remove the figure-wide native Makie title. Pass `nothing`
to remove it.
"""
function figuretitle! end

"""
    paneltitle!(plot::UIPlot, panel, title)

Set the native axis title for one logical plot panel. Pass `nothing` to clear
it.
"""
function paneltitle! end

"""
    plotwindow(callback; title, figure_title=nothing, size=(800, 400), kwargs...)

Build the standard Makie shell, pass its caller-owned content `GridLayout` to
`callback`, and return a [`UIPlot`](@ref). The callback uses ordinary Makie and
is not constrained by a renderer-independent plot specification.
"""
function plotwindow end

"""
    materialcolors(property, [range]; alpha=1.0)

Construct one reusable material colour scheme. Palette selection is separate
from [`materialscale!`](@ref), which only renders a supplied scheme.
"""
function materialcolors end

"""
    materialscale!(position, scheme; kwargs...)

Place one native Makie color scale at `position`. `scheme` supplies one label,
colormap, limits, and tick definition. The primitive never chooses or combines
material properties.
"""
function materialscale! end
