"Return the number of cable positions in a line-cable system."
function ncables end

"Return the number of distinct positive phases in a line-cable system."
function nphases end

"Return detached preview geometry for one owned cable layer."
function preview_shapes end

"Return the materials represented by one owned cable layer."
function preview_materials end

"""
$(TYPEDSIGNATURES)

Preview a cable design, a vector of cable designs, or a cable system with a
loaded Makie backend.

A vector produces one canvas with one cable per titled panel. Use
`layout = (rows, columns)` to select the grid or omit it for a near-square
layout. Material colorbars span the complete vector and no legend is added.

Requires an explicitly loaded `CairoMakie`, `GLMakie`, or `WGLMakie`
extension.
"""
function preview end

function preview(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling preview.",
    ))
end

"""
$(TYPEDSIGNATURES)

Display the resistivity, permeability, and permittivity colour scales used by
[`preview`](@ref) for visual regression tests. Requires a loaded Makie backend.
"""
function show_material_scale end

function show_material_scale(args...; kwargs...)
    throw(ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling show_material_scale.",
    ))
end
