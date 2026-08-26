"Return the number of cable positions in a line-cable system."
function ncables end

"Return the number of distinct positive phases in a line-cable system."
function nphases end

"Return the base cable parameters calculated for a materialised design."
function base_parameters end

"""
$(TYPEDSIGNATURES)

Preview a cable design or cable system with a loaded Makie backend.

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
