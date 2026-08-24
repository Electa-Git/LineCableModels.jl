"""
    plot(parameters[, quantities]; kwargs...)
    plot(first, second, rest...; legend, quantities=(), kwargs...)
    plot(impedance, frequencies[, quantities]; kwargs...)
    plot(admittance, frequencies[, quantities]; kwargs...)

Plot computed line parameters with a loaded Makie backend. `quantities` is a
tuple of accessors such as `(R, L, G, C)` or `(abs, angle)`.

With two or more positional [`LineParameters`](@ref) results, create one
matrix-grid page per selected quantity. Grid position `(i, j)` overlays the
corresponding frequency series from every result. The required `legend` tuple
must contain one label per result.

Without an explicit selection, [`LineParameters`](@ref) produces separate
series-impedance and shunt-admittance figures. Each figure places the real part
on the left and the imaginary part on the right. Every selected matrix element
is represented by one data series.

Requires an explicitly loaded `CairoMakie`, `GLMakie`, or `WGLMakie`
extension.

# Returns

- A `Vector{UIPlot}` containing one figure for each selected matrix family.
"""
function plot end

function plot(args...; kwargs...)
    throw(
        ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling plot.",
    ),
    )
end
