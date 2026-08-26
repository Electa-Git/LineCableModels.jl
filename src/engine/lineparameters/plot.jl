"""
    plot(parameters[, requests]; kwargs...)
    plot(sources::NamedTuple[, requests]; kwargs...)
    plot(impedance, frequencies[, requests]; kwargs...)
    plot(admittance, frequencies[, requests]; kwargs...)

Plot computed line parameters with a loaded Makie backend. `requests` is an
observable request or an ordered tuple of requests constructed with
`@observe`. Function selectors such as `(R, L, G, C)` and
`(abs, angle)` are normalized at this entry point.

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
