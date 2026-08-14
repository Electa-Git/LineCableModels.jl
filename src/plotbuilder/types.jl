abstract type AbstractPlotSpec end

struct AxisSpec
    dim::Symbol
    quantity::QuantityTag
    units::Units
    label::String
    scale::Symbol
end

struct SeriesSpec{X, Y, Z, A <: NamedTuple}
    kind::Symbol
    xdata::X
    ydata::Y
    zdata::Z
    label::Union{Nothing, String}
    attributes::A
end

function SeriesSpec(kind::Symbol, xdata, ydata, zdata, label; attributes = (;))
    return SeriesSpec(kind, xdata, ydata, zdata, label, attributes)
end

struct ViewSpec{A <: NamedTuple}
    xaxis::Union{Nothing, AxisSpec}
    yaxis::Union{Nothing, AxisSpec}
    zaxis::Union{Nothing, AxisSpec}
    title::String
    series::Vector{SeriesSpec}
    key::NamedTuple
    attributes::A
end

function ViewSpec(xaxis, yaxis, zaxis, title, series, key; attributes = (;))
    return ViewSpec(xaxis, yaxis, zaxis, title, SeriesSpec[series...], key, attributes)
end

struct PageSpec{K <: NamedTuple}
    title::String
    size::Tuple{Int, Int}
    layout::Symbol
    views::Vector{ViewSpec}
    kwargs::K
end

function PageSpec(title, size, layout, views::AbstractVector, kwargs::NamedTuple)
    return PageSpec(title, size, layout, ViewSpec[views...], kwargs)
end

struct RenderSpec{S <: AbstractPlotSpec}
    spec::Type{S}
    figures::Vector{PageSpec}
end

function RenderSpec(spec::Type{S}, figures::AbstractVector) where {S <: AbstractPlotSpec}
    return RenderSpec(spec, PageSpec[figures...])
end

"""
    UIPlot

Hold a backend-neutral render specification together with one built figure,
its panels, controls, and backend context. Line-parameter plotting returns a
`Vector{UIPlot}`; previews and statistical plots return one `UIPlot`.
"""
struct UIPlot{S <: AbstractPlotSpec, F, P, W, C}
    "Complete backend-neutral render specification."
    render::RenderSpec{S}
    "Page represented by this handle."
    page::PageSpec
    "Backend-built figure."
    figure::F
    "Built axes or panels."
    panels::P
    "Interactive control objects keyed by purpose."
    controls::W
    "Active backend and status context."
    context::C
end
