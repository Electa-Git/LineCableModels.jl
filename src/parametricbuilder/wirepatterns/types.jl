"""
$(TYPEDEF)

A candidate concentric, hexagonally packed stranded conductor.

$(TYPEDFIELDS)
"""
struct HexaPattern{T <: Real}
    layers::Int
    wires::Int
    wire_diameter_m::T
    total_area_m2::T
    awg::String
end

"""
$(TYPEDEF)

A candidate single-layer wire screen.

$(TYPEDFIELDS)
"""
struct ScreenPattern{T <: Real}
    wires::Int
    wire_diameter_m::T
    lay_diameter_m::T
    radius_m::T
    total_area_m2::T
    coverage_pct::T
    awg::String
end

"""
$(TYPEDEF)

Store ranked candidates from a wire-pattern search.

`feasible` states whether the retained candidates satisfy every requested
constraint. An infeasible result still contains the closest candidates that
the search could construct and explains the limiting constraints in `reasons`.
Use `estimate[Val(:match)]`, `estimate[Val(:layers)]`,
`estimate[Val(:wires)]`, or `estimate[Val(:diameter)]` to select a candidate.

$(TYPEDFIELDS)
"""
struct WireEstimate{T <: Real, P}
    target::T
    candidates::Vector{P}
    feasible::Bool
    status::Symbol
    reasons::Vector{String}

    function WireEstimate(
            target::T,
            candidates::Vector{P},
            feasible::Bool,
            status::Symbol,
            reasons::Vector{String}
    ) where {T <: Real, P}
        status in (:feasible, :infeasible) || throw(ArgumentError(
            "wire-estimate status must be :feasible or :infeasible",
        ))
        feasible == (status === :feasible) || throw(ArgumentError(
            "wire-estimate feasibility and status disagree",
        ))
        isempty(candidates) && throw(ArgumentError(
            "a wire estimate must retain at least one candidate",
        ))
        return new{T, P}(target, candidates, feasible, status, reasons)
    end
end

Base.length(estimate::WireEstimate) = length(estimate.candidates)
Base.iterate(estimate::WireEstimate, state...) = iterate(estimate.candidates, state...)

_area(pattern::Union{HexaPattern, ScreenPattern}) = pattern.total_area_m2
_diameter(pattern::Union{HexaPattern, ScreenPattern}) = pattern.wire_diameter_m
_selected(estimate::WireEstimate, key) = argmin(key, estimate.candidates)

function Base.getindex(estimate::WireEstimate, ::Val{:match})
    _selected(estimate, pattern -> (abs(_area(pattern) - estimate.target),
        _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate{<:Real, <:HexaPattern}, ::Val{:layers})
    _selected(estimate, pattern -> (pattern.layers,
        abs(_area(pattern) - estimate.target), _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate, ::Val{:wires})
    _selected(estimate, pattern -> (pattern.wires,
        abs(_area(pattern) - estimate.target), _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate, ::Val{:diameter})
    _selected(estimate, pattern -> (_diameter(pattern),
        abs(_area(pattern) - estimate.target), pattern.wires))
end

function Base.getindex(::WireEstimate, ::Val{selector}) where {selector}
    throw(ArgumentError(
        "unknown wire-estimate selector :$selector; use :match, :layers, " *
        ":wires, or :diameter",
    ))
end
