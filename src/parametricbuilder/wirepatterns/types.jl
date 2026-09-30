"""
$(TYPEDEF)

A concentric, hexagonally packed strand pattern.

$(TYPEDFIELDS)
"""
struct HexaPattern{T <: Real}
    "Number of concentric wire layers, including the central wire."
    layers::Int
    "Total wire count."
    wires::Int
    "Diameter of each wire \\[m\\]."
    wire_diameter::T
    "Total metallic cross-section \\[m²\\]."
    total_area::T
    "AWG size label."
    awg::String
end

"""
$(TYPEDEF)

A single-layer wire-screen pattern.

$(TYPEDFIELDS)
"""
struct ScreenPattern{T <: Real}
    "Wire count."
    wires::Int
    "Diameter of each wire \\[m\\]."
    wire_diameter::T
    "Diameter beneath the wire layer \\[m\\]."
    lay_diameter::T
    "Radius to the wire centers \\[m\\]."
    radius::T
    "Total metallic cross-section \\[m²\\]."
    total_area::T
    "Circumferential coverage \\[%\\]."
    coverage::T
    "AWG size or supplied-diameter label."
    awg::String
end

"""
$(TYPEDEF)

Store ranked patterns from a wire-pattern search.

For stranding, `feasible` means at least one pattern reaches the target metal
area; the ranked list can also contain undersized patterns. For screening, it
means at least one pattern satisfies the area, coverage and overshoot bounds,
and only those feasible patterns are retained. When none is feasible, the
geometrically valid alternatives remain ranked and `reasons` describes the limits.
Use `estimate[:closest_area]`, `estimate[:fewest_layers]`, `estimate[:fewest_wires]`, or
`estimate[:smallest_diameter]` to select a pattern.

$(TYPEDFIELDS)
"""
struct WireEstimate{T <: Real, P}
    "Requested metallic cross-section \\[m²\\]."
    target_area::T
    "Ranked wire patterns."
    patterns::Vector{P}
    "Whether at least one pattern meets the estimator's feasibility conditions."
    feasible::Bool
    "Feasibility label, either `:feasible` or `:infeasible`."
    status::Symbol
    "Explanations for unmet constraints."
    reasons::Vector{String}

    function WireEstimate(
            target_area::T,
            patterns::Vector{P},
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
        isempty(patterns) && throw(ArgumentError(
            "a wire estimate must retain at least one pattern",
        ))
        return new{T, P}(target_area, patterns, feasible, status, reasons)
    end
end

Base.length(estimate::WireEstimate) = length(estimate.patterns)
Base.iterate(estimate::WireEstimate, state...) = iterate(estimate.patterns, state...)

_area(pattern::Union{HexaPattern, ScreenPattern}) = pattern.total_area
_diameter(pattern::Union{HexaPattern, ScreenPattern}) = pattern.wire_diameter
_selected(estimate::WireEstimate, key) = argmin(key, estimate.patterns)

Base.getindex(estimate::WireEstimate, selector::Symbol) = estimate[Val(selector)]

function Base.getindex(estimate::WireEstimate, ::Val{:closest_area})
    _selected(estimate, pattern -> (abs(_area(pattern) - estimate.target_area),
        _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate{<:Real, <:HexaPattern}, ::Val{:fewest_layers})
    _selected(estimate, pattern -> (pattern.layers,
        abs(_area(pattern) - estimate.target_area), _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate, ::Val{:fewest_wires})
    _selected(estimate, pattern -> (pattern.wires,
        abs(_area(pattern) - estimate.target_area), _diameter(pattern)))
end
function Base.getindex(estimate::WireEstimate, ::Val{:smallest_diameter})
    _selected(estimate, pattern -> (_diameter(pattern),
        abs(_area(pattern) - estimate.target_area), pattern.wires))
end

function Base.getindex(::WireEstimate, ::Val{selector}) where {selector}
    throw(ArgumentError(
        "unknown wire-estimate selector :$selector; use :closest_area, :fewest_layers, " *
        ":fewest_wires, or :smallest_diameter",
    ))
end
