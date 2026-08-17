module WirePatterns

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import ...LineCableModels: maxfill, nominal

export WireEstimate, make_stranded, make_screened

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

Best-effort result of a wire-pattern search.

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

_wire_area(diameter::Real) = (one(diameter) * pi / 4) * diameter^2

const _AWG_BASE = 92
const _D0_MM = 0.127
const _AREA0_MM2 = 0.012668

function awg_to_d_mm(number::Real)
    oftype(float(number), _D0_MM) *
    oftype(float(number), _AWG_BASE) ^
    ((oftype(float(number), 36) - number) / oftype(float(number), 39))
end
function awg_to_area_mm2(number::Real)
    oftype(float(number), _AREA0_MM2) *
    oftype(float(number), _AWG_BASE) ^
    ((oftype(float(number), 36) - number) / oftype(float(number), 19.5))
end

function d_mm_to_awg(diameter_mm::Real)
    oftype(float(diameter_mm), 36) -
    oftype(float(diameter_mm), 39) *
    log(diameter_mm / oftype(float(diameter_mm), _D0_MM)) /
    log(oftype(float(diameter_mm), _AWG_BASE))
end
function area_mm2_to_awg(area_mm2::Real)
    oftype(float(area_mm2), 36) -
    oftype(float(area_mm2), 19.5) *
    log(area_mm2 / oftype(float(area_mm2), _AREA0_MM2)) /
    log(oftype(float(area_mm2), _AWG_BASE))
end

function awg_label(number::Integer)
    number == -3 && return "0000 (4/0)"
    number == -2 && return "000 (3/0)"
    number == -1 && return "00 (2/0)"
    number == 0 && return "0 (1/0)"
    return string(number)
end

function awg_sizes(::Type{T}, nmin::Integer = -3, nmax::Integer = 40) where {T <: Real}
    nmin <= nmax || throw(ArgumentError("nmin must not exceed nmax"))
    return [(awg_label(number), convert(T, awg_to_d_mm(number)) / T(1000))
            for number in nmin:nmax]
end

awg_sizes(nmin::Integer = -3, nmax::Integer = 40) = awg_sizes(Float64, nmin, nmax)

"Apply a fill factor to solid area to approximate stranded metallic area."
function stranded_area_mm2(number::Real; fill_factor::Real = 0.94)
    factor, area = promote(float(fill_factor), float(awg_to_area_mm2(number)))
    return factor * area
end

const _WIRE_RULES = Tuple{Int, Int, Union{Int, Nothing}}[
    (10, 6, 7), (16, 6, 7), (25, 6, 7), (35, 6, 7),
    (50, 6, 19), (70, 12, 19), (95, 15, 19), (120, 15, 37),
    (150, 15, 37), (185, 30, 37), (240, 30, 37), (300, 30, 61),
    (400, 53, 61), (500, 53, 61), (630, 53, 91), (800, 53, 91),
    (1000, 53, 91)
]

_hex_wires(layers::Int) = 1 + 3layers * (layers - 1)

function _allowed_wires(target_mm2::Real)
    for (threshold, minimum, maximum) in _WIRE_RULES
        target_mm2 <= threshold && return minimum, maximum
    end
    return 53, nothing
end

function _allowed_wires(wires::Int, minimum::Int, maximum::Union{Int, Nothing})
    return maximum === nothing ? wires >= minimum : minimum <= wires <= maximum
end

"""
$(TYPEDSIGNATURES)

Estimate hexagonally packed strand patterns for `target_mm2` of metal.

The returned [`WireEstimate`](@ref) retains ranked candidates. A valid search
that cannot reach the target returns `status == :infeasible` rather than
throwing.
"""
function make_stranded(target_mm2::Real; nmin::Integer = -3, nmax::Integer = 40)
    target_mm2 > zero(target_mm2) || throw(DomainError(
        target_mm2, "target cross-section must be positive"
    ))
    nmin <= nmax || throw(ArgumentError("nmin must not exceed nmax"))

    target_value = float(target_mm2)
    T = typeof(target_value)
    target = target_value * T(1e-6)
    minimum, maximum_wires = _allowed_wires(target_value)
    candidates = HexaPattern{T}[]

    for (awg, diameter) in awg_sizes(T, nmin, nmax)
        area = _wire_area(diameter)
        for layers in 1:300
            wires = _hex_wires(layers)
            _allowed_wires(wires, minimum, maximum_wires) && push!(
                candidates,
                HexaPattern(layers, wires, diameter, wires * area, awg)
            )
            maximum_wires !== nothing && wires > maximum_wires && break
        end
    end
    isempty(candidates) && throw(ArgumentError(
        "the AWG range and permitted wire counts produced no candidates",
    ))

    sort!(candidates;
        by = pattern -> (
            abs(pattern.total_area_m2 - target), pattern.wire_diameter_m,
            pattern.layers
        ))
    feasible = any(pattern -> pattern.total_area_m2 >= target, candidates)
    reasons = feasible ? String[] :
              [
        "no permitted strand pattern reaches the requested metallic area",
    ]
    estimate = WireEstimate(
        target, candidates, feasible, feasible ? :feasible : :infeasible, reasons
    )
    estimate[Val(:match)].wires > 271 &&
        @warn("The closest stranded pattern has a very high wire count.",
            wires=estimate[Val(:match)].wires,)
    return estimate
end

"""
$(TYPEDSIGNATURES)

Return the maximum number of screen wires that fit on `lay_radius`.

`wire_radius` and `lay_radius` are center-to-center geometric radii in the same
unit. `gap_frac` adds a fractional clearance between adjacent wires.
"""
function maxfill(
        ::Type{ScreenPattern},
        lay_radius::Real,
        wire_radius::Real;
        gap_frac::Real = 0
)
    lay_radius > zero(lay_radius) || throw(DomainError(
        lay_radius, "lay radius must be positive"
    ))
    wire_radius > zero(wire_radius) || throw(DomainError(
        wire_radius, "wire radius must be positive"
    ))
    gap_frac >= zero(gap_frac) || throw(DomainError(
        gap_frac, "gap fraction must be nonnegative"
    ))
    ratio = wire_radius * (one(gap_frac) + gap_frac) / lay_radius
    nominal_ratio = float(nominal(ratio))
    zero(nominal_ratio) < nominal_ratio < one(nominal_ratio) || return 0
    count = pi / asin(nominal_ratio)
    return max(0, floor(Int, count + 8eps(count)))
end

function _screen_candidate(
        wires::Int,
        diameter::T,
        lay_diameter::T,
        sine_angle::T,
        awg::String
) where {T <: Real}
    area = wires * _wire_area(diameter)
    coverage = T(100) * wires * diameter /
               (T(pi) * lay_diameter * sine_angle)
    return ScreenPattern(
        wires, diameter, lay_diameter, (lay_diameter + diameter) / T(2),
        area, coverage, awg
    )
end

function _screen_feasible(
        pattern::ScreenPattern,
        target,
        coverage_min,
        coverage_max,
        overshoot_max
)
    enough_area = pattern.total_area_m2 >= target
    coverage_ok = coverage_min <= pattern.coverage_pct <= coverage_max
    T = typeof(target)
    overshoot = T(100) * (pattern.total_area_m2 / target - one(T))
    overshoot_ok = !isfinite(overshoot_max) || overshoot <= overshoot_max
    return enough_area && coverage_ok && overshoot_ok
end

"""
$(TYPEDSIGNATURES)

Estimate single-layer wire-screen patterns for the required area and laying
diameter, both expressed in millimetre-based input units.

Geometrically valid candidates must meet the requested area, coverage, and
overshoot bounds to make the result feasible. Otherwise the returned
[`WireEstimate`](@ref) contains ranked best-effort candidates and reasons.
"""
function make_screened(
        required_area_mm2::Real,
        lay_diameter_mm::Real;
        alpha_deg::Real = 15,
        coverage_min_pct::Real = 85,
        gap_frac::Real = 0,
        min_wires::Integer = 6,
        extra_span::Integer = 8,
        nmin::Integer = -3,
        nmax::Integer = 40,
        coverage_max_pct::Real = 100,
        max_overshoot_pct::Real = 10,
        custom_diameters_mm::AbstractVector{<:Real} = Float64[]
)
    required_area_mm2 > zero(required_area_mm2) || throw(DomainError(
        required_area_mm2, "required cross-section must be positive"
    ))
    lay_diameter_mm > zero(lay_diameter_mm) || throw(DomainError(
        lay_diameter_mm, "laying diameter must be positive"
    ))
    zero(coverage_min_pct) < coverage_min_pct <= 100 || throw(DomainError(
        coverage_min_pct, "minimum coverage must be in (0, 100]"
    ))
    coverage_max_pct >= coverage_min_pct || throw(DomainError(
        coverage_max_pct, "maximum coverage must not be below its minimum"
    ))
    max_overshoot_pct >= zero(max_overshoot_pct) || throw(DomainError(
        max_overshoot_pct, "maximum overshoot must be nonnegative"
    ))
    gap_frac >= zero(gap_frac) || throw(DomainError(
        gap_frac, "gap fraction must be nonnegative"
    ))
    min_wires >= 1 || throw(DomainError(min_wires, "minimum wires must be positive"))
    extra_span >= 0 || throw(DomainError(extra_span, "extra span must be nonnegative"))
    nmin <= nmax || throw(ArgumentError("nmin must not exceed nmax"))
    all(>(0), custom_diameters_mm) || throw(DomainError(
        custom_diameters_mm, "custom wire diameters must be positive"
    ))

    custom_types = isempty(custom_diameters_mm) ? () : (eltype(custom_diameters_mm),)
    T = promote_type(
        typeof(float(required_area_mm2)), typeof(float(lay_diameter_mm)),
        typeof(float(alpha_deg)), typeof(float(coverage_min_pct)),
        typeof(float(coverage_max_pct)), typeof(float(max_overshoot_pct)),
        typeof(float(gap_frac)), custom_types...
    )
    target = convert(T, required_area_mm2) * T(1e-6)
    lay_diameter = convert(T, lay_diameter_mm) * T(1e-3)
    angle = deg2rad(convert(T, alpha_deg))
    sine_angle = sin(angle)
    sine_angle > zero(T) || throw(DomainError(
        alpha_deg, "lay angle must have a positive sine"
    ))
    coverage_min = convert(T, coverage_min_pct)
    coverage_max = convert(T, coverage_max_pct)
    overshoot_max = convert(T, max_overshoot_pct)
    gap = convert(T, gap_frac)

    sizes = awg_sizes(T, nmin, nmax)
    append!(sizes,
        [("custom($(round(diameter; digits=3)) mm)", convert(T, diameter) / T(1000))
         for diameter in custom_diameters_mm])

    geometric = ScreenPattern{T}[]
    for (awg, diameter) in sizes
        area = _wire_area(diameter)
        required_by_area = ceil(Int, target / area)
        required_by_coverage = ceil(
            Int, coverage_min * T(pi) * lay_diameter * sine_angle /
                 (T(100) * diameter)
        )
        required = max(Int(min_wires), required_by_area, required_by_coverage)
        lay_radius = (lay_diameter + diameter) / T(2)
        maximum_wires = maxfill(
            ScreenPattern, lay_radius, diameter / T(2); gap_frac = gap
        )
        upper = min(maximum_wires, required + Int(extra_span))
        if upper >= min_wires
            append!(geometric,
                [_screen_candidate(wires, diameter, lay_diameter, sine_angle, awg)
                 for wires in Int(min_wires):upper])
        elseif maximum_wires > 0
            push!(geometric, _screen_candidate(
                maximum_wires, diameter, lay_diameter, sine_angle, awg
            ))
        end
    end
    isempty(geometric) && throw(ArgumentError(
        "the supplied diameters cannot form even one wire on this laying radius",
    ))

    feasible_candidates = filter(
        pattern -> _screen_feasible(
            pattern, target, coverage_min, coverage_max, overshoot_max
        ),
        geometric
    )
    feasible = !isempty(feasible_candidates)
    candidates = feasible ? feasible_candidates : geometric
    sort!(candidates;
        by = pattern -> (
            abs(pattern.total_area_m2 - target), pattern.wires,
            pattern.wire_diameter_m
        ))

    reasons = String[]
    if !feasible
        maximum_area = Base.maximum(pattern.total_area_m2 for pattern in geometric)
        maximum_coverage = Base.maximum(pattern.coverage_pct for pattern in geometric)
        maximum_area < target && push!(
            reasons, "available single-layer patterns do not reach the requested area"
        )
        maximum_coverage < coverage_min && push!(
            reasons, "available single-layer patterns do not reach minimum coverage"
        )
        isempty(reasons) && push!(
            reasons, "coverage or overshoot limits reject every geometric candidate"
        )
    end
    return WireEstimate(
        target, candidates, feasible, feasible ? :feasible : :infeasible, reasons
    )
end

end
