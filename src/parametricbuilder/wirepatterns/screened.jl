function ScreenPattern(
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
        target_area,
        coverage_min,
        coverage_max,
        overshoot_max
)
    enough_area = pattern.total_area >= target_area
    coverage_ok = coverage_min <= pattern.coverage <= coverage_max
    T = typeof(target_area)
    overshoot = T(100) * (pattern.total_area / target_area - one(T))
    overshoot_ok = !isfinite(overshoot_max) || overshoot <= overshoot_max
    return enough_area && coverage_ok && overshoot_ok
end

"""
$(TYPEDSIGNATURES)

Estimate single-layer wire-screen patterns.

# Arguments

- `target_area`: required metallic cross-section \\[mm²\\].
- `lay_diameter`: diameter beneath the wire layer \\[mm\\].

# Keywords

- `lay_angle=15`: wire lay angle \\[degrees\\], with positive sine.
- `coverage_min=85`, `coverage_max=100`: circumferential coverage bounds \\[%\\].
- `gap_frac=0`: additional wire clearance as a fraction of wire diameter.
- `min_wires=6`: minimum wire count.
- `extra_span=8`: additional wire counts examined above the required count.
- `awg_min=-3`, `awg_max=40`: inclusive AWG-number limits.
- `max_area_overshoot=10`: maximum excess metallic area relative to the target \\[%\\].
- `wire_diameters=[]`: additional wire diameters \\[mm\\].

# Returns

A [`WireEstimate`](@ref) whose stored lengths and areas use \\[m\\] and \\[m²\\].

Geometrically valid patterns must meet the requested area, coverage, and
overshoot bounds to make the result feasible. Otherwise the returned
[`WireEstimate`](@ref) contains ranked best-effort patterns and reasons.
"""
function estimate_screen(
        target_area::Real,
        lay_diameter::Real;
        lay_angle::Real = 15,
        coverage_min::Real = 85,
        gap_frac::Real = 0,
        min_wires::Integer = 6,
        extra_span::Integer = 8,
        awg_min::Integer = -3,
        awg_max::Integer = 40,
        coverage_max::Real = 100,
        max_area_overshoot::Real = 10,
        wire_diameters::AbstractVector{<:Real} = Float64[]
)
    target_area > zero(target_area) || throw(DomainError(
        target_area, "required cross-section must be positive"
    ))
    lay_diameter > zero(lay_diameter) || throw(DomainError(
        lay_diameter, "laying diameter must be positive"
    ))
    zero(coverage_min) < coverage_min <= 100 || throw(DomainError(
        coverage_min, "minimum coverage must be in (0, 100]"
    ))
    coverage_max >= coverage_min || throw(DomainError(
        coverage_max, "maximum coverage must not be below its minimum"
    ))
    max_area_overshoot >= zero(max_area_overshoot) || throw(DomainError(
        max_area_overshoot, "maximum overshoot must be nonnegative"
    ))
    gap_frac >= zero(gap_frac) || throw(DomainError(
        gap_frac, "gap fraction must be nonnegative"
    ))
    min_wires >= 1 || throw(DomainError(min_wires, "minimum wires must be positive"))
    extra_span >= 0 || throw(DomainError(extra_span, "extra span must be nonnegative"))
    awg_min <= awg_max || throw(ArgumentError("awg_min must not exceed awg_max"))
    all(>(0), wire_diameters) || throw(DomainError(
        wire_diameters, "custom wire diameters must be positive"
    ))

    custom_types = isempty(wire_diameters) ? () : (eltype(wire_diameters),)
    T = promote_type(
        typeof(float(target_area)), typeof(float(lay_diameter)),
        typeof(float(lay_angle)), typeof(float(coverage_min)),
        typeof(float(coverage_max)), typeof(float(max_area_overshoot)),
        typeof(float(gap_frac)), custom_types...
    )
    target_area = convert(T, target_area) * T(1e-6)
    lay_diameter = convert(T, lay_diameter) * T(1e-3)
    angle = deg2rad(convert(T, lay_angle))
    sine_angle = sin(angle)
    sine_angle > zero(T) || throw(DomainError(
        lay_angle, "lay angle must have a positive sine"
    ))
    coverage_min = convert(T, coverage_min)
    coverage_max = convert(T, coverage_max)
    overshoot_max = convert(T, max_area_overshoot)
    gap = convert(T, gap_frac)

    sizes = awg_sizes(T, awg_min, awg_max)
    append!(sizes,
        [("custom($(round(diameter; digits=3)) mm)", convert(T, diameter) / T(1000))
         for diameter in wire_diameters])

    geometric = ScreenPattern{T}[]
    for (awg, diameter) in sizes
        area = _wire_area(diameter)
        required_by_area = ceil(Int, target_area / area)
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
                [ScreenPattern(wires, diameter, lay_diameter, sine_angle, awg)
                 for wires in Int(min_wires):upper])
        elseif maximum_wires > 0
            push!(geometric, ScreenPattern(
                maximum_wires, diameter, lay_diameter, sine_angle, awg
            ))
        end
    end
    isempty(geometric) && throw(ArgumentError(
        "the supplied diameters cannot form even one wire on this laying radius",
    ))

    feasible_patterns = filter(
        pattern -> _screen_feasible(
            pattern, target_area, coverage_min, coverage_max, overshoot_max
        ),
        geometric
    )
    feasible = !isempty(feasible_patterns)
    patterns = feasible ? feasible_patterns : geometric
    sort!(patterns;
        by = pattern -> (
            abs(pattern.total_area - target_area), pattern.wires,
            pattern.wire_diameter
        ))

    reasons = String[]
    if !feasible
        maximum_area = Base.maximum(pattern.total_area for pattern in geometric)
        maximum_coverage = Base.maximum(pattern.coverage for pattern in geometric)
        maximum_area < target_area && push!(
            reasons, "available single-layer patterns do not reach the requested area"
        )
        maximum_coverage < coverage_min && push!(
            reasons, "available single-layer patterns do not reach minimum coverage"
        )
        isempty(reasons) && push!(
            reasons, "coverage or overshoot limits reject every geometric pattern"
        )
    end
    return WireEstimate(
        target_area, patterns, feasible, feasible ? :feasible : :infeasible, reasons
    )
end
