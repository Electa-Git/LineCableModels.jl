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
