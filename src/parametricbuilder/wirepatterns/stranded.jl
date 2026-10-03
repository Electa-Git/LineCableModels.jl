"""
$(TYPEDSIGNATURES)

Estimate hexagonally packed strand patterns for a metallic cross-section.

# Arguments

- `target_area`: required cross-sectional area of metal \\[mm²\\]. Stored areas use \\[m²\\].

# Keywords

- `awg_min=-3`, `awg_max=40`: inclusive AWG-number limits.

# Returns

The returned [`WireEstimate`](@ref) retains ranked patterns, including those
below the target area. When no pattern reaches that area, the result has
`status == :infeasible`.
"""
function estimate_stranding(target_area::Real; awg_min::Integer = -3, awg_max::Integer = 40)
    target_area > zero(target_area) || throw(DomainError(
        target_area, "target cross-section must be positive"
    ))
    awg_min <= awg_max || throw(ArgumentError("awg_min must not exceed awg_max"))

    target_value = float(target_area)
    T = typeof(target_value)
    target_area = target_value * T(1e-6)
    minimum, maximum_wires = _allowed_wires(target_value)
    patterns = HexaPattern{T}[]

    for (awg, diameter) in awg_sizes(T, awg_min, awg_max)
        area = _wire_area(diameter)
        for layers in 1:300
            wires = _hex_wires(layers)
            _allowed_wires(wires, minimum, maximum_wires) && push!(
                patterns,
                HexaPattern(layers, wires, diameter, wires * area, awg)
            )
            maximum_wires !== nothing && wires > maximum_wires && break
        end
    end
    isempty(patterns) && throw(ArgumentError(
        "the AWG range and permitted wire counts produced no patterns",
    ))

    sort!(patterns;
        by = pattern -> (
            abs(pattern.total_area - target_area), pattern.wire_diameter,
            pattern.layers
        ))
    feasible = any(pattern -> pattern.total_area >= target_area, patterns)
    reasons = feasible ? String[] :
              [
        "no permitted strand pattern reaches the requested metallic area",
    ]
    estimate = WireEstimate(
        target_area, patterns, feasible, feasible ? :feasible : :infeasible, reasons
    )
    estimate[:closest_area].wires > 271 &&
        @warn("The closest stranded pattern exceeds 271 wires.",
            wires=estimate[:closest_area].wires,)
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
