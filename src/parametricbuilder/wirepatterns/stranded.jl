"""
$(TYPEDSIGNATURES)

Estimate hexagonally packed strand patterns for `target_mm2` of metal.

The returned [`WireEstimate`](@ref) retains ranked candidates. When no
candidate reaches the target, the result has `status == :infeasible`.
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
        @warn("The closest stranded pattern exceeds 271 wires.",
            wires=estimate[Val(:match)].wires,)
    return estimate
end

"""
$(TYPEDSIGNATURES)

Return the maximum number of screen wires that fit on `lay_radius`.

`wire_radius` and `lay_radius` are centre-to-centre geometric radii in the same
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
