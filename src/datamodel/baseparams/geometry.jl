"""
$(TYPEDSIGNATURES)

Return mean diameter, pitch length, and helical overlength for a radial layer.

```math
\\lambda=L_p/D_e,\\qquad
k=\\sqrt{1+(\\pi D_e/L_p)^2}.
```

# Notes

The lay ratio follows EN 50182. A zero lay ratio represents a straight layer
and therefore has unit overlength.
"""
function helix(r_in::Real, r_ex::Real, lay_ratio::Real)
    rin, rex, ratio = promote(float(r_in), float(r_ex), float(lay_ratio))
    all(isfinite, (rin, rex, ratio)) || throw(DomainError(
        (rin, rex, ratio), "helix geometry must be finite"
    ))
    rin >= zero(rin) || throw(DomainError(rin, "inner radius must be nonnegative"))
    rex >= rin ||
        throw(DomainError(rex, "outer radius must not be smaller than inner radius"))
    ratio >= zero(ratio) || throw(DomainError(ratio, "lay ratio must be nonnegative"))
    mean_diameter = rin + rex
    pitch_length = ratio * mean_diameter
    overlength = iszero(pitch_length) ? one(pitch_length) :
                 sqrt(one(pitch_length) +
                      (oftype(pitch_length, π) * mean_diameter / pitch_length)^2)
    return mean_diameter, pitch_length, overlength
end

"""
$(TYPEDSIGNATURES)

Return the centers of `num_wires` equally spaced circular wires. For wire
``i=0,\\ldots,N-1`` the coordinates are
``C+r_l(\\cos(2\\pi i/N),\\sin(2\\pi i/N))``.
"""
function wire_coordinates(
        num_wires::Integer,
        radius_wire::Real,
        r_in::Real;
        C = nothing
)
    num_wires >= 0 || throw(DomainError(num_wires, "number of wires must be nonnegative"))
    center = C === nothing ? (zero(radius_wire), zero(radius_wire)) : C
    rw, rin, cx, cy = promote(
        float(radius_wire), float(r_in), float(center[1]), float(center[2])
    )
    all(isfinite, (rw, rin, cx, cy)) || throw(DomainError(
        (rw, rin, cx, cy), "wire geometry must be finite"
    ))
    rw > zero(rw) || throw(DomainError(rw, "wire radius must be positive"))
    rin >= zero(rin) || throw(DomainError(rin, "inner radius must be nonnegative"))
    radius = num_wires == 1 ? zero(rin) : rin + rw
    num_wires == 0 && return Tuple{typeof(rin), typeof(rin)}[]
    step = 2 * (one(rin) * π) / num_wires
    return [(cx + radius * cos(index * step), cy + radius * sin(index * step))
            for index in 0:(num_wires - 1)]
end

function wire_coordinates(num_wires::Integer, radius_wire::Real, r_in::Real, C::Tuple)
    wire_coordinates(num_wires, radius_wire, r_in; C)
end

"Return GMD integration elements `(coordinates, radius, area)` for a cable part."
function gmd_elements(part::AbstractCablePart)
    T = eltype(part)
    return ([(zero(T), zero(T))], part.r_ex, part.cross_section)
end

"""
$(TYPEDSIGNATURES)

Calculate the area-weighted geometric mean distance between two cable parts:

```math
\\log GMD=\\frac{\\sum_i\\sum_j s_i s_j\\log d_{ij}}
                 {\\sum_i\\sum_j s_i s_j}.
```

# Notes

A sum of logarithms is used instead of a product, preventing numerical
underflow and overflow for large wire arrays. Coincident centers use the
outermost element radius.
"""
function gmd(first_part::AbstractCablePart, second_part::AbstractCablePart)
    coords1, radius1, area1 = gmd_elements(first_part)
    coords2, radius2, area2 = gmd_elements(second_part)
    T = promote_type(
        typeof(float(radius1)), typeof(float(radius2)),
        typeof(float(area1)), typeof(float(area2))
    )
    log_sum = zero(T)
    weights = zero(T)
    weight = convert(T, area1) * convert(T, area2)
    for first in coords1, second in coords2

        distance = hypot(first[1] - second[1], first[2] - second[2])
        effective = iszero(distance) ? max(radius1, radius2) : distance
        log_sum += weight * log(effective)
        weights += weight
    end
    return exp(log_sum / weights)
end

"""
$(TYPEDSIGNATURES)

Combine two conductor zones using the multizone GMR expression
``GMR=GMR_1^{\\beta^2}GMR_2^{(1-\\beta)^2}GMD^{2\\beta(1-\\beta)}``.
"""
function equivalent_gmr(existing::AbstractCablePart, new_layer::AbstractCablePart)
    beta = existing.cross_section / (existing.cross_section + new_layer.cross_section)
    distance = gmd(existing, new_layer)
    return existing.gmr^(beta^2) * new_layer.gmr^((one(beta) - beta)^2) *
           distance^(2 * beta * (one(beta) - beta))
end

"""
$(TYPEDSIGNATURES)

Return the helical-solenoid permeability factor
``1+2\\pi^2N^2(r_i^2-r_c^2)/\\log(r_i/r_c)``.
"""
function solenoid_factor(num_turns::Real, r_con::Real, r_ins::Real)
    turns, conductor, insulator = promote(
        float(num_turns), float(r_con), float(r_ins)
    )
    isnan(turns) && return one(turns)
    turns >= zero(turns) || throw(DomainError(turns, "number of turns must be nonnegative"))
    conductor >= zero(conductor) || throw(DomainError(
        conductor, "conductor radius must be nonnegative"
    ))
    insulator > conductor || throw(DomainError(
        insulator, "insulator radius must exceed conductor radius"
    ))
    pi_value = one(turns) * π
    return one(turns) +
           2 * pi_value^2 * turns^2 *
           (insulator^2 - conductor^2) / log(insulator / conductor)
end
