"""
$(TYPEDSIGNATURES)

Calculate the GMR of a circular wire array:
``GMR=(r_g N a^{N-1})^{1/N}``, with
``r_g=r_w\\exp(-\\mu_r/4)``.

`lay_radius` is the radius `a` of the strand-center circle in meters;
`count` is the positive integer strand count `N`; `wire_radius` is the
positive strand radius `r_w` in meters; and `mu_r` is the positive,
dimensionless relative permeability. The lay radius must be positive for
multiple wires and may be zero for one wire. The returned geometric mean
radius (GMR) is in meters.
"""
function strand_gmr(lay_radius::Real, count::Integer, wire_radius::Real, mu_r::Real)
    count > 0 || throw(DomainError(count, "wire count must be positive"))
    radius, wire, permeability = promote(
        float(lay_radius), float(wire_radius), float(mu_r)
    )
    wire > zero(wire) || throw(DomainError(wire, "wire radius must be positive"))
    radius >= zero(radius) || throw(DomainError(radius, "lay radius must be nonnegative"))
    (count == 1 || radius > zero(radius)) || throw(DomainError(
        radius, "lay radius must be positive for a multi-wire layer"
    ))
    permeability > zero(permeability) || throw(DomainError(
        permeability, "relative permeability must be positive"
    ))
    wire_gmr = wire * exp(-permeability / 4)
    return exp(log(wire_gmr * count * radius^(count - 1)) / count)
end

"""
$(TYPEDSIGNATURES)

Calculate the GMR of identical circular strands from their resolved centers:

```math
\\log GMR=\\frac{1}{N^2}\\left[
N\\log r_g+2\\sum_{i=1}^{N}\\sum_{j=i+1}^{N}\\log d_{ij}
\\right],\\qquad
r_g=r_w\\exp(-\\mu_r/4).
```

# Arguments

- `coordinates`: Strand-center coordinates \\[m\\].
- `wire_radius`: Strand radius \\[m\\].
- `mu_r`: Relative permeability \\[dimensionless\\].

# Returns

- Geometric mean radius of the complete strand set \\[m\\].
"""
function strand_gmr(coordinates, wire_radius::Real, mu_r::Real)
    isempty(coordinates) && throw(ArgumentError(
        "strand GMR requires at least one centre coordinate"
    ))
    wire, permeability = promote(float(wire_radius), float(mu_r))
    wire > zero(wire) || throw(DomainError(
        wire, "wire radius must be positive"
    ))
    permeability > zero(permeability) || throw(DomainError(
        permeability, "relative permeability must be positive"
    ))
    count = length(coordinates)
    logarithmic_sum = count * log(tubular_gmr(wire, zero(wire), permeability))
    for left in eachindex(coordinates)
        for right in (left + 1):length(coordinates)
            distance = hypot(
                coordinates[left][1] - coordinates[right][1],
                coordinates[left][2] - coordinates[right][2]
            )
            distance > zero(distance) || throw(ArgumentError(
                "strand centres must be distinct"
            ))
            logarithmic_sum += 2log(distance)
        end
    end
    return exp(logarithmic_sum / count^2)
end

"""
$(TYPEDSIGNATURES)

Calculate the GMR of an annular conductor.

```math
\\log GMR=\\log r_2-\\mu_r\\left[
\\frac{r_1^4}{(r_2^2-r_1^2)^2}\\log\\frac{r_2}{r_1}
-\\frac{3r_1^2-r_2^2}{4(r_2^2-r_1^2)}\\right].
```

`r_in` and `r_ex` are the inner radius ``r_1`` and outer radius ``r_2``
in meters, with ``0≤r_1≤r_2`` and ``r_2>0``. `mu_r` is the positive,
dimensionless relative permeability. The returned GMR is in meters.

# Notes

The expression reduces to the solid-conductor result at zero inner radius and
to the outer radius for an infinitesimally thin shell.
"""
function tubular_gmr(r_ex::Real, r_in::Real, mu_r::Real)
    rex, rin, permeability = promote(float(r_ex), float(r_in), float(mu_r))
    rin >= zero(rin) && rex > zero(rex) && rex >= rin || throw(DomainError(
        (rin, rex), "outer radius must be positive and not smaller than inner radius"
    ))
    permeability > zero(permeability) || throw(DomainError(
        permeability, "relative permeability must be positive"
    ))
    isapprox(rin, rex) && return rex
    iszero(rin) && return rex * exp(-permeability / 4)
    denominator = rex^2 - rin^2
    term1 = rin^4 / denominator^2 * log(rex / rin)
    term2 = (3 * rin^2 - rex^2) / (4 * denominator)
    return exp(log(rex) - permeability * (term1 - term2))
end

"""
$(TYPEDSIGNATURES)

Calculate the GMR of two heterogeneous conductor zones:

```math
GMR_{eq}=GMR_1^{\\beta^2}GMR_2^{(1-\\beta)^2}
          GMD^{2\\beta(1-\\beta)},\\qquad
\\beta=\\frac{A_1}{A_1+A_2}.
```

# Arguments

- `gmr1`: GMR of the accumulated conductor zone \\[m\\].
- `area1`: Cross-sectional area of the accumulated conductor zone \\[m²\\].
- `gmr2`: GMR of the added conductor zone \\[m\\].
- `area2`: Cross-sectional area of the added conductor zone \\[m²\\].
- `gmd`: Geometric mean distance between the zones \\[m\\].

# Returns

- Equivalent GMR of the combined conductor \\[m\\].
"""
function equivalent_gmr(
        gmr1::Real,
        area1::Real,
        gmr2::Real,
        area2::Real,
        gmd::Real
)
    first_gmr, first_area, second_gmr, second_area, distance = promote(
        float(gmr1), float(area1), float(gmr2), float(area2), float(gmd)
    )
    first_gmr > zero(first_gmr) && second_gmr > zero(second_gmr) ||
        throw(DomainError((first_gmr, second_gmr), "GMR values must be positive"))
    first_area > zero(first_area) && second_area > zero(second_area) ||
        throw(DomainError((first_area, second_area), "conductor areas must be positive"))
    distance > zero(distance) || throw(DomainError(
        distance, "geometric mean distance must be positive"
    ))
    fraction = first_area / (first_area + second_area)
    complement = one(fraction) - fraction
    return first_gmr^(fraction^2) * second_gmr^(complement^2) *
           distance^(2 * fraction * complement)
end

"""
$(TYPEDSIGNATURES)

Recover dimensionless relative permeability from the GMR of an annular conductor:

```math
\\mu_r=\\frac{\\log(r_{ex}/GMR)}{A-B},\\qquad
A=\\frac{r_{in}^4}{(r_{ex}^2-r_{in}^2)^2}\\log\\frac{r_{ex}}{r_{in}},\\qquad
B=\\frac{3r_{in}^2-r_{ex}^2}{4(r_{ex}^2-r_{in}^2)}.
```

`gmr`, `r_ex`, and `r_in` are the geometric mean radius, outer conductor
radius, and inner conductor radius in meters. The GMR and outer radius must
be positive, with ``0≤r_{in}≤r_{ex}``. A `NaN` input returns `NaN`.

For a solid conductor (`r_in=0`), the result is
``-4\\log(GMR/r_{ex})``. The implementation sets ``A=0`` when the radii are
approximately equal. In the thin-shell limit, GMR approaches the outer radius
regardless of permeability, so the inverse does not determine permeability
reliably.
"""
function equivalent_mu(gmr::Real, r_ex::Real, r_in::Real)
    radius, rex, rin = promote(float(gmr), float(r_ex), float(r_in))
    any(isnan, (radius, rex, rin)) && return radius + rex + rin
    radius > zero(radius) || throw(DomainError(radius, "GMR must be positive"))
    rin >= zero(rin) && rex > zero(rex) && rex >= rin || throw(DomainError(
        (rin, rex), "outer radius must be positive and not smaller than inner radius"
    ))
    is_solid = iszero(rin) || isapprox(rin, rex)
    denominator = rex^2 - rin^2
    term1 = is_solid ? zero(rin) : rin^4 / denominator^2 * log(rex / rin)
    term2 = (3 * rin^2 - rex^2) / (4 * denominator)
    return -(log(radius) - log(rex)) / (term1 - term2)
end
