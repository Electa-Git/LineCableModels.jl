@inline function _mu0(value)
    return one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)
end

"""
$(TYPEDSIGNATURES)

Calculate internal tubular inductance in the DC approximation:
``L=\\mu_r\\mu_0\\log(r_{ex}/r_{in})/(2\\pi)`` \\[H/m\\].
"""
function tubular_inductance(r_in::Real, r_ex::Real, mu_r::Real)
    rin, rex, permeability = promote(float(r_in), float(r_ex), float(mu_r))
    rin > zero(rin) && rex > rin || throw(DomainError(
        (rin, rex), "tubular inductance requires 0 < r_in < r_ex"
    ))
    permeability >= zero(permeability) || throw(DomainError(
        permeability, "relative permeability must be nonnegative"
    ))
    return permeability * _mu0(rin) * log(rex / rin) / (2 * (one(rin) * π))
end

"""
$(TYPEDSIGNATURES)

Calculate the GMR of a circular wire array:
``GMR=(r_g N a^{N-1})^{1/N}``, with
``r_g=r_w\\exp(-\\mu_r/4)``.
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

Calculate the GMR of identical circular strands from their resolved centres:

```math
\\log GMR=\\frac{1}{N^2}\\left[
N\\log r_g+2\\sum_{i=1}^{N}\\sum_{j=i+1}^{N}\\log d_{ij}
\\right],\\qquad
r_g=r_w\\exp(-\\mu_r/4).
```

# Arguments

- `coordinates`: Strand-centre coordinates \\[m\\].
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

Recover relative permeability by inverting [`tubular_gmr`](@ref).

# Notes

The solid-conductor and thin-shell limits use their analytic reductions.
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

"""
$(TYPEDSIGNATURES)

Calculate the positive-sequence inductance of a solid-bonded trefoil cable
from the CIGRE TB-531 impedance reduction
``Z_d=(Z_a-Z_x)-(Z_m-Z_x)^2/(Z_s-Z_x)``.
"""
function trefoil_inductance(
        r_in_co::Real,
        r_ex_co::Real,
        rho_co::Real,
        mu_co::Real,
        r_in_screen::Real,
        r_ex_screen::Real,
        rho_screen::Real,
        mu_screen::Real,
        spacing::Real;
        rho_e::Real = oftype(float(spacing), 100),
        f::Real = oftype(float(spacing), 50)
)
    values = promote(
        float(r_in_co), float(r_ex_co), float(rho_co), float(mu_co),
        float(r_in_screen), float(r_ex_screen), float(rho_screen),
        float(mu_screen), float(spacing), float(rho_e), float(f)
    )
    rin, rex, rhoc, muc, rins, rexs, rhos, mus, distance, earth, frequency = values
    omega = oftype(frequency, 2π) * frequency
    mu0 = _mu0(frequency)
    coefficient = mu0 / (2 * (one(frequency) * π))
    earth_depth = oftype(frequency, 659) * sqrt(earth / frequency)
    earth_resistance = omega * mu0 / 8
    core_gmr = tubular_gmr(rex, rin, muc)
    screen_gmr = tubular_gmr(rexs, rins, mus)
    Za = earth_resistance + rhoc / (oftype(rhoc, π) * (rex^2 - rin^2)) +
         im * omega * coefficient * log(earth_depth / core_gmr)
    Zs = earth_resistance + rhos / (oftype(rhos, π) * (rexs^2 - rins^2)) +
         im * omega * coefficient * log(earth_depth / screen_gmr)
    Zm = earth_resistance + im * omega * coefficient * log(earth_depth / screen_gmr)
    Zx = earth_resistance + im * omega * coefficient * log(earth_depth / distance)
    return imag((Za - Zx) - (Zm - Zx)^2 / (Zs - Zx)) / omega
end
