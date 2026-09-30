function assumptions(::Val{:wedepohl1973})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

function description(::Type{<:Formula{:wedepohl1973}}; compact::Bool = false)
    compact ? "Wedepohl" : "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

"""
$(TYPEDSIGNATURES)

Evaluate the Wedepohl–Wilcox low-order approximation for underground
self and mutual exterior impedance, in Ω/m.

# Assumptions

Buried conductors in homogeneous conductive earth, neglecting displacement
current. The expansion assumes small penetration-depth products and uses the
supplied absolute soil permeability. It is distinct from the internal-impedance
formulation with the same `:wedepohl1973` identifier.

# Expression

For the positive-time convention ``e^{j\\omega t}``, define

```math
m=\\sqrt{j\\omega\\mu_g/\\rho_g},\\qquad e_c=1.7811.
```

Here ``\\rho_g`` is resistivity [Ω·m], ``\\mu_g`` is permeability [H/m], and
``m`` is inverse penetration depth [1/m], using the principal square root.
The published rounded logarithmic constant ``e_c`` approximates the exponential
of Euler's constant.

For positive burial depths ``h_i,h_j``, horizontal separation ``x``, exterior
radius ``r_i``, and axis distance ``d=\\sqrt{x^2+(h_i-h_j)^2}`` (all lengths in m),

```math
Z_{e,ii}=\\frac{j\\omega\\mu_g}{2\\pi}
\\left[-\\ln\\left(\\frac{e_cmr_i}{2}\\right)+\\frac12-\\frac43mh_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_g}{2\\pi}
\\left[-\\ln\\left(\\frac{e_cmd}{2}\\right)+\\frac12-\\frac23m(h_i+h_j)\\right].
```

The implementation evaluates this approximation at every requested frequency. Coincident mutual axes are invalid
physical geometry; a vertical pair at distinct depths has ``d>0`` and needs
no special case.

# Reference

L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of Underground
Power-Transmission Systems: System-Model and Wave-Propagation Characteristics,”
*Proceedings of the IEE* 120, 253–260 (1973), p. 255, Eqs. (7)–(8),
[doi:10.1049/piee.1973.0056](https://doi.org/10.1049/piee.1973.0056).
The equations and rounded constant are also reproduced in
[PSCAD's earth-return impedance documentation](https://www.pscad.com/webhelp-v5-ol/EMTDC/Transmission_Lines/Mutual_Impedance_with_Earth_Return.htm),
Eqs. (8-31)–(8-32).
"""
function earth_impedance(
        ::Formula{:wedepohl1973}, ::Val{:self}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    s, μ = functor.state.jω, functor.state.mu[2]
    m = sqrt(s * μ / functor.state.rho[2])
    ec = one(real(m)) * 17811 / 10000
    return s * μ / (2 * (one(real(s)) * π)) *
           (-log(ec * m * pair.radius / 2) + one(m) / 2 - 4m * abs(pair.heights[1]) / 3)
end

function earth_impedance(
        ::Formula{:wedepohl1973}, ::Val{:mutual}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    s, μ = functor.state.jω, functor.state.mu[2]
    m = sqrt(s * μ / functor.state.rho[2])
    ec = one(real(m)) * 17811 / 10000
    hi, hj = abs.(pair.heights)
    d = hypot(pair.separation, hi - hj)
    return s * μ / (2 * (one(real(s)) * π)) *
           (-log(ec * m * d / 2) + one(m) / 2 - 2m * (hi + hj) / 3)
end

function formulation_options(::FormulaMethod{
        <:Formula{:wedepohl1973}, typeof(earth_impedance)})
    FormulationOptions()
end

:wedepohl1973
