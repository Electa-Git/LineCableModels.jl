function assumptions(::Val{:lucca1994})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

function description(::Type{<:Formula{:lucca1994}}; compact::Bool = false)
    compact ? "Lucca" : "Lucca mixed-pair homogeneous-earth impedance (1994)"
end

"""
$(TYPEDSIGNATURES)

Evaluate Lucca's mutual exterior impedance between an aerial and a buried
conductor, in Ω/m.

# Assumptions

Homogeneous, nonmagnetic conductive earth, neglecting displacement current.
Both ordered aerial/buried interactions use the same reciprocal expression.
The formula supplies no self or same-half-space interactions; those remain
separate selections.

# Expression

With the positive-time convention ``e^{j\\omega t}``, let ``h_a>0`` be aerial
height, ``h_g>0`` burial depth, and ``x`` horizontal separation, all in m.
Signed package heights are converted to these positive physical distances.
With soil conductivity ``\\sigma_g`` [S/m], define

```math
h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},\\qquad
H=h_a+h_g+2h_e,\\qquad
S=\\sqrt{H^2+x^2},\\qquad
D=\\sqrt{(h_a+h_g)^2+x^2},
```

where ``\\mu_0=4\\pi\\,10^{-7}`` H/m and square roots use their principal value.
The complete mixed coefficient is

```math
Z_{e,ag}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[\\ln\\frac{S}{D}
-\\frac23\\left(\\frac{h_e}{S^2}\\right)^3 H(H^2-3x^2)\\right].
```

The correction uses ``(h_e/S^2)^3``, not ``h_e^3/S^3``.

# Reference

G. Lucca, “Mutual Impedance Between an Overhead and a Buried Line with Earth
Return,” *9th International Conference on Electromagnetic Compatibility* (1994),
[doi:10.1049/cp:19940679](https://doi.org/10.1049/cp:19940679).
The implemented expression is reproduced in
[PSCAD's earth-return impedance documentation](https://www.pscad.com/webhelp-v5-ol/EMTDC/Transmission_Lines/Mutual_Impedance_with_Earth_Return.htm),
Eq. (8-34), with the geometric definitions immediately below it.
"""
function earth_impedance(
        ::Formula{:lucca1994}, ::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace
)
    s = functor.state.jω
    μ0 = 4 * (one(real(s)) * π) * (one(real(s)) * 10)^(-7)
    he = inv(sqrt(s * μ0 * functor.state.sigma[2]))
    vertical = abs(pair.heights[1]) + abs(pair.heights[2])
    x = pair.separation
    H = vertical + 2he
    S2 = H^2 + x^2
    D = hypot(vertical, x)
    return s * μ0 / (2 * (one(real(s)) * π)) *
           (log(sqrt(S2) / D) - 2 * (he / S2)^3 * H * (H^2 - 3x^2) / 3)
end

# The reciprocal mixed coefficient has the same physical distances in either direction.
function earth_impedance(
        selected::Formula{:lucca1994}, kind::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace
)
    return earth_impedance(selected, kind, Val(1), Val(2), functor, pair, workspace)
end

function formulation_options(::FormulaMethod{
        <:Formula{:lucca1994}, typeof(earth_impedance)})
    FormulationOptions()
end

:lucca1994
