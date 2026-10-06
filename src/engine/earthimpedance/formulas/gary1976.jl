function description(::Type{<:Formula{:gary1976}}; compact::Bool = false)
    compact ? "Gary" : "Gary complex-depth approximation (1976)"
end

"""
$(TYPEDSIGNATURES)

Evaluate Gary's complex-depth approximation for overhead self and mutual
exterior impedance, in Ω/m.

# Assumptions

Aerial conductors above homogeneous, nonmagnetic conductive earth, with an
image plane at a complex depth. Displacement current is neglected.

# Expression

For the positive-time convention ``e^{j\\omega t}``, define

```math
h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},\\qquad
\\mu_0=4\\pi\\,10^{-7}\\ \\mathrm{H/m}.
```

Soil conductivity ``\\sigma_g`` is in S/m. The principal square root gives
complex depth ``h_e`` in m. With aerial heights ``h_i,h_j`` and horizontal
separation ``x``, the direct and complex-image distances are

```math
d=\\sqrt{x^2+(h_i-h_j)^2},\\qquad
S=\\sqrt{x^2+(h_i+h_j+2h_e)^2}.
```

The complete exterior mutual impedance is

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{S}{d}.
```

For self impedance, the image distance is ``2(h_i+h_e)`` and the direct distance
is the exterior radius ``r_i``:

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{2(h_i+h_e)}{r_i}.
```

All geometric lengths are in m. The expression includes the ideal-earth exterior
term. Returning only the lossy-ground correction ``\\ln(S/D)`` with real-image
distance ``D`` would omit that term from the package's assembly.

# Reference

C. Gary, “Approche complète de la propagation multifilaire en haute fréquence
par utilisation des matrices complexes,” *EDF Bulletin de la Direction des
Études et Recherches*, série B (1976). The implementation follows the
complex-image expression reproduced in
[PSCAD's earth-return impedance documentation](https://www.pscad.com/webhelp-v5-ol/EMTDC/Transmission_Lines/Mutual_Impedance_with_Earth_Return.htm),
Eq. (8-28) and its complex-depth diagram. PSCAD calls its native implementation
Deri-Semlyen.
"""
function earth_impedance(
        ::Formula{:gary1976}, ::Val{:self}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    s = functor.state.jω
    μ0 = vacuum_permeability(typeof(real(s)))
    he = inv(sqrt(s * μ0 * functor.state.sigma[2]))
    return s * μ0 / (2 * (one(real(s)) * π)) * log(2 * (pair.heights[1] + he) / pair.radius)
end

function earth_impedance(
        ::Formula{:gary1976}, ::Val{:mutual}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    s = functor.state.jω
    μ0 = vacuum_permeability(typeof(real(s)))
    he = inv(sqrt(s * μ0 * functor.state.sigma[2]))
    hi, hj = pair.heights
    x = pair.separation
    S = sqrt((hi + hj + 2he)^2 + x^2)
    d = hypot(x, hi - hj)
    return s * μ0 / (2 * (one(real(s)) * π)) * log(S / d)
end

function formulation_options(::FormulaMethod{<:Formula{:gary1976}, typeof(earth_impedance)})
    FormulationOptions()
end

:gary1976
