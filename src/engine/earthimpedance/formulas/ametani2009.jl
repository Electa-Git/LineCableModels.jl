function assumptions(::Val{:ametani2009})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth approximation for mixed overhead-underground pairs.

**Availability.** Registered scientific identity; the coaxial implementation is
not yet implemented. No numerical fallback is provided.

**Expression.** Its distinctive mixed term is

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{2\\pi}e^{-h_g/h_e}\\ln\\frac{S}{D},
\\quad h_e=(j\\omega\\mu_0\\sigma_g)^{-1/2},
```

```math
S=\\sqrt{(h_a+h_g+2h_e)^2+y_{ij}^2},\\qquad
D=\\sqrt{(h_a+h_g)^2+y_{ij}^2}.
```

**Reference.** A. Ametani, T. Yoneda, Y. Baba, and N. Nagaoka, “An
Investigation of Earth-Return Impedance Between Overhead and Underground
Conductors and Its Approximation,” *IEEE Transactions on Electromagnetic
Compatibility*, 51, 860–867, 2009.
DOI: 10.1109/TEMC.2009.2019953.
PSCAD's help lists this journal article with a 2005 date; its journal volume
and DOI identify the 2009 publication used by this registration.
"""
function description(::Type{<:Formula{:ametani2009}}; compact::Bool = false)
    compact ? "Ametani" : "Ametani mixed-pair homogeneous-earth impedance (2009)"
end

function earth_impedance(
        ::Formula{:ametani2009}, kind::Val{:mutual}, ::Val{1}, ::Val{2},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :ametani2009 ($kind), source layer 1, target layer 2: not yet implemented for the coaxial backend"))
end

function earth_impedance(
        ::Formula{:ametani2009}, kind::Val{:mutual}, ::Val{2}, ::Val{1},
        functor, pair, workspace
)
    throw(ArgumentError("earth_impedance :ametani2009 ($kind), source layer 2, target layer 1: not yet implemented for the coaxial backend"))
end

formulation_options(::FormulaMethod{<:Formula{:ametani2009}, typeof(earth_impedance)}) =
    FormulationOptions()

:ametani2009
