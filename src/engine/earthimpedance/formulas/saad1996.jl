function assumptions(::Val{:saad1996})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

function description(::Type{<:Formula{:saad1996}}; compact::Bool = false)
    compact ? "Saad" : "Saad underground closed form (1996)"
end

"""
$(TYPEDSIGNATURES)

Evaluate Saad, Gaba, and Giroux's closed-form approximation for underground
self and mutual earth-return impedance, in Ω/m.

# Assumptions

Parallel horizontal cables in homogeneous, linear, isotropic, nonmagnetic
conductive earth below air. Displacement current is neglected, and the
wavelength is long compared with the transverse dimensions. Neither air
propagation nor a nonzero imposed longitudinal propagation constant is retained.
The earth-return calculation treats cables as filaments, substituting the
exterior cable radius in the self term. Internal conductor and insulation
contributions are separate.

# Expression

With the package's positive-time convention ``e^{j\\omega t}``, define

```math
m=\\sqrt{\\frac{j\\omega\\mu_0}{\\rho_g}},\\qquad
\\ell=h_i+h_j,\\qquad
d=\\sqrt{x^2+(h_i-h_j)^2},\\qquad
D=\\sqrt{x^2+(h_i+h_j)^2}.
```

Here ``\\rho_g`` is earth resistivity [Ω·m], ``\\mu_0=4\\pi\\,10^{-7}`` H/m,
and ``m`` has units of inverse meters. Depths ``h_i,h_j`` are positive downward.
``x`` is horizontal separation, ``d`` is axis distance, and ``D`` is image
distance, all in meters. The implementation uses the principal square root.
Its branch and a separate time factor remain unspecified in the source.

The proposed mutual and self expressions, Eqs. (5)-(6), repeated as (26)-(27), are

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[K_0(md)+\\frac{2e^{-m\\ell}}{4+m^2x^2}\\right],
```

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[K_0(mr_i)+\\frac{2e^{-2mh_i}}{4+m^2r_i^2}\\right].
```

``K_0`` is the modified Bessel function of the second kind, order zero.
``r_i`` is the exterior cable radius [m], written ``R`` in the paper.
The source prefactor ``\\rho_g m^2/(2\\pi)`` equals ``j\\omega\\mu_0/(2\\pi)``.
The self expression follows the mutual geometry with ``x=r_i`` and
``h_i=h_j``. The mutual correction denominator uses horizontal separation,
not axis distance.

# Approximation and limitations

The source starts from the Pollaczek/Wedepohl representations, Eqs. (1)-(4):

```math
Z_m=\\frac{\\rho_gm^2}{2\\pi}
\\left[K_0(md)-K_0(mD)+J_m\\right],\\qquad
J_m=\\int_{-\\infty}^{\\infty}
\\frac{e^{-\\ell\\sqrt{\\gamma^2+m^2}}}
{|\\gamma|+\\sqrt{\\gamma^2+m^2}}e^{j\\gamma x}\\,d\\gamma,
```

```math
Z_s=\\frac{\\rho_gm^2}{2\\pi}
\\left[K_0(mR)-K_0\\!\\left(m\\sqrt{R^2+4h^2}\\right)+J_s\\right],\\qquad
J_s=\\int_{-\\infty}^{\\infty}
\\frac{e^{-2h\\sqrt{\\gamma^2+m^2}}}
{|\\gamma|+\\sqrt{\\gamma^2+m^2}}e^{j\\gamma R}\\,d\\gamma.
```

Here ``\\gamma`` is the Fourier integration variable [1/m], not ``m`` or a
longitudinal line propagation constant. The source declares these parent
expressions without attributing new kernels to Saad et al. The paper cites
Pollaczek's 1931 French publication separately from the 1926 overhead result.

After contour deformation, Eq. (15) approximates
``\\sqrt{\\delta^2+1}/(\\delta+\\sqrt{\\delta^2+1})`` by
``(1+e^{-2\\delta})/2``. Equations (20)-(21) then use
``\\sqrt{\\delta^2+1}\\simeq1`` in the rapidly decaying part. The approximated
interface integral cancels the separate image Bessel term. The further
small-argument reduction discussed in the paper is not this implementation.

The contour proof is restricted to ``x/\\ell<1``. The paper reports about 3%
maximum relative error for the first kernel approximation, errors below about
1.5% for typical ``x/\\ell<1`` geometries, and negligible error through 10 kHz
for its ``x/\\ell=5`` examples. These reported test cases do not establish universal
bounds or runtime acceptance conditions.

# Reference

O. Saad, G. Gaba, and M. Giroux, “A Closed-Form Approximation for Ground Return
Impedance of Underground Cables,” *IEEE Transactions on Power Delivery*,
11(3), 1536–1545, 1996. Model and parent expressions: pp. 1536–1537,
Eqs. (1)-(4). Proposed equations: p. 1537, Eqs. (5)-(6). Derivation:
pp. 1537–1539, Eqs. (7)-(27). Error discussion: p. 1540 and Figs. 5–6.
The transcription was checked against the original publication's page images.
"""
function earth_impedance(
        ::Formula{:saad1996}, ::Val{:self}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    s = functor.state.jω
    μ0 = vacuum_permeability(typeof(real(s)))
    m = sqrt(s * μ0 / functor.state.rho[2])
    h, r = abs(pair.heights[1]), pair.radius
    return s * μ0 / (2 * (one(real(s)) * π)) *
           (special_besselk(0, m * r) + 2exp(-2m * h) / (4 + (m * r)^2))
end

function earth_impedance(
        ::Formula{:saad1996}, ::Val{:mutual}, ::Val{2}, ::Val{2},
        functor, pair, workspace
)
    s = functor.state.jω
    μ0 = vacuum_permeability(typeof(real(s)))
    m = sqrt(s * μ0 / functor.state.rho[2])
    hi, hj = abs.(pair.heights)
    x = pair.separation
    d = hypot(x, hi - hj)
    return s * μ0 / (2 * (one(real(s)) * π)) *
           (special_besselk(0, m * d) + 2exp(-m * (hi + hj)) / (4 + (m * x)^2))
end

function formulation_options(::FormulaMethod{<:Formula{:saad1996}, typeof(earth_impedance)})
    FormulationOptions()
end

:saad1996
