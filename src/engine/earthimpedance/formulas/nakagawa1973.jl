function routes(identifier::Val{:Nakagawa1973})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Nakagawa1973})
    (
        air = _lossless,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Nakagawa1973}) = Val(:zero)
media(::Formula{:Nakagawa1973}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel line currents above stratified earth; conductor skin, insulation, and a general finite-radius self term are excluded. |
| Calculated quantities | Mutual external and earth-return impedance, own-image self correction, and two-layer and homogeneous reductions. |
| Earth structure | Air above two finite earth layers and a lower earth half-space; ``d_2`` is cumulative depth. No arbitrary-layer recursion is supplied. |
| Model and approximation | Integral representation with the source's longitudinal prescription ``\\gamma_0=jk`` and Bessel-transform reduction. The two- and one-layer formulae are geometric limits; no approximation error bound is supplied. |
| Main source | M. Nakagawa, A. Ametani, and K. Iwamoto (1973). |
| Citation key(s) | Primary: `:Nakagawa1973`; later witness: `:Ametani1975` |
| Evidence status | All 1973 and 1975 pages checked. Appendix grouping, derivation defects, the root branch, and the finite-radius self prescription remain unresolved. |

**Expression.**

The source's final coefficient (14) and reductions (15) are evaluated
as a bounded reflection ratio ``r=c_2/c_1``:

```math
B_2=\\frac{1+r}{\\lambda(1+r)+\\mu_0b_1(1-r)}.
```

The two finite thicknesses are ``d_1`` and ``d_2-d_1``.
One to three earth layers are supported; no arbitrary-layer formula is
attributed to this paper. The positive-real transverse root supplies
decaying fields. The source's own-image correction uses zero lateral
separation for self. The standard thin-wire ideal-ground logarithm
completes the external self term for matrix assembly; it is not a
new finite-radius result attributed to Nakagawa.

**Reference.** M. Nakagawa, A. Ametani, and K. Iwamoto, “Further Studies on
Wave Propagation in Overhead Lines with Earth Return: Impedance of Stratified
Earth,” *Proceedings of the IEE*, 120, 1521–1528, 1973.
"""
function description(::Formula{:Nakagawa1973})
    "Nakagawa et al. three-layer overhead impedance and reductions (1973)"
end

function propagation_constant(
        ::Val{:Nakagawa1973}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Nakagawa1973})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    2<=length(rho)<=4 || throw(DimensionMismatch(
        ":Nakagawa1973 covers one to three earth layers, not an arbitrary-layer recursion"
    ))
    isinf(first(rho)) || throw(ArgumentError(":Nakagawa1973 requires the lossless-air longitudinal prescription"))
    return _stratified_functor(
        Val(:Nakagawa1973), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate the final Nakagawa coefficient and its homogeneous/two-layer
reductions using finite-layer thicknesses and bounded reflection ratios.
"""
function earth_impedance(
        ::Val{:Nakagawa1973}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    N=length(state.rho)-1
    self=pair.row==pair.column
    lateral=self ? zero(geometry.H) : geometry.y_ij
    ideal=self ? log(geometry.H/geometry.y_ij) : log(geometry.D_ij/geometry.d_ij)
    integral=_height_quadrature(state,geometry.H) do lambda
        a=ntuple(N) do m
            spectral_root(lambda^2+state.gamma_medium_squared[m+1]-
                state.gamma_medium_squared[1],state.jω)
        end
        b=ntuple(m->a[m]/state.mu[m+1],N)
        # Divide the printed c2/c1 by their common scale. The middle
        # thickness is d2-d1, not cumulative d2 or d1+d2.
        reflection=zero(state.jω)
        for m in (N-1):-1:1
            sumterm=b[m]+b[m+1]
            contrast=b[m]-b[m+1]
            reflection=(contrast+sumterm*reflection)/
                (sumterm+contrast*reflection)*exp(-2a[m]*state.thickness[m+1])
        end
        B=(1+reflection)/
            (lambda*(1+reflection)+state.mu[1]*b[1]*(1-reflection))
        B*exp(-lambda*geometry.H)*cos(lambda*lateral)
    end
    return _complex_result(state.jω,state.jω*state.mu[1]/
        (2*(one(geometry.H)*π))*(ideal+2integral))
end

:Nakagawa1973
