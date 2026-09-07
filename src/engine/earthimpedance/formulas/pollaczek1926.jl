function routes(identifier::Val{:Pollaczek1926})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        overhead = FormulaMethod(identifier, earth_impedance, Val(:overhead)),
        underground = FormulaMethod(
            identifier, earth_impedance, Val(:underground)
        ),
        mixed = FormulaMethod(identifier, earth_impedance, Val(:mixed)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Pollaczek1926})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pollaczek1926}) = Val(:zero)

function Formula(::Val{:Pollaczek1926};evaluation::Symbol=:spectral,
        dcim_samples::Int=100,dcim_terms::Int=14,kwargs...)
    evaluation in (:spectral,:finite_integral,:bessel_product,:single_bessel,
        :hypergeometric,:recursive,:dcim,:dcim_two_level) ||
        throw(ArgumentError("unknown Pollaczek integral evaluation"))
    identifier=Val(:Pollaczek1926); defaults=routes(identifier)
    values=assumptions(identifier)
    if evaluation in (:dcim,:dcim_two_level)
        1<=dcim_terms<=dcim_samples÷2 && dcim_samples>=4 ||
            throw(ArgumentError("DCIM needs at least four samples and 1 <= terms <= samples/2"))
        levels=evaluation===:dcim ? 1 : 2
        values=merge(values,(dcim_levels=levels,dcim_samples,dcim_terms))
        defaults=merge(defaults,(;
            overhead=FormulaMethod(identifier,earth_impedance,Val(:dcim),Val(:overhead)),
            underground=FormulaMethod(identifier,earth_impedance,Val(:dcim),Val(:underground)),
            mixed=FormulaMethod(identifier,earth_impedance,Val(:dcim),Val(:mixed))))
    elseif evaluation !== :spectral
        defaults=merge(defaults,(underground=FormulaMethod(
            identifier,earth_impedance,Val(:series),Val(evaluation)),))
    end
    overrides=(;kwargs...)
    isempty(setdiff(keys(overrides),keys(defaults))) || throw(ArgumentError(
        "unknown routes for earth-impedance formula :Pollaczek1926"))
    return Formula(identifier,merge(defaults,overrides),values)
end
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinitely long straight parallel filamentary conductors. A physical round wire of radius ``\\rho`` is treated only through the stated small-radius self averaging; a buried wire is galvanically and dielectrically isolated from earth. No finite insulation impedance/admittance is supplied. |
| Calculated quantities | Complex generalized mutual-induction coefficient for air–air, earth–earth, and both mixed source/observation placements; source-prescribed finite-radius self coefficient and series-impedance assembly |
| Earth structure | Plane homogeneous conducting half-space (medium 1) below homogeneous air (medium 2). |
| Model and approximation | Not an analytical approximation within the explicitly reduced filamentary, two-dimensional physical model. The reductions ``\\Gamma=0``, ``\\varepsilon_1=0``, ``\\mu_1=\\mu_2=1`` and ``k_2=0`` precede application. The finite-radius self rule is an additional small-radius approximation; separate small- and large-``\\|k\\eta\\|`` asymptotic self expressions (59a)–(59b) are not substituted here. |
| Main source | F. Pollaczek (1926) |
| Citation key(s) | `:Pollaczek1926` |
| Evidence status | Original-publication scan image verified for the displayed equations and definitions |

**Numerical scope.** The implementation evaluates the conductive-earth specialization with the air spectral constant set to zero. It does not solve the retained-propagation parent.

**Expression.** The underground and mixed terms are

```math
Z_{e,ij}^{11}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-H\\sqrt{\\lambda^2+\\gamma_1^2}}}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_1^2}}
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
Z_{e,ij}^{01}=\\frac{j\\omega\\mu_0}{\\pi}\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda|h_i|-a_1|h_j|}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda,
\\qquad a_1=\\sqrt{\\lambda^2+\\gamma_1^2}.
```

**Reference.** F. Pollaczek, “Über das Feld einer unendlich langen
wechselstromdurchflossenen Einfachleitung,” *Elektrische Nachrichtentechnik*,
3, 339–360, 1926.

## Theodoulidis alternative integral evaluations

## Bessel-product series

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel buried conductors at depths ``h_1,h_2`` and horizontal separation ``x``; the source uses a filamentary mutual geometry and the stated radius substitution for self. No insulation region appears. |
| Calculated quantities | Exact modified-Bessel series for the Pollaczek interface integral, used in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous earth half-space below air. |
| Model and approximation | Not an analytical approximation of ``J_{\\mathrm{Pollaczek}}``: the infinite series is derived exactly from (2). Finite truncation is a numerical approximation; the author's rule of thumb is roughly ``20x/H`` terms for ``x>H``. The physical parent remains the TEM Pollaczek model, with optional Sunde constitutive replacement. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked |

## Single-Bessel series

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary mutual geometry at depths ``h_1,h_2`` and separation ``x``; source-prescribed radius substitution for self; no insulation term. |
| Calculated quantities | Exact single-modified-Bessel series for the Pollaczek interface integral in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous half-space. |
| Model and approximation | Not an analytical approximation to the parent integral: the infinite series follows from the exact Bessel identity used in (14)–(16). Any finite truncation is numerical. The source supplies no universal fixed term count for (5); it says convergence is poor for large separation and good for ``x<H``. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked |

## Hypergeometric series

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Depths ``h_1,h_2``, horizontal separation ``x``; radius used for the self ``x``; no insulation region. |
| Calculated quantities | Exact confluent-hypergeometric series for the Pollaczek interface integral in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous earth half-space. |
| Model and approximation | The infinite series is an exact representation of (2). It comes from the exact finite-range form (22) and Taylor expansions about ``t=1``; only finite truncation introduces approximation. The source's remainder estimate gives the stated 10-, 21-, and 40-term guarantees over ``t\\in[0,1]``. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked for final series and parent |

The physical formula retains the identifier `Pollaczek1926`.
Select `evaluation=:bessel_product`, `:single_bessel`, or
`:hypergeometric` for the three alternative evaluations.
`:finite_integral` selects the common finite parent (22).
The default remains `:spectral`; overhead and mixed routes are unchanged.

The numerical single-Bessel sum includes `n=0`, required by (9a)
and its nonzero limit at zero lateral spacing. Its denominator is
``2^n n!(2n-1)``, as printed, not the earlier corpus transcription
``2^{2n}n!(2n-1)``. The hypergeometric sum includes the common
``e^{-kH}`` factor obtained from (22). The printed (5),(6),(16)
are retained as source witnesses; these computational selections agree
with the independently integrated parent and the first Bessel series.

Working-precision arithmetic protects cancellation in (3).
The finite integral removes the endpoint square root through `t=cos(θ)`.
It also supplies the slow-convergence fallback. For large complex
arguments, the existing spectral parent is used. These changes concern
numerical evaluation, not displacement current, new mutual physics,
or a full-wave correction.

**Reference.** [Theodoulidis2012](@cite), equations (1)–(6), (9a),
(16), (22), and (28), pp. 807–810.

## Recursive finite-integral evaluation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent; mutual uses cable depths ``h_1,h_2`` and horizontal spacing ``x``; self substitutes conductor radius ``R`` for ``x`` and ``h_1=h_2``. |
| Calculated quantities | Recursive infinite-series evaluator for the finite integral in the Wedepohl–Wilcox decomposition of buried-cable earth-return impedance |
| Earth structure | Homogeneous soil half-space below air. |
| Model and approximation | The exponential in the finite integral is expanded in its everywhere-convergent power series, accumulated recursively, and truncated by tolerance. The source recommends the series only for ``\\|D/p\\|\\le30`` and approximates ``I_w=0`` beyond that threshold. |
| Main source | R. Iracheta-Cortez (2015) for the recursion; Pollaczek (1926) and Wedepohl–Wilcox (1973) for the physical kernel and decomposition |
| Citation key(s) | `:Iracheta2015` |
| Evidence status | Original publication page images checked |

Select `evaluation=:recursive` for the moment-series evaluator.
The moment definitions (5a)–(5b) fix the initialization and parity
recurrences by integration by parts:

```math
A_n=\\frac{a^{n-1}(1-a^2)^{3/2}+(n-1)A_{n-2}}{n+2},
\\qquad
B_n=\\frac{a^{n-1}\\sqrt{1-a^2}+(n-1)B_{n-2}}{n}.
```

Here ``a=H/D``, ``A_0=(\\theta-a\\sqrt{1-a^2})/2``,
``A_1=(1-a^2)^{3/2}/3``, ``B_0=\\theta``,
``B_1=\\sqrt{1-a^2}``, and ``\\theta=\\arccos(a)``.
This resolves the duplicated initializer label and the inconsistent
printed parity recurrences without altering the source-defined integrals.

The resulting residual is inserted into the checked finite decomposition
(22) of Theodoulidis. This preserves the Pollaczek normalization rather
than copying inconsistent prefactors from Iracheta's restatement (3c).
The source's prescription to omit the residual for ``|D/p|>30``
is retained as a distinct approximation. It has no universal error bound.

**Reference.** [Iracheta2015](@cite), equations (4)–(7) and the
cutoff discussion on p. 37; [Theodoulidis2012](@cite), equation (22).

## Rallis discrete complex images

## Overhead Carson correction

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent at heights ``h_1,h_2`` and horizontal separation ``x``. |
| Calculated quantities | Finite rational complex-image sum for the homogeneous-earth overhead Carson correction |
| Earth structure | Homogeneous conductive half-space. |
| Model and approximation | The denominator is fitted by complex exponentials with GPOF and each is integrated using the elementary Laplace–cosine identity. The exact Struve/Bessel form (4.18), separately reproduced by the thesis, is not attributed to this DCIM method and is already represented elsewhere in the corpus. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF evaluator; Carson for the parent kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

## Mixed Pollaczek correction

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filaments at overhead height ``h_1``, burial depth ``h_2`` and horizontal spacing ``x``. |
| Calculated quantities | Finite rational complex-image sum for mutual earth-return impedance between one overhead and one buried conductor |
| Earth structure | Homogeneous conductive half-space below air. |
| Model and approximation | GPOF fits the spectral factor over a finite sampling interval and (4.14), ``\\int_0^\\infty e^{-p\\lambda}\\cos(q\\lambda)d\\lambda=p/(p^2+q^2)``, integrates each term. The two-level sampling variant is an empirical refinement. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF representation; Pollaczek for the parent kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

## Underground Pollaczek correction

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent; depths ``h_1,h_2``, horizontal separation ``x``. |
| Calculated quantities | Finite discrete-complex-image sum for the homogeneous-earth buried/buried Pollaczek correction |
| Earth structure | Homogeneous conductive half-space below air. |
| Model and approximation | Complex poles ``s_n`` and residues ``c_n`` are fitted with GPOF after mapping ``\\gamma=P_at+P_b`` over ``t\\in[0,1]``; (4.8) then integrates each exponential analytically. This is a data-fitted finite-sum approximation, not an identity. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF representation; Pollaczek for the physical kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

## Numerical interpretation

Select `EarthImpedance.Formula(:Pollaczek1926;evaluation=:dcim)` for
one-level GPOF, or `evaluation=:dcim_two_level` for two-level fitting.
These selections change all three placement evaluators, not the physical
Pollaczek assumptions. The default spectral evaluator is unchanged.

The fit uses 100 uniformly spaced samples and at most 14 retained singular
directions per level. Appendix A defines the shifted sample matrices
``Y_1,Y_2``; for ``Y_1=UDV^H``, the image poles follow from the
eigenvalues of ``D^{-1}U^HY_2V``. A least-squares solve determines the
residues. Numerically unresolved directions and growing exponentials are
excluded before the residue solve. `dcim_samples` and `dcim_terms`
expose the sample and rank limits.

The sampling variable is divided by ``|k|`` before fitting. The
two-level implementation first fits the far interval, then fits its
residual on the near interval. Each retained residue is stored with its
sampling origin; this avoids forming overflowing coefficients before
multiplication by a decaying Bessel function. Coefficients are generated
in double precision, including for higher-precision output arithmetic.
Neither higher output precision nor the sample count constitutes an
error bound for the finite fit.

The fitted spectral function is integrated to infinity as prescribed by
the image sum. No parent-integral fallback is hidden behind this selection.
Tests separately verify the matrix-pencil recovery, integration of the
fitted function, parent comparisons, and matrix assembly. In particular,
large lateral separation can expose substantial buried-conductor fit
error; the spectral or exact-series selections remain available.

The one-level interval is ``0\\le\\lambda/|k|\\le28``, following
the construction of §4.3.2 used by §4.4. Two-level intervals are
``[5,100]`` and ``[0,5]``. The conductor radius enters only
the perfect-ground self logarithm; the self correction sets ``x=0``.

The one-level interval is ``0\\le\\lambda/|k|\\le28``; the two-level
intervals are ``[5,100]`` and ``[0,5]``. Burial depth remains
inside the fitted spectral numerator, so mixed coefficients are
regenerated for that depth. Exchanging source and receiver preserves
the same height and depth assignments.

The one-level path is
``\\gamma/|k|=\\kappa+t(\\sqrt{100+\\kappa^2}-\\kappa)``,
where ``\\kappa=k/|k|``. Two-level endpoints correspond to
``T_{02}/|k|=0.22`` and ``T_{01}/|k|=10``.
The affine spans are the differences between successive endpoints:
``\\sqrt{100+\\kappa^2}-\\sqrt{0.22^2+\\kappa^2}`` and
``\\sqrt{0.22^2+\\kappa^2}-\\kappa``.
These differences implement the connected paths in Fig. 4.4. The
unnumbered printed spans on p. 72 use a plus sign in the first path
and omit the subtraction of ``k`` in the second; they would not
connect the stated endpoints.

For ``H_n=H-s_n``, the branch of
``\\sqrt{H_n^2+x^2}`` is continuous from positive ``H_n``.
The shifted residue and scaled ``K_1`` are multiplied before their
exponentials are expanded. The original ``K_0(kd)-K_0(kD)``
term is unchanged.

**Reference.** [Rallis2013](@cite), Chapter 4 and Appendix A.

"""
function description(::Formula{:Pollaczek1926})
    "Pollaczek homogeneous-earth overhead, underground, and mixed impedance (1926)"
end

function propagation_constant(
        ::Val{:Pollaczek1926}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Pollaczek1926})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    base=_homogeneous_functor(
        Val(:Pollaczek1926), formula, rho, epsilon, mu, jω, Γ, segments
    )
    return hasproperty(formula.assumptions,:dcim_levels) ?
        _rallis_functor(formula,base) : base
end

"""
$(TYPEDSIGNATURES)

Select Pollaczek's leaf impedance from the physical conductor placement.

The registered recipe contains the homogeneous-earth overhead, underground,
and overhead-underground interactions. Its public identity therefore does not
encode which leaf the pair requires. The `overhead`, `underground`, and
`mixed` routes remain individually replaceable when composing an experiment.
"""
function earth_impedance(
        ::Val{:Pollaczek1926}, ::Val{:mutual}, functor, pair
)
    placement = _placement(pair)
    typeof(placement) === Val{:overhead} &&
        return functor.routes.overhead(functor, pair)
    typeof(placement) === Val{:underground} &&
        return functor.routes.underground(functor, pair)
    return functor.routes.mixed(functor, pair)
end

function earth_impedance(
        ::Val{:Pollaczek1926}, ::Val{:overhead}, functor, pair
)
    return earth_impedance(Val(:Carson1926), Val(:mutual), functor, pair)
end

raw"""
Evaluate Pollaczek's homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[K_0(\gamma_1d_{ij})-
K_0(\gamma_1D_{ij})+2\int_0^\infty
\frac{e^{-(h_i+h_j)\sqrt{\lambda^2+\gamma_1^2}}}
{\lambda+\sqrt{\lambda^2+\gamma_1^2}}
\cos(y_{ij}\lambda)\,d\lambda\right].
```
"""
function earth_impedance(
        ::Val{:Pollaczek1926}, ::Val{:underground}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    integral = _quadrature(state) do lambda
        u_1 = sqrt(lambda^2 + gamma^2)
        exp(-geometry.H * u_1) * cos(geometry.y_ij * lambda) /
        (lambda + u_1)
    end
    direct = _complex_result(
        state.jω,
        special_besselk(0, gamma * geometry.d_ij) -
        special_besselk(0, gamma * geometry.D_ij)
    )
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) * (direct + 2 * integral)
end

raw"""
Evaluate Pollaczek's overhead-underground mutual impedance:

```math
Z_{e,ij}^{01}=\frac{j\omega\mu_0}{\pi}\int_0^\infty
\frac{\mu_1e^{-\lambda|h_i|-a_1|h_j|}}
{\lambda\mu_1+a_1\mu_0}\cos(y_{ij}\lambda)d\lambda,
\qquad a_1=\sqrt{\lambda^2+\gamma_1^2}.
```
"""
function earth_impedance(
        ::Val{:Pollaczek1926}, ::Val{:mixed}, functor, pair
)
    state = functor.state
    air = pair.layers[1] == 1 ? 1 : 2
    earth = air == 1 ? 2 : 1
    h_air = abs(pair.heights[air])
    h_earth = abs(pair.heights[earth])
    integral = _quadrature(state) do lambda
        a_1 = sqrt(lambda^2 + state.gamma_medium_squared[2])
        state.mu[2] * exp(-lambda * h_air - a_1 * h_earth) /
        (lambda * state.mu[2] + a_1 * state.mu[1]) *
        cos(pair.separation * lambda)
    end
    πT = one(h_air) * π
    return state.jω * state.mu[1] / πT * integral
end


function earth_impedance(identifier::Val{:Pollaczek1926},::Val{:series},method,functor,pair)
    _require(pair,Val(:underground))
    value=_pollaczek_series_coefficient(functor.state,pair,method)
    return isnothing(value) ? earth_impedance(identifier,Val(:underground),functor,pair) : value
end
:Pollaczek1926
