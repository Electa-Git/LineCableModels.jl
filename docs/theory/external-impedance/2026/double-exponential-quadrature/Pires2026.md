# Pires–Moreira double-exponential cable-integral evaluator

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Cable self/mutual entries evaluated independently. Geometry enters parent integrands, not the transform. |
| Calculated quantities | Numerical quadrature for cable earth-return impedance integrals |
| Earth structure | The quadrature is structure-agnostic; applications use the stated cable media. |
| Model and approximation | Finite implementation truncates the formally infinite node sum; the physical parent kernel is unchanged. |
| Main source | D. L. Pires, F. A. Moreira, M. G. Soares, and F. M. Vasconcellos (2026) |
| Citation key(s) | `:Pires2026` |
| Evidence status | Page-image verified |

**Description.** Double-exponential variable transformation used as a stable evaluator for semi-infinite cable-constant integrals.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Inherited from each parent cable integral. | Application sections. |
| Air propagation constant ``γ_air`` | Inherited. | Parent formulas. |
| Earth propagation constant ``γ_earth`` | Inherited. | Parent formulas. |
| Earth permittivity and displacement current | Not altered by the evaluator. | Method. |
| Range of validity | Integrals transformable to a finite interval with decaying endpoint behavior. | §2. |
| Earth permeability ``μ_earth`` | Inherited. | Parent formulas. |
| Arrangement | Cable self/mutual entries evaluated independently. | Applications. |
| Earth structure | The quadrature is structure-agnostic; applications use the stated cable media. | Method. |
| Conductor and insulation geometry | Geometry enters parent integrands, not the transform. | Method. |
| Constitutive and field assumptions | Numerical evaluator only. | Contribution statement. |
| Conventions | Symmetric integer nodes ``k`` and step ``h``. | (12)–(14). |

**Expression.**

```math
\int_{-1}^{1}f(x)\,dx\simeq\sum_{k=-\infty}^{\infty}w_kf(x_k),
\qquad\text{(12)}
```

```math
\begin{aligned}
x_k&=\tanh\!\left[\frac\pi2\sinh(kh)\right], &&\text{(13)}\\
w_k&=\frac{(\pi/2)h\cosh(kh)}{\cosh^2[(\pi/2)\sinh(kh)]}. &&\text{(14)}
\end{aligned}
```

**Approximation.** Finite implementation truncates the formally infinite node sum; the physical parent kernel is unchanged.

**Limitations.** This is an evaluator record, not a new earth constitutive model.

**Reference.** [Pires2026](@cite), equations (12)–(14) and impedance examples.

**Transcription source.** Accepted and final source PDFs cross-checked.

## Source transcription

The same transform also evaluates admittance/potential kernels; that taxonomic output is recorded separately.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``h`` | unchanged | quadrature mesh size | dimensionless |
| ``x_k,w_k`` | unchanged | transformed nodes and weights | dimensionless |
| ``f`` | unchanged | parent impedance integrand | inherited units |

## Evidence and approximation sources

The source states the formula in its numerical-integration method.

## Limitations and discrepancies

No universal truncation count or error bound is inferred beyond the paper's tests.

## Numerical interpretation and implementation

The evaluator is selected as `quadrature=:double_exponential` on
`Xue2018b`, without introducing another physical formula identifier.
The paper's equations (5) and (6) are algebraic rearrangements of the
same quasi-TEM buried-cable kernels. Both half-spaces may retain
conductivity and displacement current; the application to seawater
above seabed does not require a new integral.

The printed appendix map has denominator ``x-1`` and maps the interval
to negative spectral coordinates. For the stated positive
semi-infinite integral, the numerical interpretation uses

```math
\lambda=\frac{1+x}{1-x},\qquad
\frac{d\lambda}{dx}=\frac{2}{(1-x)^2}.
```

Composing this map with (13) before numerical evaluation gives

```math
\begin{aligned}
\lambda_k&=\exp[\pi\sinh(kh)],\\
W_k&=h\pi\cosh(kh)\lambda_k,\\
\int_0^\infty F(\lambda)\,d\lambda
&\simeq\sum_k W_kF(\lambda_k).
\end{aligned}
```

This form avoids subtracting rounded endpoints. The mesh is halved,
and two successive refinements must satisfy the requested tolerance
and a check on the retained tail terms. Endpoint sampling is bounded
by working precision. These are implementation stopping controls,
not a universal error bound or a claim to reproduce the source's
flowchart step for step. Failure to establish convergence is reported.

The implementation is checked against the separate source impedance
and potential integrals, including lossy seawater and seabed, and
against assembled cable matrices. Potential coefficients are assembled
with the radial insulation terms before computing
``Y=j\omega P^{-1}``.
