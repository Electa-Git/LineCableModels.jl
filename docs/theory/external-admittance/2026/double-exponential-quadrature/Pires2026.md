# Pires–Moreira double-exponential potential/admittance evaluator

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Self/mutual matrix entries evaluated individually. Geometry enters parent integrand. |
| Calculated quantities | Numerical quadrature for earth potential/admittance cable integrals |
| Earth structure | Structure enters parent integrand. |
| Model and approximation | Finite node truncation only; no change to the selected admittance physics. |
| Main source | D. L. Pires, F. A. Moreira, M. G. Soares, and F. M. Vasconcellos (2026) |
| Citation key(s) | `:Pires2026` |
| Evidence status | Page-image verified |

**Description.** Double-exponential quadrature applied to the potential/admittance integrals in cable-constants programs.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Inherited from parent potential kernel. | Applications. |
| Air propagation constant ``γ_air`` | Inherited. | Parent kernel. |
| Earth propagation constant ``γ_earth`` | Inherited. | Parent kernel. |
| Earth permittivity and displacement current | Unchanged by numerical transform. | Method. |
| Range of validity | Suitable endpoint-decaying transformed integrals. | §2. |
| Earth permeability ``μ_earth`` | Inherited. | Parent kernel. |
| Arrangement | Self/mutual matrix entries evaluated individually. | Application. |
| Earth structure | Structure enters parent integrand. | Method. |
| Conductor and insulation geometry | Geometry enters parent integrand. | Method. |
| Constitutive and field assumptions | Pure numerical evaluator. | Contribution. |
| Conventions | Symmetric nodes ``k`` and mesh ``h``. | (12)–(14). |

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

**Approximation.** Finite node truncation only; no change to the selected admittance physics.

**Limitations.** Evaluator rather than a new material/interface kernel.

**Reference.** [Pires2026](@cite), equations (12)–(14) and admittance examples.

**Transcription source.** Accepted and final PDF versions cross-checked.

## Source transcription

Separate series-impedance and shunt-admittance records prevent application labels from merging the outputs.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``h`` | unchanged | quadrature mesh | dimensionless |
| ``x_k,w_k`` | unchanged | DE nodes and weights | dimensionless |
| ``f`` | unchanged | parent potential integrand | inherited units |

## Evidence and approximation sources

The same evaluator applies to the impedance and admittance integrals stated by the source.

## Limitations and discrepancies

No universal stopping rule is supplied by the source.

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
