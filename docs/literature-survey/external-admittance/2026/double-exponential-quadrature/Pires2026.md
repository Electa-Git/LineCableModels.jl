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
\quad x_k=\tanh\!\left[\frac\pi2\sinh(kh)\right],
\quad w_k=\frac{(\pi/2)h\cosh(kh)}{\cosh^2[(\pi/2)\sinh(kh)]}.
\qquad\text{(12–14)}
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
