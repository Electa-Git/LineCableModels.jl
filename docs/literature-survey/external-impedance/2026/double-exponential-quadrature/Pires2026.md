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
\quad x_k=\tanh\!\left[\frac\pi2\sinh(kh)\right],
\quad w_k=\frac{(\pi/2)h\cosh(kh)}{\cosh^2[(\pi/2)\sinh(kh)]}.
\qquad\text{(12–14)}
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
