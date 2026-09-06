# Maaouni–Amri two-half-space potential/admittance approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Scalar self/mutual wire coefficients. Thin bare wire; no coating. |
| Calculated quantities | Thin-wire potential coefficient and qTEM admittance near a lossy interface |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | qTEM/asymptotic analytical approximation; not a full-wave solved ``Γ``. |
| Main source | A. Maaouni, A. Amri, and N. Zouhir (2001) |
| Citation key(s) | `:Maaouni2001` |
| Evidence status | Original publication page images checked |

**Description.** Potential/admittance companion of the source's two-half-space qTEM wire formulation.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | qTEM approximation of exact spectral parent. | §3. |
| Air propagation constant ``γ_air`` | ``k_0``. | Definitions. |
| Earth propagation constant ``γ_earth`` | Complex medium ratio ``n``. | Parent equations. |
| Earth permittivity and displacement current | Retained through complex material constants. | Definitions. |
| Range of validity | Infinite thin wire parallel to a planar interface. | Geometry. |
| Earth permeability ``μ_earth`` | Medium ratio retained. | Definitions. |
| Arrangement | Scalar self/mutual wire coefficients. | Formulation. |
| Earth structure | Two homogeneous half-spaces. | Figure 1. |
| Conductor and insulation geometry | Thin bare wire; no coating. | Model. |
| Constitutive and field assumptions | Linear isotropic harmonic media. | Derivation. |
| Conventions | Matrix admittance requires inversion after coefficient assembly. | Source TL assembly. |

**Expression.** The source's potential kernel uses the same analytical helper

```math
G(X,Y)\simeq\frac{n^2}{2(n^4-1)}[Q(bz)+Q(b\bar z)]
-\frac{P(b,z)+P(b,\bar z)-P(-b,z)-P(-b,\bar z)-n^2b[Q(-bz)+Q(-b\bar z)]}{2b(n^4-1)},
\qquad\text{(26)}
```

with ``z=k_0(Y+iX)``, ``b=i/\sqrt{1+n^2}``, ``Q(z)=e^{-z}E_1(-z)``, and ``P`` from (25). The assembled potential-coefficient matrix is inverted as a matrix to obtain the shunt parameters.

**Approximation.** qTEM/asymptotic analytical approximation; not a full-wave solved ``Γ``.

**Limitations.** Two homogeneous half-spaces and thin-wire geometry.

**Reference.** [Maaouni2001](@cite), equations (2)–(7), (22)–(26).

**Transcription source.** Original page-image verification of the full ``G`` expression.

## Source transcription

The impedance and potential outputs share helper functions but remain separate taxonomic records.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``G`` | unchanged | electric/potential correction kernel | source normalized |
| ``z,\bar z`` | unchanged | complex separation arguments | dimensionless |
| ``P,Q`` | unchanged | exponential-integral helpers | source convention |

## Evidence and approximation sources

The formulation was found during the broad lossy-dielectric/interface sweep.

## Limitations and discrepancies

The source's scalar exposition is not converted into entrywise inversion of matrix coefficients.
