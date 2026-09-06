# Pawlik–Woodhouse thin-wire lossy coating admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Scalar single conductor. Core radius ``a``, insulation outer radius ``b``. |
| Calculated quantities | Coaxial insulation admittance and its scalar series combination with external admittance |
| Earth structure | External admittance handled separately. |
| Model and approximation | Coaxial qTEM coating relation; source's one-wire scalar assembly is retained. |
| Main source | Brent Pawlik, Darren J. Woodhouse, and Terrence J. Summers (2020) |
| Citation key(s) | `:Pawlik2020` |
| Evidence status | Original publication page images checked |

**Description.** Complex-permittivity radial insulation admittance for a thin covered conductor, explicitly kept in series with the external potential path.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | qTEM extraction from full modal parent. | §III. |
| Air propagation constant ``γ_air`` | Not part of radial coating term. | (13). |
| Earth propagation constant ``γ_earth`` | Not part of radial coating term. | (13). |
| Earth permittivity and displacement current | Earth absent here; coating conductivity and displacement retained. | (13). |
| Range of validity | Thin infinite covered circular conductor. | Figure 1. |
| Earth permeability ``μ_earth`` | Not applicable. | Shunt dielectric term. |
| Arrangement | Scalar single conductor. | Derivation. |
| Earth structure | External admittance handled separately. | (11). |
| Conductor and insulation geometry | Core radius ``a``, insulation outer radius ``b``. | Figure 1. |
| Constitutive and field assumptions | Linear homogeneous isotropic coating, coaxial radial field. | §II. |
| Conventions | Scalar series combination ``1/Y_s=1/Y_i+1/Y_e``. | (11). |

**Expression.**

```math
Y_i=2\pi(\sigma_i+j\omega\epsilon_i)\left[\ln\frac ba\right]^{-1},
\tag{13}
```

```math
\frac1{Y_s}=\frac1{Y_i}+\frac1{Y_e}.
\tag{11}
```

**Approximation.** Coaxial qTEM coating relation; source's one-wire scalar assembly is retained.

**Limitations.** One homogeneous coating; no multilayer insulation recursion or multiconductor potential matrix.

**Reference.** [Pawlik2020](@cite), equations (11), (13).

**Transcription source.** Original page images; complex material factor and logarithmic inverse checked.

## Source transcription

The coating output is separated from the paper's external admittance record by taxonomy.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``a,b`` | unchanged | core and coating radii | m |
| ``σ_i,ε_i`` | unchanged | coating conductivity/permittivity | ``S/m``, ``F/m`` |
| ``Y_i`` | unchanged | radial insulation admittance per length | ``S/m`` |

## Evidence and approximation sources

The lossy dielectric term is stated separately from the cable application.

## Limitations and discrepancies

The source's scalar series relation must not be applied entrywise to matrices.
