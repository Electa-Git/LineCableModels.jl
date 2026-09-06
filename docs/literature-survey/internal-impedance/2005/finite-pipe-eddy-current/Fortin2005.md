# Fortin–Yang finite-pipe eddy-current impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Arbitrary core locations within pipe. Circular pipe, core axes; coatings do not carry pipe current. |
| Calculated quantities | Pipe-type cable self and mutual impedance from return and eddy currents in a finite wall |
| Earth structure | None in displayed contribution. |
| Model and approximation | Infinite harmonic sums require truncation; vector potential is two-dimensional and pipe material is unsaturated. |
| Main source | Simon Fortin, Y. Yang, J. Ma, and F. P. Dawalibi (2005) |
| Citation key(s) | `:Fortin2005` |
| Evidence status | Original-page image verified |

**Description.** Vector-potential integration of finite-pipe return and eddy currents, yielding direct terminal self/mutual voltage drops.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in two-dimensional pipe cross-section. | Formulation. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal pipe formula. |
| Earth propagation constant ``γ_earth`` | Earth term is separate. | Scope. |
| Earth permittivity and displacement current | Not part of pipe term. | Scope. |
| Range of validity | Multiple cables inside a circular finite-thickness pipe. | Title/model. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal formula. |
| Arrangement | Arbitrary core locations within pipe. | Figure 1. |
| Earth structure | None in displayed contribution. | Scope. |
| Conductor and insulation geometry | Circular pipe, core axes; coatings do not carry pipe current. | Geometry. |
| Constitutive and field assumptions | Linear nonsaturated magnetic or nonmagnetic pipe. | Text before §3. |
| Conventions | Unit source current may be used for impedance extraction. | §2. |

**Expression.** With the vector potential from (3)–(11),

```math
Z_{self}=j\omega[A(a,0)-A(b,0)]-E_z(a,0),
\tag{12}
```

```math
Z_{mutual}=j\omega[A(a,0)-A(r_j,\alpha)]-E_z(a,0).
\tag{13}
```

The finite-wall current density uses ``p=\sqrt{j\omega\mu\sigma}`` and the Bessel/harmonic expansion in (1)–(5).

**Approximation.** Infinite harmonic sums require truncation; vector potential is two-dimensional and pipe material is unsaturated.

**Limitations.** External earth return and cable-core internal impedance are appended separately.

**Reference.** [Fortin2005](@cite), equations (1)–(13).

**Transcription source.** Original-author conference PDF images; signs and evaluation points verified.

## Source transcription

The terminal formulas are retained with their dependent current/vector-potential equations.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``A`` | unchanged | longitudinal magnetic vector potential | ``Wb/m`` |
| ``E_z`` | unchanged | longitudinal pipe electric field | ``V/m`` |
| ``p`` | unchanged | pipe diffusion root | ``m^{-1}`` |

## Limitations and discrepancies

The inspected conference publication gives the year as 2005.
