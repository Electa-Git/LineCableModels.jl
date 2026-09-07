# Brito–Machado multipole sheath skin/proximity impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Arbitrary phase positions inside sheath. Solid circular cores and circular conducting sheath. |
| Calculated quantities | Three-phase cable series-impedance matrix including core/sheath skin and proximity effects |
| Earth structure | None. |
| Model and approximation | Multipole series is truncated numerically; perfect-sheath and finite-sheath assemblies are distinct source cases. |
| Main source | Ana Isabel Brito, V. Maló Machado, M. E. Almeida, and M. Guerreiro das Neves (2016) |
| Citation key(s) | `:Brito2016` |
| Evidence status | Original publication page images checked |

**Description.** Magnetic multipole expansion for arbitrary phase positions inside a cable sheath, with perfect- and nonperfect-sheath assemblies.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in cross-section diffusion. | Formulation. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal cable model. |
| Earth propagation constant ``γ_earth`` | External earth term is separate. | Assembly. |
| Earth permittivity and displacement current | Not part of conductor/sheath diffusion. | Scope. |
| Range of validity | Parallel round phase conductors within circular sheath. | Geometry. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal formula. |
| Arrangement | Arbitrary phase positions inside sheath. | Abstract/Fig. 1. |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Solid circular cores and circular conducting sheath. | Geometry. |
| Constitutive and field assumptions | Linear harmonic conductors; multipole Bessel solutions. | (4)–(11). |
| Conventions | Perfect-sheath formula precedes nonperfect-sheath corrections. | §§2–3. |

**Expression.** For a perfect sheath,

```math
Z_{kk}^{ps}=-\frac{j\omega\mu_k}{2\pi}
\frac{J_0(\bar\chi_k)}{\bar\chi_kJ_1(\bar\chi_k)}
+\frac{j\omega\mu_0}{2\pi}\ln\frac{r_0}{r_k}
+j\omega\frac{P_{kk}}{I_k},
```

```math
Z_{ki}^{ps}=\frac{j\omega\mu_0}{2\pi}\ln\frac{r_0}{|\bar w_{ki}|}
+j\omega\frac{P_{ki}}{I_i}.
```

The multipole potential ``P_k`` is the sum printed in (6), and (30)–(37) assemble the nonperfect-sheath terms.

**Approximation.** Multipole series is truncated numerically; perfect-sheath and finite-sheath assemblies are distinct source cases.

**Limitations.** Circular enclosure/cores; external earth and dielectric shunt parameters are outside this record.

**Reference.** [Brito2016](@cite), equations (4)–(37).

**Transcription source.** Original images; the mutual logarithmic term omitted in an earlier working summary is included here after direct verification.

## Source transcription

The harmonic coefficients ``C_p,D_p,U_{ki},V_k`` are defined by (4)–(11) and remain required dependencies.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``r_0,r_k`` | unchanged | sheath/core radii | m |
| ``\bar w_{ki}`` | unchanged | complex core separation | m |
| ``P_{ki}`` | unchanged | proximity multipole potential | source normalized |

## Limitations and discrepancies

No closed universal multipole truncation order is supplied.
