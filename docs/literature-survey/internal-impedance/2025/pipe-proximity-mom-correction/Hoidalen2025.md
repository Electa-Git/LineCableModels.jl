# Høidalen–Høyer-Hansen pipe proximity correction

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Multiple cores in pipe; self/mutual entries. Circular core(s), circular finite-wall pipe. |
| Calculated quantities | Incremental proximity correction for cores inside a finite-thickness pipe |
| Earth structure | None. |
| Model and approximation | Finite MoM harmonic count ``N_p``; the subtraction isolates the incremental proximity part. |
| Main source | Hans Kristian Høidalen, Martin Høyer-Hansen, Felipe Camara Neto, and Claus Leth Bak (2025) |
| Citation key(s) | `:Hoidalen2025` |
| Evidence status | Original page image verified |

**Description.** MoM-derived incremental proximity correction added to an existing finite-pipe impedance without double-counting the zero-harmonic base.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in pipe cross-section extraction. | Method. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal pipe term. |
| Earth propagation constant ``γ_earth`` | Earth return separate. | Assembly. |
| Earth permittivity and displacement current | Not part of pipe correction. | Scope. |
| Range of validity | Cores inside a circular conducting pipe represented by MoM harmonics. | Geometry. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal formula. |
| Arrangement | Multiple cores in pipe; self/mutual entries. | Equations. |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Circular core(s), circular finite-wall pipe. | Figure. |
| Constitutive and field assumptions | Linear harmonic conductor diffusion and surface harmonics. | Method. |
| Conventions | ``N_p=0`` is the zero-harmonic reference. | (19)–(21). |

**Expression.**

```math
Z(N_p)=Z_i+\Delta Z_{prox},
\qquad
Z(N_p=0)=Z_i+\frac{j\omega\mu_0}{2\pi}\ln\frac{r_{p1}}r,
\qquad\text{(19–20)}
```

```math
Z_{p,i}=Z_{pi}+Z(N_p)-Z(N_p=0)=Z_{pi}+\Delta Z_{prox},
\qquad\text{(21)}
```

```math
Z_{pi}=\frac{j\omega\mu_0}{2\pi}\frac{\mu_{rp}}{x_1}
\frac{I_0(x_1)K_1(x_2)+I_1(x_2)K_0(x_1)}
{I_1(x_2)K_1(x_1)-I_1(x_1)K_1(x_2)}.
\qquad\text{(22)}
```

**Approximation.** Finite MoM harmonic count ``N_p``; the subtraction isolates the incremental proximity part.

**Limitations.** Requires an existing compatible base pipe model and circular geometry.

**Reference.** [Hoidalen2025](@cite), equations (19)–(22).

**Transcription source.** Original equation images; signs, Bessel ordering, and subtraction verified.

## Source transcription

The record does not rename ``ΔZ_prox`` by cable application or earth model.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``N_p`` | unchanged | pipe MoM harmonic order | integer |
| ``x_1,x_2`` | unchanged | pipe diffusion arguments | dimensionless |
| ``Z_{pi}`` | unchanged | finite-pipe base input term | ``Ω/m`` |

## Evidence and approximation sources

Found during the newer-papers proximity/pipe sweep.

## Limitations and discrepancies

The 2025 source is represented by the accessible preprint/publication version; later pagination may differ.
