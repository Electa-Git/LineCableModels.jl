# Moghram three-layer earth DCFT impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Self and mutual overhead conductors. Filament source with radius used for self evaluation. |
| Calculated quantities | Overhead-line self and mutual impedances over three-layer earth, with two-layer and homogeneous limits |
| Earth structure | Three conducting layers below air; reductions to two/one. |
| Model and approximation | Displacement and longitudinal variation are neglected; within that model the DCFT layer solution is not a fitted approximation. |
| Main source | I. S. Moghram (1998), extending Wedepohl–Wasley |
| Citation key(s) | `:Moghram1998` |
| Evidence status | Original publication equations and text checked |

**Description.** Double-complex-Fourier-transform derivation of a three-layer-earth overhead impedance, explicitly reduced to Wedepohl's two-layer and Carson's homogeneous cases.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Earth currents/fields do not vary along propagation direction. | Stated assumptions. |
| Air propagation constant ``γ_air`` | Longitudinal variation neglected. | §2. |
| Earth propagation constant ``γ_earth`` | ``m_k^2=jωμ_kσ_k``. | (10). |
| Earth permittivity and displacement current | Neglected. | Stated assumptions. |
| Range of validity | Infinite overhead thin conductor, power-frequency conductive earth. | §2. |
| Earth permeability ``μ_earth`` | Independent layer permeabilities retained. | Abstract/(10). |
| Arrangement | Self and mutual overhead conductors. | (16)–(17). |
| Earth structure | Three conducting layers below air; reductions to two/one. | Figure 1, §§3–4. |
| Conductor and insulation geometry | Filament source with radius used for self evaluation. | (1), (16). |
| Constitutive and field assumptions | Linear layers; earth electric field parallel to conductor. | Stated assumptions. |
| Conventions | Sinusoidal source and source DCFT definitions (4)–(5). | §2. |

**Expression.** The air field is obtained in (13) from the DCFT solution with the three-layer interface constant ``A`` defined by (14)–(15). The terminal formulas are

```math
Z_s=-\frac{E_1(s-r,h)}{I},
\qquad
Z_m=-\frac{E_1(s_2,h_2)}{I},
\tag{16–17}
```

where ``E_1`` contains the printed direct logarithm plus the semi-infinite three-layer spectral integral. Setting layers 3 and 4 equal gives the source's two-layer coefficient (18); making all earth layers equal gives (19) and Carson's limit.

**Approximation.** Displacement and longitudinal variation are neglected; within that model the DCFT layer solution is not a fitted approximation.

**Limitations.** Three earth layers are printed; arbitrary-layer extension is described but not given as a complete recursion.

**Reference.** [Moghram1998](@cite), equations (1)–(19).

**Transcription source.** Original pages; self/mutual definitions and layer reductions checked. The long coefficient is retained by locator rather than risk OCR repair.

## Source transcription

Moghram's ``A`` is an independent three-layer extension, not merely the two-layer Wedepohl equation.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``m_k`` | unchanged | layer diffusion constant | ``m^{-1}`` |
| ``E_1`` | unchanged | longitudinal air electric field | ``V/m`` |
| ``Z_s,Z_m`` | unchanged | self/mutual p.u.l. impedance | ``Ω/m`` |

## Evidence and approximation sources

The source gives the three-layer DCFT construction explicitly.

## Limitations and discrepancies

OCR of equations (13)–(15) is poor; page images remain the authority for implementation.
