# Mouhaidali–Chadebec multilayer HVDC cable FEM impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | ``k+1`` conductors with earth/reference reduction. Finite multizone cable geometry; conductive layers meshed. |
| Calculated quantities | Total cable series matrix in multilayer earth from loop-based FEM, including skin/proximity and earth return |
| Earth structure | Arbitrary meshed underground/seawater/seabed layers. |
| Model and approximation | Finite-element mesh and finite outer domain; source requires refinement across skin depth and reports a three-to-four-element practical rule. |
| Main source | A. Mouhaidali, O. Chadebec, S. Silvant, D. Tromeur-Dervout, and J.-M. Guichon (2018) |
| Citation key(s) | `:Mouhaidali2018` |
| Evidence status | Original-publication page image/text checked |

**Description.** Loop-excitation FEM extraction for HVDC cables embedded in multilayer underground/submarine media.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Cross-section is longitudinally homogeneous. | §III-A. |
| Air propagation constant ``γ_air`` | Quasistatic finite domain. | Formulation. |
| Earth propagation constant ``γ_earth`` | Conductive diffusion solved spatially. | (9). |
| Earth permittivity and displacement current | Neglected for series extraction. | §II-A. |
| Range of validity | Parallel cable geometry resolved by mesh and outer-boundary convergence. | §III. |
| Earth permeability ``μ_earth`` | Material-region permeability retained. | (7)–(9). |
| Arrangement | ``k+1`` conductors with earth/reference reduction. | (12)–(13). |
| Earth structure | Arbitrary meshed underground/seawater/seabed layers. | Abstract/case study. |
| Conductor and insulation geometry | Finite multizone cable geometry; conductive layers meshed. | §III. |
| Constitutive and field assumptions | Linear materials, harmonic magnetodynamic FEM. | (7)–(9). |
| Conventions | Reference-conductor/NRC reduction. | (13). |

**Expression.** The series solve uses

```math
\nabla\!\left(\frac1\mu\nabla\times\mathbf A\right)
+\sigma(j\omega\mathbf A+\nabla V)=\mathbf J,
\tag{9}
```

with loop extraction ``Z_{loop,i,j}=V/I`` and reference reduction

```math
Z_{ij}^{reduced}=Z_{ij}+Z_{k+1,k+1}-Z_{i,k+1}-Z_{j,k+1}.
\tag{13}
```

**Approximation.** Finite-element mesh and finite outer domain; source requires refinement across skin depth and reports a three-to-four-element practical rule.

**Limitations.** The FEM terminal matrix combines conductor and external fields; this record is placed here because the novel installation feature is multilayer earth return, not because internal loss disappears.

**Reference.** [Mouhaidali2018](@cite), equations (7)–(14).

**Transcription source.** Original paper, formulation and extraction sections.

## Source transcription

No analytical earth kernel is inferred from the numerical solve.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``A,V`` | unchanged | magnetic vector/electric scalar potentials | SI |
| ``k+1`` | unchanged | reference conductor/ground domain | index |
| ``Z_reduced`` | unchanged | terminal cable matrix | ``Ω/m`` |

## Evidence and approximation sources

The paper gives explicit parameter-extraction equations.

## Limitations and discrepancies

Conductor and earth contributions cannot be separated from the printed total FEM result without additional solves.
