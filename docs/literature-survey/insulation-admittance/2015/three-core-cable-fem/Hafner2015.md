# Hafner–Ferreira da Luz three-core cable FEM admittance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Three cores/sheaths plus common armor reference. Nonconcentric trefoil with homogenized core insulation/screen. |
| Calculated quantities | Core–sheath, sheath–sheath, and sheath–armor admittance matrix of a nonconcentric three-core cable |
| Earth structure | None. |
| Model and approximation | 2-D nodal FEM plus homogenized screen/insulation permittivity; analytical annulus used only where concentric. |
| Main source | Angelo A. Hafner, Mauricio V. Ferreira da Luz, and Walter P. Carpes Jr. (2015) |
| Citation key(s) | `:Hafner2015` |
| Evidence status | Original-paper text/page images checked |

**Description.** Electrostatic FEM shunt extraction for nonconcentric cable dielectric paths, with an analytical core–sheath annulus and full matrix assembly.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not part of cross-sectional electrostatic extraction. | §VII–VIII. |
| Air propagation constant ``γ_air`` | Not applicable. | Cable insulation. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Extracted formula. |
| Earth permittivity and displacement current | Earth absent; insulation displacement retained. | (11), (16). |
| Range of validity | Three-core cable dielectric geometry representable by 2-D mesh. | Abstract. |
| Earth permeability ``μ_earth`` | Not applicable. | Electric solve. |
| Arrangement | Three cores/sheaths plus common armor reference. | (12). |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Nonconcentric trefoil with homogenized core insulation/screen. | §§V,VII. |
| Constitutive and field assumptions | Linear quasistatic scalar-potential FEM. | (16). |
| Conventions | Maxwell/nodal matrix signs use negative mutual admittances. | (12)–(15). |

**Expression.** The analytical core–sheath branch is

```math
y_{12}=\frac{2\pi(\sigma_{cs}+j\omega\varepsilon'_{cs})}{\ln(r_s/r_c)},
\qquad\text{(11)}
```

and the source assembles the six-node core/sheath matrix in (12), with

```math
y_{cc}=-y_{cs}=y_{a1},\qquad
y_{ss}=y_{11}=y_{a1}+y_{1g}+y_{12}+y_{13},\qquad
y'_{ss}=-y_{12}=-y_{13}.
\qquad\text{(13–15)}
```

**Approximation.** 2-D nodal FEM plus homogenized screen/insulation permittivity; analytical annulus used only where concentric.

**Limitations.** Specific three-core/armor topology; external earth potential is not derived.

**Reference.** [Hafner2015](@cite), equations (11)–(16).

**Transcription source.** Original IPST PDF equations checked.

## Source transcription

The full nonconcentric paths are numerical, not inferred coaxial substitutions.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``y_cc,y_cs`` | unchanged | core self/core–sheath mutual | ``S/m`` |
| ``y_ss,y'_ss`` | unchanged | sheath self/mutual | ``S/m`` |
| ``ε'_cs`` | unchanged | homogenized insulation permittivity | ``F/m`` |

## Evidence and approximation sources

The source derives both impedance and admittance.

## Limitations and discrepancies

The source's homogenization must be recomputed for other screen/insulation radii.
