# Hafner–Ferreira da Luz three-core cable FEM impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Three trefoil cores, sheaths, and common armor. Homogenized stranded cores/sheaths and finite nonconcentric metallic regions. |
| Calculated quantities | Self/mutual and sequence impedances of nonconcentric three-core cable metallic elements |
| Earth structure | None in the extracted formula. |
| Model and approximation | 2-D finite-element mesh and material homogenization; source explicitly contrasts the full coupling with concentric analytical simplifications. |
| Main source | Angelo A. Hafner, Mauricio V. Ferreira da Luz, and Walter P. Carpes Jr. (2015) |
| Citation key(s) | `:Hafner2015` |
| Evidence status | Original-paper text/page images checked |

**Description.** Magnetodynamic FEM extraction of all metallic-element couplings in a trefoil three-core submarine cable, including proximity and sheath/armor currents.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally invariant 2-D cross-section. | Abstract/model. |
| Air propagation constant ``γ_air`` | Not applicable. | Cable cross-section. |
| Earth propagation constant ``γ_earth`` | External earth is not the focus of this solve. | Model. |
| Earth permittivity and displacement current | Not included in series magnetodynamic solve. | (17). |
| Range of validity | Parallel three-core cable geometry resolved by mesh. | Abstract. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal matrix. |
| Arrangement | Three trefoil cores, sheaths, and common armor. | Figure 2/(8). |
| Earth structure | None in the extracted formula. | Scope. |
| Conductor and insulation geometry | Homogenized stranded cores/sheaths and finite nonconcentric metallic regions. | §§IV–VI. |
| Constitutive and field assumptions | Linear harmonic FEM; homogenized resistivity. | (1)–(2), (17). |
| Conventions | Unit-current excitation; induced voltage/current ratios define matrix entries. | §VI-B. |

**Expression.** Equation (8) is the full ``7×7`` metallic-element impedance matrix. For bonded sheaths,

```math
z_+=z_-=z_{aa}-z_{ab}-\frac{(z_{a1}-z_{a2})^2}{z_{11}-z_{12}},
\qquad
z_0=z_{aa}+2z_{ab}-\frac{z_{a1}+2z_{a2}}{z_{11}+2z_{a2}},
\qquad\text{(9)}
```

and for noncirculating sheath current, ``z_+=z_-=z_{aa}-z_{ab}``, ``z_0=z_{aa}+2z_{ab}`` (10). The FEM weak form is printed in (17).

**Approximation.** 2-D finite-element mesh and material homogenization; source explicitly contrasts the full coupling with concentric analytical simplifications.

**Limitations.** The formula covers a specific three-core topology and excludes the external earth-return component.

**Reference.** [Hafner2015](@cite), equations (8)–(10), (17).

**Transcription source.** Original IPST PDF text and equation layout checked.

## Source transcription

The denominator in the printed zero-sequence relation appears as shown and is not silently repaired.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``z_{aa},z_{ab}`` | unchanged | core self/mutual terms | ``Ω/m`` |
| ``z_{a1},z_{a2}`` | unchanged | core–sheath couplings | ``Ω/m`` |
| ``z_+,z_-,z_0`` | unchanged | sequence impedances | ``Ω/m`` |

## Evidence and approximation sources

The source provides the numerical formulation directly.

## Limitations and discrepancies

For implementation, derive the sequence quantities by Kron reduction of the published ``7\times7`` impedance matrix, as specified in `gaps.md`.
