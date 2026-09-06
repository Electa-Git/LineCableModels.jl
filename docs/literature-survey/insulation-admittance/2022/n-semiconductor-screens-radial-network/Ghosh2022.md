# Ghosh–Das N-screen radial admittance network

## Identity and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Loop coordinates across consecutive radial layers. Concentric semiconductor and insulation annuli. |
| Calculated quantities | Diagonal loop-admittance entries for arbitrary numbers of semiconducting screens and insulation annuli |
| Earth structure | None. |
| Model and approximation | Exact series-network algebra within the concentric homogeneous-annulus model. |
| Main source | Swarnankur Ghosh and Supriyo Das (2022) |
| Citation key(s) | `:Ghosh2022` |
| Evidence status | Original publication page images checked |

**Description.** General radial shunt network for ``N`` semiconductor screens, extending the earlier fixed two-screen series combination.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Parameters feed a TL propagation solve; not in radial formula. | Context. |
| Air propagation constant ``γ_air`` | Not applicable. | Insulation formula. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Insulation formula. |
| Earth permittivity and displacement current | Earth absent; semiconductor conduction and dielectric displacement retained. | (13)–(16). |
| Range of validity | Concentric cable with ``N`` screens. | Title/model. |
| Earth permeability ``μ_earth`` | Not applicable. | Shunt model. |
| Arrangement | Loop coordinates across consecutive radial layers. | (17). |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Concentric semiconductor and insulation annuli. | Figures. |
| Constitutive and field assumptions | Homogeneous annuli, quasistatic radial field. | (13)–(16). |
| Conventions | ``p=2,...,N-1`` indexes intermediate loops. | (17). |

**Expression.**

```math
Y_{l,11}=\left(\frac1{y_{sem,1}}+\frac1{y_{ins,1}}+\frac1{y_{sem,2}}\right)^{-1},
```

```math
Y_{l,pp}=\left(\frac1{y_{ins,p}}+\frac1{y_{sem,p+1}}\right)^{-1},
\quad p=2,\ldots,N-1,
\qquad Y_{l,NN}=y_{ins,N}.
\qquad\text{(17)}
```

```math
y_{sem,i}=\frac{j\omega\epsilon_0\epsilon_{sem,i}}{\ln(r_{out}/r_{in})},
\quad \epsilon_{sem,i}=\epsilon_{r,sem,i}+\frac1{j\omega\rho_{sem,i}},
\quad y_{ins,i}=\frac{j\omega\epsilon_0\epsilon_{r,ins,i}}{\ln(r_{out}/r_{in})}.
```

**Approximation.** Exact series-network algebra within the concentric homogeneous-annulus model.

**Limitations.** Diagonal radial-loop network; eccentric fields and insulation anisotropy are excluded.

**Reference.** [Ghosh2022](@cite), equations (13)–(17).

**Transcription source.** Original equation images; range of ``p`` and terminal ``N`` entry verified.

## Source transcription

The screen conductance is represented through complex relative permittivity as printed.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``y_sem,y_ins`` | unchanged | annular screen/insulation admittance | ``S/m`` |
| ``ρ_sem`` | unchanged | screen resistivity | ``Ω m`` |
| ``Y_l`` | unchanged | loop-coordinate shunt matrix | ``S/m`` |

## Evidence and approximation sources

This is distinct from the 2019 fixed-screen record because it supplies the explicit arbitrary-``N`` recursion.

## Limitations and discrepancies

Radius subscripts are abbreviated here; each occurrence must use the corresponding annulus radii from the source geometry.
