# De Silva–Shafieipour finite-pipe appendix formulation

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Self/mutual core pairs inside pipe. Circular pipe radii ``r_{p1},r_{p2}``; core axes. |
| Calculated quantities | Finite-pipe self/mutual component and cable-impedance block assembly |
| Earth structure | External component independent. |
| Model and approximation | Infinite harmonic sum and the paper's classical/numerical hybrid assembly; numerical truncation is required. |
| Main source | H. M. J. De Silva and M. Shafieipour (2019) |
| Citation key(s) | `:DeSilva2019` |
| Evidence status | Original equation-page verified |

**Description.** Classical/numerical hybrid paper's explicit correction and matrix assembly for cores inside a finite conducting pipe.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in pipe cross-section. | Appendix. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal term. |
| Earth propagation constant ``γ_earth`` | External earth matrix separate. | (A3)–(A5). |
| Earth permittivity and displacement current | Not part of pipe diffusion. | Scope. |
| Range of validity | Multiple circular cores in circular finite-thickness pipe. | Appendix geometry. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal term. |
| Arrangement | Self/mutual core pairs inside pipe. | (A1)–(A5). |
| Earth structure | External component independent. | Assembly. |
| Conductor and insulation geometry | Circular pipe radii ``r_{p1},r_{p2}``; core axes. | (A6). |
| Constitutive and field assumptions | Linear harmonic cylindrical diffusion. | Appendix. |
| Conventions | Matrix blocks assembled exactly as (A3)–(A5). | Source. |

**Expression.**

```math
Z_{pjk}=Q_{jk}+J_{jk},
\tag{A1}
```

```math
J_{pjk}=2\mu\sum_{n=1}^{\infty}
\frac{C_n}{n(1+\mu_p)+x_1I_{n-1}(x_1)/K_n(x_1)},
\tag{A2}
```

with pipe surface term

```math
Z_1=\frac{m_p\rho_p}{2\pi r_{p1}}
\frac{K_1(m_pr_{p2})I_0(m_pr_{p1})+K_0(m_pr_{p1})I_1(m_pr_{p2})}
{K_1(m_pr_{p1})I_1(m_pr_{p2})-K_1(m_pr_{p2})I_1(m_pr_{p1})}.
\tag{A6}
```

**Approximation.** Infinite harmonic sum and the paper's classical/numerical hybrid assembly; numerical truncation is required.

**Limitations.** Appendix notation depends on definitions in the cited pipe model; earth return is separate.

**Reference.** [DeSilva2019](@cite), appendix (A1)–(A6).

**Transcription source.** Original IPST appendix images.

## Source transcription

The matrix blocks (A3)–(A5) are dependencies but are not expanded here to avoid obscuring the new pipe correction.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``m_p`` | unchanged | pipe diffusion root | ``m^{-1}`` |
| ``r_{p1},r_{p2}`` | unchanged | pipe radii | m |
| ``J_{jk}`` | unchanged | harmonic pipe correction | ``Ω/m`` |

## Limitations and discrepancies

Attribution is to the inspected 2019 paper; parent classical terms retain their own earlier attributions.
