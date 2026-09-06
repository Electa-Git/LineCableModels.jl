# Pawlik–Woodhouse–Summers full-spectrum mixed mutual impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel thin conductors placed in opposite homogeneous half-spaces. |
| Calculated quantities | Full-spectrum cross-boundary mutual per-unit-length impedance |
| Earth structure | Two homogeneous half-spaces with independent conductivity, permittivity, and permeability. |
| Model and approximation | Integral representation within the full TM/TE model retaining ``Γ``; thin-wire and infinite-line assumptions remain. |
| Main source | B. Pawlik, D. Woodhouse, and T. J. Summers (2018) |
| Citation key(s) | `:Pawlik2018` |
| Evidence status | Original publication page images checked; apparent printed LHS index defect in (88) retained |

**Description.** Full-spectrum per-unit-length mutual impedance for two infinitely long, parallel thin conductors placed on opposite sides of an interface between homogeneous media with independently specified conductivity, permittivity, and permeability.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained in the source current ``I=I_Ae^{j\omega t-\Gamma z}``; free modes follow by solving the assembled system for ``\Gamma``. | Stated — p. 267 and discussion after (83), p. 272. |
| Air propagation constant ``γ_air`` | For medium 1, ``\gamma_1^2=\Gamma^2+k_1^2`` and ``u_1=\sqrt{\lambda^2-\gamma_1^2}``; medium 1 becomes air only under the later specialization. | Stated — (9), (39), and §V.A. |
| Earth propagation constant ``γ_earth`` | For medium 2, ``\gamma_2^2=\Gamma^2+k_2^2`` and ``u_2=\sqrt{\lambda^2-\gamma_2^2}``. | Stated — (21), (40). |
| Earth permittivity and displacement current | Retained through ``k_i=\omega\sqrt{\mu_i(\varepsilon_i-j\sigma_i/\omega)}``. | Stated — nomenclature and (5), (17). |
| Range of validity | Infinite-wire and thin-wire assumptions; the source requires decay before physical ends for an infinite-wire approximation and warns that homogeneous earth is limiting at power-system frequencies. No universal numerical frequency bound is supplied. | Stated — p. 273 after (91). |
| Earth permeability ``μ_earth`` | Independent ``\mu_2`` retained. | Stated — abstract, definitions, and (88)–(89). |
| Arrangement | Cross-boundary mutual term; ``j`` lies in medium 1 and ``p`` in medium 2. Reciprocal placements are stated equivalent after swapping medium and conductor subscripts. | Stated — paragraph before (88) and after (91), p. 273. |
| Earth structure | Two homogeneous half-spaces separated by a plane interface. | Stated — Fig. 1 and §II.A. |
| Conductor and insulation geometry | Infinitely long circular thin wires reduced to filaments for the source field; finite conductor conductance is included elsewhere in the full system. The cross-boundary external kernel contains no insulation region. | Stated — abstract and §II.A. |
| Constitutive and field assumptions | Linear, isotropic, homogeneous media; full TM and TE components; interface boundary conditions; wire finite conductance allowed in total system. | Stated — abstract and §§II–III. |
| Conventions | ``e^{j\omega t-\Gamma z}``; ``h_j,h_p`` are positive perpendicular distances into their respective media; ``d_{jp}`` is horizontal separation; roots are chosen for decay. | Stated — p. 267, (39)–(40), and (89). |

**Expression.** The source prints

```math
Z_{11}^{jk}=\frac{j\omega\mu_1}{\pi}
\left[Q_{12}^{jp}-jP_{12}^{jp}\right],
\qquad\text{(88)}
```

```math
Q_{12}^{jp}-jP_{12}^{jp}
=\int_0^\infty
\frac{e^{-(u_1h_j+u_2h_p)}}
{u_1+(\mu_1/\mu_2)u_2}
\cos(\lambda d_{jp})\,d\lambda,
\qquad\text{(89)}
```

with

```math
k_i=\omega\sqrt{\mu_i\left(\varepsilon_i-j\frac{\sigma_i}{\omega}\right)},
\qquad
\gamma_i^2=\Gamma^2+k_i^2,
\qquad
u_i=\sqrt{\lambda^2-\gamma_i^2}.
```

Equation (88)'s left-hand side is transcribed exactly; the surrounding paragraph and right-hand side identify the output as the cross-boundary ``j,p`` mutual term.

**Approximation.** Not an analytical approximation within the source's infinite thin-wire, two-homogeneous-half-space full-spectrum model. Numerical evaluation and solution of the modal system are still required. Thin-wire and infinite-line assumptions are physical reductions.

**Limitations.** Homogeneous half-spaces, parallel infinite wires, no end effects, and no layered earth. Equation (88) appears to retain the same-medium ``Z_{11}^{jk}`` left-hand side from (84), despite the cross-boundary paragraph and ``jp/12`` right-hand side; it is preserved as a suspected published index defect, not repaired.

**Reference.** [Pawlik2018](@cite), equations (88)–(89), printed p. 273.

**Transcription source.** Original IEEE page image. Medium ratios, exponent signs, conductor indices on the kernel, cosine argument, and the apparently inconsistent left-hand side of (88) were visually checked.

## Source transcription

The source first assembles an ``n``-wire system in (74)–(83). For ``j`` in medium 1 and ``p`` in medium 2, it then prints (88)–(89). It states after (91) that the reversed mixed quantities are equivalent after interchanging both material and conductor subscripts.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{11}^{jk}`` in (88) | unchanged | printed LHS; surrounding context identifies a cross-boundary mutual impedance | ``\Omega/\mathrm m``; index discrepancy retained |
| ``Q_{12}^{jp}-jP_{12}^{jp}`` | unchanged | mixed magnetic spectral coefficient | source normalization |
| ``h_j,h_p`` | unchanged | positive distances into media 1 and 2 | m |
| ``d_{jp}`` | unchanged | horizontal separation | m |
| ``k_i,\gamma_i,u_i`` | unchanged | material, longitudinally modified, and transverse wave numbers | ``\mathrm m^{-1}`` |

No notation was renamed.

## Evidence and approximation sources

The retained ``\Gamma`` distinguishes this full-spectrum formulation from the low-frequency Dawalibi–Southey and Pollaczek reductions. Martins-Britto et al. 2024 explicitly identify their mixed impedance equation with Pawlik's under common restrictions; this is corroboration, not a replacement transcription.

## Limitations and discrepancies

- **Suspected published defect:** (88) prints ``Z_{11}^{jk}``, while its introducing sentence, bracketed coefficient, and (89) use the cross-boundary indices ``jp`` and ``12``. No correction is supplied here.
- The source specifies decay selection for the transverse roots in the derivation, but no compact branch rule is repeated beside (89).
- Steady-state end effects of a horizontal buried wire are explicitly left for further investigation.
