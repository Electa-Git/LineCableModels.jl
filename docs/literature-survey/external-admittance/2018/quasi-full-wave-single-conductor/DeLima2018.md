# De Lima et al. quasi-full-wave single-conductor admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite circular conductor of radius ``r`` at height/depth ``h``; no finite insulation layer in the displayed formula. |
| Calculated quantities | Full-wave parent and quasi-full-wave per-unit-length admittance of one overhead or bare buried conductor, normalized as the scalar line admittance ``Y=\gamma/Z_c`` |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | qFW replaces the unknown longitudinal root by a predefined image-derived value inside the full-wave spectral quantities. It does not apply scalar reciprocals to potential coefficients, and it does not replace the remaining integral by a closed form. |
| Main source | A. C. S. de Lima, A. P. C. Magalhães, P. E. D. Rocha, R. A. Meyberg, and M. T. C. de Barros (2018) |
| Citation key(s) | `:DeLima2018` |
| Evidence status | Original publication page images checked; qFW image dependency unresolved |

**Description.** Scalar per-unit-length admittance of a single thin conductor parallel to a planar interface between two lossy media, defined by the source through its characteristic impedance and longitudinal modal constant and evaluated with either the full-wave solution or the prescribed quasi-full-wave estimate.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source ``\gamma`` is solved by the full-wave modal equation or replaced in qFW by an image-derived ``\bar\gamma``; longitudinal dependence is ``e^{-\gamma z}``. | Stated — (1), (7)–(9), (15)–(16), pp. 1874–1875. |
| Air propagation constant ``γ_air`` | Indexed ``\gamma_i=\sqrt{j\omega\mu_i(\sigma_i+j\omega\epsilon_i)}``; both media may be lossy. | Stated — p. 1874. |
| Earth propagation constant ``γ_earth`` | Same indexed definition; earth is medium 2 overhead and medium 1 in the buried application. | Stated — p. 1874 and §III-B. |
| Earth permittivity and displacement current | Retained in every bulk constant and in the printed admittance prefactor. | Stated — (13) and definitions, pp. 1874–1875. |
| Range of validity | Thin infinite single conductor, one interface, ``|\gamma_c|\gg|\gamma|`` when including conductor loss. Numerical tests are configuration-specific and do not state a universal frequency limit. | Stated — pp. 1874–1878. |
| Earth permeability ``μ_earth`` | ``\mu_1=\mu_2=\mu_0``. | Stated — p. 1874. |
| Arrangement | Single overhead or bare buried conductor; self only. | Stated — abstract, Fig. 1, and §III. |
| Earth structure | Two homogeneous half-spaces. | Stated — Fig. 1. |
| Conductor and insulation geometry | Infinite circular conductor of radius ``r`` at height/depth ``h``; no finite insulation layer in the displayed formula. | Stated — Fig. 1 and §II. |
| Constitutive and field assumptions | Linear isotropic lossy media; full-wave scalar/vector potentials. qFW approximates only the modal constant supplied to spectral functions. | Stated — §§II–III. |
| Conventions | ``e^{j\omega t}``, ``e^{-\gamma z}``; source voltage integrates transverse field from interface to conductor and keeps the interface potential reference; scalar per-unit-length output. | Stated — (1), (10)–(13), pp. 1874–1875. |

**Expression.** Source equation (13) with (14), evaluated for qFW with (15)–(16), printed p. 1875.

```math
Y=\frac{\gamma}{Z_c}=2\pi(\sigma_1+j\omega\epsilon_r\epsilon_1)[\Lambda_1-S_4],
\tag{13}
```

```math
S_4=2\int_0^\infty\frac{u_2}{u_1}
\frac{e^{-hu_1}-e^{-2hu_1}}{n^2u_1+u_2}\cos(r\lambda)\,d\lambda,
\qquad n=\frac{\gamma_2}{\gamma_1},
\tag{14}
```

```math
\Lambda_1=\Lambda(r,h),\qquad
\Lambda=K_0(\eta_1d)-K_0(\eta_1D),qquad
u_i=\sqrt{\lambda^2+\gamma_i^2-\gamma^2},
\quad \eta_1^2=\gamma_1^2-\gamma^2.
```

For qFW,

```math
u_i\approx\bar u_i=\sqrt{\lambda^2+\gamma_i^2-\bar\gamma^2},
\qquad
\eta_1\approx\bar\eta=\sqrt{\lambda^2-\bar\gamma^2},
\tag{15}
```

with ``\bar\gamma`` supplied by the image approximation and constrained by source equation (16), transcribed in the companion impedance record.

**Approximation.** qFW replaces the unknown longitudinal root by a predefined image-derived value inside the full-wave spectral quantities. It does not apply scalar reciprocals to potential coefficients, and it does not replace the remaining integral by a closed form.

**Limitations.** This is a one-conductor scalar admittance, not a multiconductor Maxwell-potential matrix. The source's image expression for ``\bar\gamma`` is not reprinted, so the qFW dependency is incomplete. The factor ``\epsilon_r\epsilon_1`` in (13) is reproduced exactly as printed even though the preceding general material definition is ``\epsilon_i=\epsilon_{ri}\epsilon_0``; the apparent notation inconsistency is not repaired. Root branches are not stated beside the formula.

**Reference.** [DeLima2018](@cite).  A. C. S. de Lima et al., “A Noniterative Approximation of a Full-Wave Model of Thin Wire Above and Buried in a Lossy Ground,” *IEEE Transactions on Electromagnetic Compatibility*, 60(6), 1873–1881 (2018), DOI `10.1109/TEMC.2017.2762241`, equations (10)–(16), p. 1875.

**Transcription source.** Original IEEE page image, printed p. 1875. Equation (13)'s literal material factor, the two separate exponentials in (14), ratios, roots, and qFW bars were visually checked. No inferred normalization or typographic correction is substituted.

## Source transcription

The source defines

```math
U=\int_0^h E_{y1}(r,\xi)\,d\xi,
\qquad
Z_c=\frac UI=\frac1I\left(\varphi_{1h}-\varphi_{10}+j\omega\int_0^hA_{1y}(r_j,\xi)\,d\xi\right),
\tag{10--11}
```

before stating ``Y=\gamma/Z_c``. Thus (13) is a source-defined admittance, not an inferred inversion of individual kernel entries.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Y`` | unchanged | scalar per-unit-length admittance | ``\mathrm S/\mathrm m`` |
| ``Z_c`` | unchanged | source characteristic impedance | ``\Omega`` |
| ``\gamma,\bar\gamma`` | unchanged | solved full-wave and prescribed qFW longitudinal constants | ``\mathrm m^{-1}`` |
| ``\gamma_i`` | unchanged | bulk medium constants | ``\mathrm m^{-1}`` |
| ``\Lambda_1,S_4`` | unchanged | source Bessel and Sommerfeld terms | source normalized |
| ``r,h`` | unchanged | conductor radius and interface distance | m |
| ``\lambda`` | unchanged | spectral variable | ``\mathrm m^{-1}`` |

No notation was renamed.

## Evidence and approximation sources

The voltage and normalization are explicit in (10)–(13). The complete full-wave parent and qFW prescription are equations (1)–(16); numerical sections compare full-wave, qFW, qTEM, and image variants without redefining the equation.

## Limitations and discrepancies

- The literal ``\epsilon_r\epsilon_1`` token in (13) conflicts with the surrounding indexed-permittivity notation; source clarification or erratum is needed.
- The equation is scalar. It supplies no mutual coefficient or matrix-assembly rule.
