# De Lima et al. quasi-full-wave single-conductor admittance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite circular conductor of radius ``r`` at height/depth ``h``; no finite insulation layer in the displayed formula. |
| Calculated quantities | Full-wave parent and quasi-full-wave per-unit-length admittance of one overhead or bare buried conductor, normalized as the scalar line admittance ``Y=\gamma/Z_c`` |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | qFW replaces the unknown longitudinal root by a predefined image-derived value inside the full-wave spectral quantities. It does not apply scalar reciprocals to potential coefficients, and it does not replace the remaining integral by a closed form. |
| Main source | A. C. S. de Lima, A. P. C. Magalhães, P. E. D. Rocha, R. A. Meyberg, and M. T. C. de Barros (2018) |
| Citation key(s) | `:DeLima2018` |
| Evidence status | Original publication page images and Appendix D image expressions checked |

**Description.** Scalar per-unit-length admittance of a single thin conductor parallel to a planar interface between two lossy media, defined by the source through its characteristic impedance and longitudinal modal constant and evaluated with either the full-wave solution or the prescribed quasi-full-wave estimate.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source ``\gamma`` is solved by the full-wave modal equation or replaced in qFW by an image-derived ``\bar\gamma``; longitudinal dependence is ``e^{-\gamma z}``. | Stated — (1), (7)–(9), (15)–(16), pp. 1874–1875. |
| Air propagation constant ``γ_air`` | Indexed ``\gamma_i=\sqrt{j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)}``; both media may be lossy. | Stated — p. 1874. |
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
Y=\frac{\gamma}{Z_c}=2\pi(\sigma_1+j\omega\varepsilon_r\varepsilon_1)[\Lambda_1-S_4],
\qquad\text{(13)}
```

```math
\begin{aligned}
S_4&=2\int_0^\infty\frac{u_2}{u_1}
\frac{e^{-hu_1}-e^{-2hu_1}}{n^2u_1+u_2}\cos(r\lambda)\,d\lambda \\
n&=\frac{\gamma_2}{\gamma_1},
\end{aligned}\qquad\text{(14)}
```

```math
\begin{aligned}
\Lambda_1&=\Lambda(r,h) \\
\Lambda&=K_0(\eta_1d)-K_0(\eta_1D) \\
u_i&=\sqrt{\lambda^2+\gamma_i^2-\gamma^2} \\
\eta_1^2&=\gamma_1^2-\gamma^2.
\end{aligned}
```

For qFW,

```math
\begin{aligned}
u_i&\approx\bar u_i=\sqrt{\lambda^2+\gamma_i^2-\bar\gamma^2} \\
\eta_1&\approx\bar\eta=\sqrt{\lambda^2-\bar\gamma^2},
\end{aligned}\qquad\text{(15)}
```

with ``\bar\gamma`` supplied by the image approximation and constrained by source equation (16), transcribed in the companion impedance record.

**Approximation.** qFW replaces the unknown longitudinal root by a predefined image-derived value inside the full-wave spectral quantities. It does not apply scalar reciprocals to potential coefficients, and it does not replace the remaining integral by a closed form.

**Limitations.** This is a one-conductor scalar admittance, not a multiconductor Maxwell-potential matrix. Appendix D supplies the image estimate; the numerical interpretation below fixes the inverse and material factor from the voltage definition without requiring an erratum.

**Reference.** [DeLima2018](@cite).  A. C. S. de Lima et al., “A Noniterative Approximation of a Full-Wave Model of Thin Wire Above and Buried in a Lossy Ground,” *IEEE Transactions on Electromagnetic Compatibility*, 60(6), 1873–1881 (2018), DOI `10.1109/TEMC.2017.2762241`, equations (10)–(16), p. 1875.

**Transcription source.** Original IEEE page image, printed p. 1875. Equation (13)'s literal material factor, the two separate exponentials in (14), ratios, roots, and qFW bars were visually checked. No inferred normalization or typographic correction is substituted.

## Source transcription

The source defines

```math
\begin{aligned}
U&=\int_0^h E_{y1}(r,\xi)\,d\xi \\
Z_c&=\frac UI=\frac1I\left(\varphi_{1h}-\varphi_{10}+j\omega\int_0^hA_{1y}(r_j,\xi)\,d\xi\right),
\end{aligned}\qquad\text{(10--11)}
```

before stating ``Y=\gamma/Z_c``. Thus (13) is a source-defined admittance, not an inferred inversion of individual kernel entries.

## Numerical evaluation

Appendix D supplies the image expressions needed by qFW. With
``N=n^2=\gamma_2^2/\gamma_1^2``,
``\beta=\sqrt{\gamma_2^2-\gamma_1^2}``,
``D=\sqrt{4h^2+r^2}``, and
``A=(N+1)/(\beta D)``, the author writes

```math
\begin{aligned}
Z_{\mathrm{image}}&=z_i+\frac{j\omega\mu_0}{2\pi}
\left[\ln\frac{2h}{r}+\bar S_1-\bar S_4-\bar S_2\right],\\
Y_{\mathrm{image}}&=2\pi(\sigma_1+j\omega\varepsilon_1)
\left[\ln\frac{2h}{r}-\bar S_4\right]^{-1}.
\end{aligned}\qquad\text{(23)}
```

```math
\begin{aligned}
\bar S_1&=\ln\left(1+\frac{2}{\beta D}\right),\\
\bar S_2&=\frac{2}{N+1}\ln(1+A),\\
\bar S_4&=2\ln2+\frac{2N}{N+1}\ln\frac{1+A}{1+2A}.
\end{aligned}\qquad\text{(24)}
```

The prescribed estimate is ``\bar\gamma=\sqrt{Z_{\mathrm{image}}Y_{\mathrm{image}}}``,
with nonnegative attenuation. The internal contribution ``z_i`` is
taken from the selected conductor formula. The external matrix term is
``Z_{\mathrm{image}}-z_i``; internal impedance is assembled separately.

For integral evaluation, definition (2) fixes
``\eta_1^2=\gamma_1^2-\bar\gamma^2``. The spectral variable in the
printed radial expression (15) is not used. Independently integrating
the potentials in (3) according to the voltage definition (10)–(11) gives

```math
\begin{aligned}
Z_c&=\frac{\gamma}{2\pi\kappa_1}(\Lambda_1-S_4),\\
Y&=\frac{2\pi\kappa_1}{\Lambda_1-S_4},\\
P&=\frac{j\omega}{2\pi\kappa_1}(\Lambda_1-S_4),\qquad
\kappa_1=\sigma_1+j\omega\varepsilon_1.
\end{aligned}
```

This inverse and material normalization also agree with (23).
The printed equations (13) and (15) remain separate from this
numerical interpretation. At ``\eta_1=0``, the Bessel difference
is evaluated as ``\ln(D/r)``.

The default implementation evaluates (23)–(24).
The `approximation=:quasi_full_wave` selection evaluates (12) and
the potential coefficient above at the explicitly supplied
``k_x=j\bar\gamma``. The image estimate need not satisfy the exact
modal equation (7); no automatic full-wave root selection is implied.
Each coefficient is scalar, and a second exterior conductor is unsupported.

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

The equation is scalar and supplies no mutual coefficient. The implemented potential follows the source voltage definition and is assembled as ``Y=j\omega P^{-1}``; it does not prescribe entrywise inversion of a multiconductor matrix.
