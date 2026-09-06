# Dawalibi–Southey low-frequency mixed external impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Long parallel thin conductors on opposite sides of the air/soil interface. |
| Calculated quantities | Low-frequency mixed mutual per-unit-length external impedance |
| Earth structure | Homogeneous soil half-space below air. |
| Model and approximation | Integral representation after the source neglects longitudinal propagation in the observation zone. The technical-manual reproduction's dimensionally conflicting printed ``\gamma_i`` definition is retained. |
| Main source | F. P. Dawalibi and R. D. Southey (1989) |
| Citation key(s) | `:Dawalibi1989` |
| Evidence status | Original-author technical-manual reproduction checked against page images; original IEEE equation page not independently matched |

**Description.** Low-frequency per-unit-length mutual external impedance for two parallel thin conductors on opposite sides of a planar air/soil interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Propagation effects within the observation zone are neglected. | Stated — §4.3 immediately before (7)–(9), Technical Manual 3.2 p. 53. |
| Air propagation constant ``γ_air`` | The inspected reproduction literally prints ``\gamma_0=j\omega\mu_0\theta_0``, with ``\theta_0=\sigma_0+j\omega\varepsilon_0``; its dimensional conflict with (3) remains unresolved. | Stated — (4)–(5), p. 53. |
| Earth propagation constant ``γ_earth`` | The inspected reproduction literally prints ``\gamma_1=j\omega\mu_1\theta_1``, with ``\theta_1=\sigma_1+j\omega\varepsilon_1``; no square or square root is supplied on that page. | Stated — (4)–(5), p. 53. |
| Earth permittivity and displacement current | Retained through ``\theta_1=\sigma_1+j\omega\varepsilon_1``. | Stated — (4)–(5), p. 53. |
| Range of validity | The source calls neglect of propagation effects legitimate at low frequencies; no numerical frequency bound is supplied. | Stated — §4.3, p. 53. |
| Earth permeability ``μ_earth`` | Independent ``\mu_1`` is retained. | Equation-implied — prefactor, denominator, and ``\gamma_1`` in (8). |
| Arrangement | Mixed mutual interaction: one conductor in air and one in soil. | Stated — heading before (8), p. 53. |
| Earth structure | One homogeneous soil half-space below air. | Stated — Fig. 2 referenced by §4.3 and medium indices in (7)–(9). |
| Conductor and insulation geometry | Long cylindrical conductors reduced to axial filaments for the external interaction; no insulation region enters (8). | Stated — opening paragraphs of §4.3, p. 53. |
| Constitutive and field assumptions | Linear homogeneous media; thin-wire external field; azimuthally symmetric longitudinal currents; local soil changes and leakage losses are excluded from the kernel. | Stated — §4.3, p. 53. |
| Conventions | ``j=\sqrt{-1}``; ``h_i`` and ``h_k`` are positive distances into air and soil; ``x`` is lateral separation; result is per unit length. | Equation-implied — (8) and its mixed-arrangement heading. |

**Expression.** The source prints the mixed mutual external impedance as equation (8):

```math
Z_{ik}^{e}=\frac{j\omega\mu_0\mu_1}{\pi}
\int_0^\infty
\frac{e^{-\alpha_0h_i}e^{-\alpha_1h_k}}
{\alpha_0\mu_1+\alpha_1\mu_0}
\cos(x\lambda)\,d\lambda,
\qquad\text{(8)}
```

with

```math
\begin{aligned}
\alpha_i&=\left(\lambda^2+\gamma_i^2\right)^{1/2} \\
\theta_i&=\sigma_i+j\omega\varepsilon_i \\
\gamma_i&=j\omega\mu_i\theta_i.
\end{aligned}\qquad\text{(3--5)}
```

**Approximation.** Not an analytical approximation after the source imposes its thin-wire, homogeneous-medium, and low-frequency zero-longitudinal-propagation model. Those model reductions precede equation (8).

**Limitations.** The source supplies a mixed mutual term, not a finite-radius self prescription or a layered-earth result. Its low-frequency qualification has no printed quantitative bound. The original IEEE equation page was not available as the controlling image; fidelity is established only against the explicitly identified original-author reproduction.

**Reference.** [Dawalibi1989](@cite), equation (8), as reproduced in Technical Manual 3.2 p. 53; journal DOI `10.1109/61.32680`.

**Transcription source.** Original-author later reproduction. The prefactor, two exponential factors, permeability-weighted denominator, cosine argument, and definitions (3)–(5) were checked against the page image. Equation (5) visibly lacks a square on ``\gamma_i``; that token is retained rather than repaired. This verifies the reproduction, not token-for-token identity with an unseen original IEEE page.

## Source transcription

The published reproduction separates self impedance into internal and external parts before presenting three placement cases. The mixed formula is equation (8), reproduced literally above. Section 4.3 explicitly says that propagation effects are neglected in the observation zone for the low-frequency use of (7)–(9).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{ik}^{e}`` | unchanged | mixed mutual external impedance | ``\Omega/\mathrm m`` |
| ``h_i,h_k`` | unchanged | distances from the interface into air and soil | m, positive |
| ``x`` | unchanged | horizontal separation | m |
| ``\alpha_0,\alpha_1`` | unchanged | transverse spectral roots | ``\mathrm m^{-1}`` |
| ``\theta_i`` | unchanged | complex conductivity | ``\mathrm S/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

The inspected source is an equation-preserving technical-manual reproduction carrying the paper title/byline and the authors' formula. The DOI and journal metadata identify the 1989 publication. Martins-Britto et al. later state that their generalized equation (5) is identical to Dawalibi–Southey after the corresponding restrictions; that later statement is corroborative and is not used to rewrite equation (8).

## Limitations and discrepancies

- Original IEEE page-image equivalence remains unverified.
- **Suspected published/reproduction defect:** (3) uses ``\lambda^2+\gamma_i^2`` while (5) prints ``\gamma_i=j\omega\mu_i\theta_i``. This is dimensionally inconsistent with the usual propagation-constant definition, but no square or square root is inserted here.
- The source gives no square-root branch beside (3)–(5); decay is equation-implied by the negative exponentials.
- No mixed external-admittance expression is supplied in this source.
