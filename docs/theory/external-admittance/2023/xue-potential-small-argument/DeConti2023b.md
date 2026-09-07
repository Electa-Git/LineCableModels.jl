# De Conti et al. small-argument underground potential coefficients

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Circular insulated cables; self uses total cable radius and equal depths. |
| Calculated quantities | Self and mutual Bessel-free small-argument entries of ``P_g``; assembled conversion ``Y_g=j\omega P_g^{-1}`` |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | Both Bessel functions in parent ``P_{g(m,n)}=j\omega[K_0(\gamma_1d)+\alpha K_0(\gamma_1D)]/[2\pi(\sigma_1+j\omega\varepsilon_1)]`` receive the leading small-argument expansion (7); logarithms are then collected to give (9). |
| Main source | A. De Conti, N. Duarte, R. Alipio, and O. E. Leal (2023) |
| Citation key(s) | `:DeConti2023b` |
| Evidence status | Original publication page images checked |

**Description.** Bessel-free small-argument approximation of each self or mutual ground-return potential coefficient for underground cables, assembled into ``P_g`` before matrix inversion.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed zero in the inherited compact Xue parent. | Equation-implied — (4)–(9), p. 2. |
| Air propagation constant ``γ_air`` | ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``. | Stated — (2), p. 2. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)}``. | Stated — (3), p. 2. |
| Earth permittivity and displacement current | Retained in both ``\gamma_1`` and the potential prefactor. | Stated — (3), (9). |
| Range of validity | Requires the relevant ``\gamma_1d`` and ``\gamma_1D`` arguments to satisfy ``0<z\ll1``. The paper reports accuracy comparable to the parent only up to about 1 MHz for the tested range. | Stated — (7) and §§3–5. |
| Earth permeability ``μ_earth`` | ``\mu_1=\mu_0``. | Stated — below (3). |
| Arrangement | Underground; self and mutual matrix entries. | Stated — Fig. 1 and §2.1. |
| Earth structure | Homogeneous earth below air. | Stated — §2. |
| Conductor and insulation geometry | Circular insulated cables; self uses total cable radius and equal depths. | Stated — below (3). |
| Constitutive and field assumptions | Linear homogeneous nonmagnetic earth; quasi-TEM parent. | Stated — §2.1. |
| Conventions | ``P_g`` is assembled entrywise and only then inverted through ``Y_g=j\omega P_g^{-1}``; ``j\omega`` convention; per-unit-length relation. | Stated — (6), p. 2. |

**Expression.** Proposed small-argument potential coefficient and assembly, equations (9) and (6), article p. 2.

```math
P_{g(m,n)}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left\{
\ln\left(\frac Dd\right)-(\alpha+1)
\left[\gamma_E+\ln\left(\frac{\gamma_1D}{2}\right)\right]
\right\},
\qquad\text{(9)}
```

```math
\begin{aligned}
\alpha&=\frac{\gamma_1^2-\gamma_0^2}{\gamma_1^2+\gamma_0^2} \\

D&=\sqrt{(h_m+h_n)^2+r^2} \\

\mathbf Y_g&=j\omega\mathbf P_g^{-1}.
\end{aligned}\qquad\text{(5--6)}
```

The entire matrix ``P_g`` must be assembled before inversion; no isolated mutual admittance equals ``1/P_{g(m,n)}``.

**Approximation.** Both Bessel functions in parent ``P_{g(m,n)}=j\omega[K_0(\gamma_1d)+\alpha K_0(\gamma_1D)]/[2\pi(\sigma_1+j\omega\varepsilon_1)]`` receive the leading small-argument expansion (7); logarithms are then collected to give (9).

**Limitations.** Mathematical small-argument restrictions apply to both direct and image distances. The tested admittance approximation loses accuracy above roughly 1 MHz. The formula covers only a homogeneous nonmagnetic earth and parallel infinite cables. Matrix inversion is required, and complex logarithm/root branches are not stated.

**Reference.** [DeConti2023b](@cite). A. De Conti, N. Duarte, R. Alipio, and O. E. Leal, “Small-argument analytical expressions for the calculation of the ground-return impedance and admittance of underground cables,” *Electric Power Systems Research* 220, 109299 (2023), equations (4)–(7), (9), article p. 2.

**Transcription source.** Original publication page image. The ``\ln(D/d)`` sign/order, ``-(\alpha+1)`` factor, bracket content, frequency normalization and matrix conversion were visually checked.

## Source transcription

The parent is

```math
P_{g(m,n)}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
[K_0(\gamma_1d)+\alpha K_0(\gamma_1D)],
\qquad\text{(4)}
```

and the only new operation is ``K_0(z)\approx-\ln(z/2)-\gamma_E``. This record does not independently invert any scalar entry.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\gamma`` in (7) | ``\gamma_E`` | Euler–Mascheroni constant | exact renaming |
| ``P_{g(m,n)}`` | unchanged | potential-matrix entry | source normalized |
| ``\mathbf P_g,\mathbf Y_g`` | bold in prose | assembled matrices | matrix inversion |
| ``\alpha`` | unchanged | air/earth bulk-constant contrast | dimensionless |
| ``d,D`` | unchanged | direct and image distances | m |

## Evidence and approximation sources

Equations (4)–(7) provide parent, material contrast and matrix conversion; (9) is the direct small-argument result. Frequency- and time-domain sections test rather than redefine it.

## Limitations and discrepancies

-
- The source calls these potential coefficients required for admittance, not individual admittance entries; this distinction is retained.
