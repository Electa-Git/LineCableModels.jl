# De Conti et al. small-argument underground impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular insulated cable positions; self replaces ``r`` by total cable radius including jacket and sets equal depths. |
| Calculated quantities | Self and mutual Bessel-free small-argument approximation of the 2023 De Conti–Duarte–Alipio closed form |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | The only additional operation on parent equation (1) is ``K_0(z)\approx-\ln(z/2)-\gamma_E`` for ``0<z\ll1``. All other terms are retained unchanged. |
| Main source | A. De Conti, N. Duarte, R. Alipio, and O. E. Leal (2023) |
| Citation key(s) | `:DeConti2023b` |
| Evidence status | Original publication page images checked |

**Description.** Bessel-free small-argument approximation of the self or mutual per-unit-length ground-return impedance for insulated cables buried in homogeneous earth.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed zero, inherited from the compact Xue parent. | Stated by parent source and equation-implied — (1)–(8), p. 2. |
| Air propagation constant ``γ_air`` | ``\gamma_0=j\omega\sqrt{\mu_0\epsilon_0}``. | Stated — (2), p. 2. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_1(\sigma_1+j\omega\epsilon_1)}``. | Stated — (3), p. 2. |
| Earth permittivity and displacement current | Retained through ``\epsilon_1`` in ``\gamma_1``. | Stated — (3). |
| Range of validity | Mathematical condition ``0<\gamma_1d\ll1``. Frequency-domain comparisons find performance comparable to the parent mainly up to 1–2 MHz and state strict practical validity “up to 1 MHz or so”; higher-frequency transient agreement is configuration-specific. | Stated — (7), p. 2 and §§3–5. |
| Earth permeability ``μ_earth`` | ``\mu_1=\mu_0``. | Stated — text below (3), p. 2. |
| Arrangement | Underground; self and mutual multiple-cable entries. | Stated — §2.1 and Fig. 1. |
| Earth structure | Homogeneous earth below air. | Stated — §2. |
| Conductor and insulation geometry | Circular insulated cable positions; self replaces ``r`` by total cable radius including jacket and sets equal depths. | Stated — below (3), p. 2. |
| Constitutive and field assumptions | Linear, isotropic, nonmagnetic ground; quasi-TEM parent. | Stated — §2.1. |
| Conventions | ``j\omega`` as printed; positive burial depths; per-unit-length impedance. | Equation-implied — (1)–(8). |

**Expression.** Proposed small-argument impedance, equation (8), article p. 2.

```math
Z_{g(m,n)}=\frac{j\omega\mu_0}{2\pi}\left[
-\ln\left(\frac{\gamma_1d}{2}\right)-\gamma_E
+\frac{\gamma_1-\gamma_0}{\gamma_0+\gamma_1}
e^{-(h_m+h_n)\gamma_1}\left(\frac{2}{4+\gamma_1^2r^2}\right)
\right],
\qquad\text{(8)}
```

```math
d=\sqrt{(h_m-h_n)^2+r^2},\qquad \gamma_E=0.5772\ldots,
```

with ``\gamma_0,\gamma_1`` as above.

**Approximation.** The only additional operation on parent equation (1) is ``K_0(z)\approx-\ln(z/2)-\gamma_E`` for ``0<z\ll1``. All other terms are retained unchanged.

**Limitations.** The small-argument condition controls the formula, and the paper's authors restrict dependable wide-configuration use to roughly 1 MHz. Homogeneous nonmagnetic earth, parallel infinite cables and quasi-TEM remain. Branches of the complex logarithm and square roots are unstated.

**Reference.** [DeConti2023b](@cite). A. De Conti, N. Duarte, R. Alipio, and O. E. Leal, “Small-argument analytical expressions for the calculation of the ground-return impedance and admittance of underground cables,” *Electric Power Systems Research* 220, 109299 (2023), equations (1)–(3), (7)–(8), article p. 2.

**Transcription source.** Original publication page image. The leading minus sign, logarithm argument, Euler constant, reflection ratio, exponential and rational term were visually checked.

## Source transcription

Parent equation (1) is transcribed in the 2023 closed-form companion record. This record retains the source expansion ``K_0(z)\approx-\ln(z/2)-\gamma_E`` and its declared domain without adding higher terms.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\gamma`` in (7) | ``\gamma_E`` | Euler–Mascheroni constant | renamed only to avoid collision with propagation constants |
| ``\gamma_0,\gamma_1`` | unchanged | air and earth bulk constants | ``\mathrm m^{-1}`` |
| ``h_m,h_n,r,d`` | unchanged | depths, separation/self radius, direct distance | m |

The single ``\gamma\mapsto\gamma_E`` renaming is exact and reversible.

## Evidence and approximation sources

Equations (1) and (7) are the complete parent and mathematical operation; (8) is their direct substitution. The empirical range is reported separately in §§3–5.

## Limitations and discrepancies

-
- The paper's use of the same glyph family for Euler's constant and propagation constants is disambiguated only in the notation map, not mathematically changed.
