# De Conti–de Lima Padé approximation of Sunde's underground-cable impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular insulated cables represented externally by their cable coordinates; ``h_m,h_n`` are burial depths and ``r`` is horizontal separation (or the source's self-distance substitution). |
| Calculated quantities | Self and mutual per-unit-length ground-return impedance of underground cables; closed-form 1/1 Padé approximation of Sunde's residual integral |
| Earth structure | Homogeneous earth; the displayed Sunde parent has no upper-medium spectral term. |
| Model and approximation | Starting from Sunde's exact-in-model rearrangement (3)–(4), the source expands about ``t=1`` and replaces only ``e^{-t\gamma D}`` by the 1/1 Padé form ``[(2-\gamma D(t-1))/(2+\gamma D(t-1))]e^{-\gamma D}`` in (6). Integrating that rational replacement gives (7)–(10), which is inserted into (5). No additional term is discarded in the printed construction. |
| Main source | A. De Conti and A. C. S. de Lima (2024) |
| Citation key(s) | `:DeConti2024` |
| Evidence status | Original publication page images checked |

**Description.** Closed-form approximation of the self or mutual per-unit-length ground-return impedance between insulated cables buried in a homogeneous earth, obtained by applying a 1/1 Padé approximation to the remaining finite integral in Sunde's expression.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` appears; the parent is Sunde's quasi-TEM underground-cable expression. | Equation-implied — (1)–(5), p. 994; parent characterization, p. 993. |
| Air propagation constant ``γ_air`` | Not present in Sunde's displayed parent or the Padé approximation. | Equation-implied — (1)–(10), p. 994. |
| Earth propagation constant ``γ_earth`` | ``\gamma=\sqrt{j\omega\mu(\sigma+j\omega\epsilon)}``. | Stated — text below (1), p. 994. |
| Earth permittivity and displacement current | Retained through ``j\omega\epsilon`` in ``\gamma``; constant-parameter and a cited frequency-dependent soil model are both tested. | Stated — pp. 994–995. |
| Range of validity | Tested from 1 Hz to 10 MHz. NRMSE is below about 2.5% for ``H/D>0.5`` (approximately ``r/H<3`` for horizontally aligned cables); ``H/D>0.6`` (``r/H<4/3``) is the stated stricter negligible-error criterion up to 10 MHz. Both may be relaxed for lower frequency limits. | Stated — (11) and Figs. 2–5, pp. 994–996. |
| Earth permeability ``μ_earth`` | ``\mu`` is retained in ``\gamma`` while the impedance prefactor is printed as ``\mu_0``. | Stated/equation-implied — (1), (5), and definition below (1), p. 994. |
| Arrangement | Underground; self and mutual terms for cables at arbitrary burial depths and horizontal separation. | Stated — Fig. 1 and §II, p. 994; self-error statement, p. 996. |
| Earth structure | Homogeneous earth; the displayed Sunde parent has no upper-medium spectral term. | Stated — §II and (1)–(5), p. 994. |
| Conductor and insulation geometry | Circular insulated cables represented externally by their cable coordinates; ``h_m,h_n`` are burial depths and ``r`` is horizontal separation (or the source's self-distance substitution). | Stated — Fig. 1 and definitions below (1), p. 994. |
| Constitutive and field assumptions | Linear isotropic earth and quasi-TEM Sunde parent. Internal conductor and insulation fields are not included. | Stated — pp. 993–994. |
| Conventions | ``j\omega`` time-harmonic convention as printed; positive burial depths, horizontal separation ``r``, per-unit-length output in ``\Omega/\mathrm m``. | Equation-implied — Fig. 1 and (1)–(5), p. 994. |

**Expression.** The proposed impedance is equation (5) with the Padé-evaluated residual (7)–(10), printed p. 994.

```math
Z_g=\frac{j\omega\mu_0}{2\pi}\left\{
K_0(\gamma d)+\frac{H^2-r^2}{D^2}K_2(\gamma D)
-2\frac{H^2-r^2}{\gamma^2D^4}e^{-\gamma H}(1+\gamma H)
-\frac{2rH}{D^2}I_S
\right\},
\tag{5}
```

```math
I_S=(I_1+I_2+I_3)e^{-\gamma D},
\tag{7}
```

```math
I_1=\left(H-\frac{8}{\gamma}\right)\frac{r}{D^2},
\qquad
I_2=16\frac{2-\gamma D}{\gamma^2D^2}
\arctan\left(\frac{r}{H+D}\right),
\tag{8--9}
```

```math
I_3=-4\frac{8-8\gamma D+\gamma^2D^2}
{\gamma^2D^2\sqrt{1-\gamma D}}
\arctan\left(\frac{r\sqrt{1-\gamma D}}{H+D}\right).
\tag{10}
```

```math
d=\sqrt{(h_m-h_n)^2+r^2},\qquad
H=h_m+h_n,\qquad D=\sqrt{H^2+r^2},\qquad
\gamma=\sqrt{j\omega\mu(\sigma+j\omega\epsilon)}.
```

**Approximation.** Starting from Sunde's exact-in-model rearrangement (3)–(4), the source expands about ``t=1`` and replaces only ``e^{-t\gamma D}`` by the 1/1 Padé form ``[(2-\gamma D(t-1))/(2+\gamma D(t-1))]e^{-\gamma D}`` in (6). Integrating that rational replacement gives (7)–(10), which is inserted into (5). No additional term is discarded in the printed construction.

**Limitations.** Accuracy worsens as cables become widely separated relative to total burial depth; the numerical limits above are tested bounds, not universal proofs. The displayed parent omits the upper-medium propagation constant, so the approximation can only approximate that Sunde operator directly. The source notes that a more general Xue–Magalhães kernel includes ``\gamma_0`` and tests (5) against it separately. Branch choices for ``\sqrt{1-\gamma D}``, the arctangent, Bessel functions, and ``\gamma`` are not stated beside the formula. The source prints ``\mu_0`` in the impedance prefactor but arbitrary ``\mu`` in ``\gamma``; this is preserved, not harmonized.

**Reference.** [DeConti2024](@cite).  A. De Conti and A. C. S. de Lima, “A Closed-Form Approximation for Sunde's Ground-Return Expression Based on a Padé Approximant,” *IEEE Transactions on Electromagnetic Compatibility*, 66(3), 993–1000 (2024), DOI `10.1109/TEMC.2024.3377553`, equations (1)–(10), p. 994.

**Transcription source.** Original final-publication page image, printed p. 994. The four terms of (5), Padé numerator/denominator, exponential factor, coefficients, square root and arctangent arguments in (6)–(10) were visually checked. Text extraction was used only for navigation and surrounding accuracy statements.

## Source transcription

The parent source order is retained by recording the Sunde integral and its exact rearrangement:

```math
Z_g=\frac{j\omega\mu_0}{2\pi}[K_0(\gamma d)-K_0(\gamma D)+2J_S],
\qquad
J_S=\int_0^\infty
\frac{e^{-H\sqrt{\lambda^2+\gamma^2}}}
{\lambda+\sqrt{\lambda^2+\gamma^2}}
\cos(r\lambda)\,d\lambda,
\tag{1--2}
```

```math
J_S=\left(\frac HD\right)^2K_0(\gamma D)
+\frac1{\gamma D}\left[2\left(\frac HD\right)^2-1\right]K_1(\gamma D)
-\frac{H^2-r^2}{\gamma^2D^4}e^{-\gamma H}(1+\gamma H)
-\frac{rH}{D^2}I_S,
\tag{3}
```

```math
I_S=\int_{H/D}^1\left(2\sqrt{1-t^2}-\frac1{\sqrt{1-t^2}}\right)e^{-t\gamma D}\,dt.
\tag{4}
```

The approximation operation is

```math
e^{-t\gamma D}\approx
\frac{2-\gamma D(t-1)}{2+\gamma D(t-1)}e^{-\gamma D}.
\tag{6}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_g`` | unchanged | Ground-return impedance | ``\Omega/\mathrm m`` |
| ``h_m,h_n`` | unchanged | Burial depths | m |
| ``r`` | unchanged | Horizontal cable separation; source self substitution is not restated beside (1) | m |
| ``d,D,H`` | unchanged | Direct distance, image distance, and depth sum | m |
| ``\lambda`` | unchanged | Spectral variable | ``\mathrm m^{-1}`` |
| ``\gamma`` | unchanged | Earth propagation constant | ``\mathrm m^{-1}`` |
| ``J_S`` | unchanged | Sunde semi-infinite integral | dimensionless |
| ``I_S`` | unchanged | Residual finite integral / Padé approximation | dimensionless after accompanying geometric factors |
| ``K_\nu`` | unchanged | Modified Bessel function of second kind | dimensionless argument |
| ``t`` | unchanged | Finite-integral variable | dimensionless |

No notation was renamed.

## Evidence and approximation sources

The formula's complete parent-to-child chain is printed contiguously as (1)–(10) on p. 994. Accuracy was assessed at 101 logarithmically spaced frequencies over 1 Hz–10 MHz, ground resistivities 10–10,000 ``\Omega\,\mathrm m``, several depths and separations, and constant and cited frequency-dependent soil models. The stated ``H/D`` and ``r/H`` criteria appear on p. 996. These tests establish the reported domain only.

## Limitations and discrepancies

- The paper attributes equations (1)–(4) to Sunde and Theodoulidis. Only the Padé operation and its integrated result are the 2024 contribution.
- The source's more-general-kernel and admittance comparisons do not turn (5) into a full-wave or admittance formula; they are validation/application evidence.
