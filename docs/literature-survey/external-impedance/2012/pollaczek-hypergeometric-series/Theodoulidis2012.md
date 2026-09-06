# Theodoulidis uniformly convergent hypergeometric series for Pollaczek's integral

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Depths ``h_1,h_2``, horizontal separation ``x``; radius used for the self ``x``; no insulation region. |
| Calculated quantities | Exact confluent-hypergeometric series for the Pollaczek interface integral in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous earth half-space. |
| Model and approximation | The infinite series is an exact representation of (2). It comes from the exact finite-range form (22) and Taylor expansions about ``t=1``; only finite truncation introduces approximation. The source's remainder estimate gives the stated 10-, 21-, and 40-term guarantees over ``t\in[0,1]``. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked for final series and parent |

**Description.** Confluent-hypergeometric infinite series for exact evaluation of the homogeneous-earth Pollaczek interface integral. The author gives fixed truncation guarantees independent of geometry and frequency for the series' mathematical error.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` is retained; the parent is explicitly restricted to the TEM/transmission-line limit. | Stated — p. 807; Equation-implied — (1)–(2). |
| Air propagation constant ``γ_air`` | Not independently retained; the air interface appears through image geometry. | Equation-implied — (1)–(2), p. 807. |
| Earth propagation constant ``γ_earth`` | ``k=\sqrt{j\omega\mu_0\sigma}=(1+j)/\delta``; optional Sunde replacement ``k=\sqrt{j\omega\mu_0\sigma+\omega^2\mu_0\varepsilon_0\varepsilon_r}``. | Stated — pp. 807–808. |
| Earth permittivity and displacement current | Classical formula neglects it; optional Sunde formula retains ``\varepsilon_0\varepsilon_r``. | Stated — p. 808. |
| Range of validity | Infinite series valid for all parameter ranges of the parent integral. Truncation at ``n\le20`` is stated to guarantee relative error below order ``10^{-7}``; 10 terms give ``10^{-4}``, and 40 terms ``10^{-11}``. These are series-error statements, not physical-model validation. | Stated — pp. 807 and 811. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — (1), p. 807. |
| Arrangement | Two parallel buried conductors; mutual impedance and prescribed self substitution. Mixed placement is mentioned but not printed. | Stated — Fig. 1 and pp. 807–808. |
| Earth structure | Homogeneous earth half-space. | Stated — §II, p. 807. |
| Conductor and insulation geometry | Depths ``h_1,h_2``, horizontal separation ``x``; radius used for the self ``x``; no insulation region. | Stated — Fig. 1 and p. 808. |
| Constitutive and field assumptions | Linear isotropic earth, TEM Pollaczek parent; exact mathematical evaluation of its interface integral only. | Stated — pp. 807 and 812–813. |
| Conventions | ``j`` imaginary unit; ``\omega=2\pi f``; positive-downward depths; impedance per unit length; root ``k=(1+j)/\delta``. | Stated — (1)–(2) and Fig. 1, p. 807. |

**Expression.** Third exact series for ``I_{\mathrm{Pollaczek}}``, source equation (6), inserted into (3) and (1).

```math
Z(j\omega)=\frac{j\omega\mu_0}{2\pi}
\left[K_0(kr)-K_0(kR)+2J_{\mathrm{Pollaczek}}\right],
\qquad\text{(1)}
```

```math
J_{\mathrm{Pollaczek}}
=\left(\frac{H}{R}\right)^2K_0(kR)
+\frac{1}{kR}\left[2\left(\frac{H}{R}\right)^2-1\right]K_1(kR)
-\frac{1}{k^2}I_{\mathrm{Pollaczek}},
\qquad\text{(3)}
```

```math
\begin{aligned}
I_{\mathrm{Pollaczek}}
={}&\frac{H^2-x^2}{R^4}e^{-kH}(1+kH)\\
&+2\sqrt{2}\,\frac{kxH}{R^3}
\sum_{n=0}^{\infty}
\frac{(2n-3)!!}{2^{2n}n!}
\left(1-\frac{H}{R}\right)^{n+1/2}\\
&\quad\times\left\{
{}_1F_1\!\left[1,n+\frac32,-k(R-H)\right]
\left[1-\frac{kR}{2}\frac{2n-1}{2n+1}\right]-1
\right\}.
\end{aligned}
\qquad\text{(6)}
```

Definitions:

```math
J_{\mathrm{Pollaczek}}=\int_0^\infty
\frac{e^{-H\sqrt{\lambda^2+k^2}}}{\lambda+\sqrt{\lambda^2+k^2}}
\cos(\lambda x)\,d\lambda,
\qquad\text{(2)}
```

```math
r=\sqrt{x^2+(h_1-h_2)^2},\qquad
R=\sqrt{x^2+H^2},\qquad H=h_1+h_2,
\qquad k=\sqrt{j\omega\mu_0\sigma}=\frac{1+j}{\delta}.
```

``{}_1F_1`` is the confluent hypergeometric function. The source defines ``(-3)!!=-1`` and ``(-1)!!=1`` for the first two series coefficients.

**Approximation.** The infinite series is an exact representation of (2). It comes from the exact finite-range form (22) and Taylor expansions about ``t=1``; only finite truncation introduces approximation. The source's remainder estimate gives the stated 10-, 21-, and 40-term guarantees over ``t\in[0,1]``.

**Limitations.** Uniform mathematical convergence is not a full-wave validity claim. The parent Pollaczek formula can produce negative mutual resistance at sufficiently high frequency; the source retains this as a defect of the physical formula, not numerical integration. The paper prints no mixed-arrangement series. The 21-term statement uses indices ``n\le20``; it must not be paraphrased as an arbitrary 20-term convention.

**Reference.** [Theodoulidis2012](@cite), equations (1)–(3), (6), and (22)–(28), printed pp. 807–811 (PDF pages 2–6), with physical-validity discussion on pp. 812–813.

**Transcription source.** Original IEEE publication. The page images were used to verify the first term, all powers, ``2\sqrt2`` coefficient, ``(2n-3)!!`` convention, hypergeometric arguments and nested bracket in (6). The Markdown conversion severely corrupts this equation and was not equation authority.

## Source transcription

The exact finite-interval parent selected by the author is Wedepohl and Wilcox's form

```math
I_{\mathrm{Pollaczek}}
=\frac{H^2-x^2}{R^4}e^{-kH}(1+kH)
+\frac{k^2xH}{R^2}\int_{H/R}^{1}
\left(2\sqrt{1-t^2}-\frac{1}{\sqrt{1-t^2}}\right)e^{-tkR}\,dt.
\qquad\text{(22)}
```

The source expands the algebraic factors about ``t=1`` in (23)–(24), evaluates the resulting integrals through an incomplete-gamma identity (26), maps it to ``{}_1F_1`` in (27), and obtains the single series (28), which combined with (22) is (6). These operations are retained as source basis and do not alter the printed final series.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z(j\omega)`` | unchanged | Earth-return impedance per unit length | ``\Omega/\mathrm m`` |
| ``J_{\mathrm{Pollaczek}}`` | unchanged | Parent interface integral | ``\mathrm m^{-1}`` |
| ``I_{\mathrm{Pollaczek}}`` | unchanged | Decomposition auxiliary | ``\mathrm m^{-2}`` |
| ``x,h_1,h_2,H,r,R`` | unchanged | Separation, depths, depth sum and distances | metres |
| ``\lambda,k`` | unchanged | Spectral variable and earth constant | ``\mathrm m^{-1}`` |
| ``{}_1F_1`` | unchanged | Confluent hypergeometric function | dimensionless |
| ``\sigma,\varepsilon_0\varepsilon_r,\mu_0`` | unchanged | Earth constitutive quantities | SI |
| ``\delta`` | unchanged | Skin depth | metres |

No notation was renamed.

## Evidence and approximation sources

- Parent, geometry and TEM restriction: p. 807.
- Series (6), hypergeometric definition and shared variables: p. 808.
- Finite-range expression and Taylor/incomplete-gamma derivation: pp. 810–811, (22)–(28).
- Double-factorial terminal definitions and truncation guarantees: p. 811.
- Parent physical limitation and critical-frequency discussion: pp. 812–813.

## Limitations and discrepancies

- The existing BibTeX entry omits the source DOI and is copied exactly as found.
- The paper's phrase “exact solution” applies to the integral evaluation. It does not erase the TEM model restriction or the published negative-resistance finding.
