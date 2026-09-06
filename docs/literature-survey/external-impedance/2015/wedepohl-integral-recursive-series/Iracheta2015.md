# Iracheta-Cortez recursive series for the Wedepohl integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent; mutual uses cable depths ``h_1,h_2`` and horizontal spacing ``x``; self substitutes conductor radius ``R`` for ``x`` and ``h_1=h_2``. |
| Calculated quantities | Recursive infinite-series evaluator for the finite integral in the Wedepohl–Wilcox decomposition of buried-cable earth-return impedance |
| Earth structure | Homogeneous soil half-space below air. |
| Model and approximation | The exponential in the finite integral is expanded in its everywhere-convergent power series, accumulated recursively, and truncated by tolerance. The source recommends the series only for ``\|D/p\|\le30`` and approximates ``I_w=0`` beyond that threshold. |
| Main source | R. Iracheta-Cortez (2015) for the recursion; Pollaczek (1926) and Wedepohl–Wilcox (1973) for the physical kernel and decomposition |
| Citation key(s) | `:Iracheta2015` |
| Evidence status | Original publication page images checked |

**Description.** Numerically recursive evaluation of the residual finite integral ``I_w`` in the Wedepohl–Wilcox representation of Pollaczek's homogeneous-earth buried-cable impedance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected by the inherited Pollaczek model. | Stated through adoption of Pollaczek's integral, pp. 34–35. |
| Air propagation constant ``γ_air`` | Not retained; air is nonconducting with ``\mu_1\simeq\mu_0`` and the reduced complex-depth form is used. | Stated — Fig. 1 and (1)–(3). |
| Earth propagation constant ``γ_earth`` | Represented by complex depth ``p=\sqrt{\rho_{soil}/(j\omega\mu_0)}``, so ``1/p`` is the conduction-only earth constant. | Stated — below (3a), p. 35. |
| Earth permittivity and displacement current | Neglected. | Stated — homogeneous soil model and complex-depth definition, pp. 34–35. |
| Range of validity | Series tested for ``|D/p|\le30``; source sets its residual to zero above 30. Convergence is stated for ``|t|<1``; tested 1 Hz–1 MHz and broad geometry/conductivity ranges. | Stated — pp. 37–38. |
| Earth permeability ``μ_earth`` | ``\mu_{ground}=\mu_0``. | Stated — below (1), p. 35. |
| Arrangement | Underground; self and mutual parallel cables. | Stated — Fig. 1 and below (3b). |
| Earth structure | Homogeneous soil half-space below air. | Stated — Fig. 1 and surrounding text. |
| Conductor and insulation geometry | Filament parent; mutual uses cable depths ``h_1,h_2`` and horizontal spacing ``x``; self substitutes conductor radius ``R`` for ``x`` and ``h_1=h_2``. | Stated — below (3b), p. 35. |
| Constitutive and field assumptions | Linear homogeneous nonmagnetic soil, classical quasi-static Pollaczek kernel; recursion changes only numerical evaluation. | Stated/inherited — introduction and §Earth return impedance. |
| Conventions | ``h=(h_1+h_2)/2``, ``D=\sqrt{x^2+(h_1+h_2)^2}``, ``t=2h/D``; per-unit-length impedance. | Stated — (3)–(4), Fig. 1. |

**Expression.** The source embeds its proposed recursive residual in

```math
Z_T=\frac{j\omega\mu_0}{2\pi}
\left[K_0(d/p)-K_0(D/p)+J\right],
\qquad\text{(3a)}
```

```math
J=\frac{4h^2}{D^2}K_0(D/p)
+\frac{(4h^2-2x^2)p}{D^3}
\left[K_1(D/p)-(2h+p)\frac{e^{-2h/p}}{D}\right]-I_w,
\qquad\text{(3c)}
```

with

```math
I_w\approx S_{w,\infty}=-\frac{2h|x|}{D^2}
\left[2S_{w1,\infty}-S_{w2,\infty}\right],
\qquad\text{(6a)}
```

```math
S_{w1,n}=DP_{2n}A_{2n}+DP_{2n+1}A_{2n+1}+S_{w1,n-1},
\qquad
S_{w2,n}=DP_{2n}B_{2n}+DP_{2n+1}B_{2n+1}+S_{w2,n-1},
\qquad\text{(6j,6k)}
```

```math
DP_{2n}=\frac{(-D/p)DP_{2n-1}}{2n},\qquad
DP_{2n+1}=\frac{(-D/p)DP_{2n}}{2n+1},
\qquad\text{(6l)}
```

```math
A_{2n}=\frac{(2h)^{2n-1}|x|^3/D^{2n+2}+(2n-1)A_{2n-2}}{2n+2},
```

```math
A_{2n+1}=\frac{(2h)^{2n-1}|x|/D^{2n}+(2n-1)A_{2n-1}}{2n},
\qquad n=1,2,3,\ldots,
\qquad\text{(6m)}
```

```math
B_{2n}=\frac{(2h)^{2n}|x|^3/D^{2n+3}+(2n)B_{2n-2}}{2n+3},
```

```math
B_{2n+1}=\frac{(2h)^{2n-1}|x|/D^{2n+1}+(2n)B_{2n-1}}{2n+1},
\qquad n=1,2,3,\ldots.
\qquad\text{(6n)}
```

The initial sums are ``S_{w1,0}=A_0+DP_1A_1`` and ``S_{w2,0}=B_0+DP_1B_1``, with ``DP_1=-D/p``. The stopping ratios are printed in (7).

**Approximation.** The exponential in the finite integral is expanded in its everywhere-convergent power series, accumulated recursively, and truncated by tolerance. The source recommends the series only for ``|D/p|\le30`` and approximates ``I_w=0`` beyond that threshold.

**Limitations.** This is an evaluator for an inherited physical kernel, not new electromagnetic physics. It is conduction-only, homogeneous-earth and filamentary. The source's printed initialization contains a label discrepancy described below, so a faithful transcription is not directly executable without resolving that token from another authoritative witness.

**Reference.** [Iracheta2015](@cite), equations (3a)–(8), printed pp. 35–37 (PDF pages 2–4).

**Transcription source.** Original journal page images. Recursion indices, powers, denominators, absolute values and the two ``DP`` recurrences were visually verified; the duplicated ``A_0`` label was not repaired.

## Source transcription

The series being accumulated is

```math
I_{w1}=\sum_{n=0}^{\infty}\left(-\frac Dp\right)^n\frac1{n!}
\underbrace{\int_{2h/D}^{1}t^n\sqrt{1-t^2}\,dt}_{A_n},
```

```math
I_{w2}=\sum_{n=0}^{\infty}\left(-\frac Dp\right)^n\frac1{n!}
\underbrace{\int_{2h/D}^{1}\frac{t^n}{\sqrt{1-t^2}}\,dt}_{B_n}.
\qquad\text{(5a,5b)}
```

The source prints ``A_0`` both for the integral of ``\sqrt{1-t^2}`` in (6f) and, again, for ``\int t\sqrt{1-t^2}dt=|x/D|^3/3`` in (6g), even though the initialization and flowchart require ``A_0,A_1``.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_T`` | unchanged | total earth-return self/mutual term | ``\Omega/\mathrm m`` |
| ``p`` | unchanged | complex depth | m |
| ``d,D`` | unchanged | direct and image distances | m |
| ``I_w,S_w`` | unchanged | finite Wedepohl integral and its recursive approximation | dimensionless |
| ``A_n,B_n`` | unchanged | definite moment integrals | dimensionless |
| ``DP_n`` | unchanged | recursively generated ``(-D/p)^n/n!`` factor | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Equations (3a)–(3d) are explicitly credited to Pollaczek and Wedepohl–Wilcox. Equations (6a)–(6n) are the paper's original implementation contribution and are therefore retained as an approximation/evaluator formula, without asserting a new physical earth-return model.

## Limitations and discrepancies

- Equation (6g) visibly labels the second initial moment ``A_0`` although (6d), the recurrence and Fig. 2 require ``A_1``. The corpus preserves the printed label and does not silently repair it.
