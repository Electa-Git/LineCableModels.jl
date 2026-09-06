# Theodoulidis second exact series for Pollaczek's integral

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary mutual geometry at depths ``h_1,h_2`` and separation ``x``; source-prescribed radius substitution for self; no insulation term. |
| Calculated quantities | Exact single-modified-Bessel series for the Pollaczek interface integral in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous half-space. |
| Model and approximation | Not an analytical approximation to the parent integral: the infinite series follows from the exact Bessel identity used in (14)–(16). Any finite truncation is numerical. The source supplies no universal fixed term count for (5); it says convergence is poor for large separation and good for ``x<H``. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked |

**Description.** Exact convergent series containing one modified Bessel function per term for the interface integral in the homogeneous-earth Pollaczek impedance of parallel buried conductors.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` occurs; the selected parent is explicitly confined to the TEM/transmission-line limit. | Stated — p. 807; Equation-implied — (1)–(2). |
| Air propagation constant ``γ_air`` | Not retained independently; air enters through the interface/image geometry. | Equation-implied — (1)–(2), p. 807. |
| Earth propagation constant ``γ_earth`` | ``k=\sqrt{j\omega\mu_0\sigma}=(1+j)/\delta``; optional Sunde formula ``k=\sqrt{j\omega\mu_0\sigma+\omega^2\mu_0\varepsilon_0\varepsilon_r}``. | Stated — pp. 807–808. |
| Earth permittivity and displacement current | Omitted in classical Pollaczek ``k`` and optionally retained through the separately stated Sunde replacement. | Stated — p. 808. |
| Range of validity | The infinite series is mathematically convergent, works well for ``x<H``, and is comparatively slow for large conductor separation. Its exactness does not remove the parent formula's high-frequency physical limitation. | Stated — pp. 809 and 812–813. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — (1) and §II, p. 807. |
| Arrangement | Parallel underground conductors; printed mutual term and self substitution ``h_1=h_2`` with ``x`` equal to conductor radius. | Stated — Fig. 1 and pp. 807–808. |
| Earth structure | Homogeneous half-space. | Stated — §II and Fig. 1, p. 807. |
| Conductor and insulation geometry | Filamentary mutual geometry at depths ``h_1,h_2`` and separation ``x``; source-prescribed radius substitution for self; no insulation term. | Stated — pp. 807–808. |
| Constitutive and field assumptions | Linear isotropic earth and TEM Pollaczek parent. Internal conductor effects and insulation are excluded. | Stated — pp. 807 and 812–813. |
| Conventions | ``j`` imaginary unit; ``\omega=2\pi f``; depths positive downward; impedance per unit length; classical branch ``k=(1+j)/\delta``. | Stated — Fig. 1 and (1)–(2), p. 807. |

**Expression.** Second exact series for ``I_{\mathrm{Pollaczek}}`` used by equations (1)–(3), source equation (5), printed p. 808.

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
I_{\mathrm{Pollaczek}}
=-k^2\sqrt{\frac{2}{\pi}}\,kH
\sum_{n=1}^{\infty}
\frac{1}{2^{2n}n!(2n-1)}
\frac{(kx)^{2n}}{(kR)^{n+3/2}}
K_{n+3/2}(kR).
\qquad\text{(5)}
```

Definitions retained from (2) and below (6):

```math
\begin{aligned}
r&=\sqrt{x^2+(h_1-h_2)^2} \\
R&=\sqrt{x^2+H^2} \\
H&=h_1+h_2 \\
u&=\sqrt{\lambda^2+k^2},
\end{aligned}
```

```math
\begin{aligned}
J_{\mathrm{Pollaczek}}&=\int_0^\infty
\frac{e^{-Hu}}{\lambda+u}\cos(\lambda x)\,d\lambda \\
k&=\sqrt{j\omega\mu_0\sigma}=\frac{1+j}{\delta}.
\end{aligned}\qquad\text{(2)}
```

``K_\nu`` is the modified Bessel function of the second kind. The self formula sets ``h_1=h_2`` and ``x`` to conductor radius.

**Approximation.** Not an analytical approximation to the parent integral: the infinite series follows from the exact Bessel identity used in (14)–(16). Any finite truncation is numerical. The source supplies no universal fixed term count for (5); it says convergence is poor for large separation and good for ``x<H``.

**Limitations.** Relative to series (4), this expression is less convenient for large separation. Extreme order or argument can require scaled Bessel evaluation. Exact integral evaluation does not fix negative mutual resistance of the TEM Pollaczek parent at sufficiently high frequency. The mentioned mixed overhead/underground extension is not printed.

**Reference.** [Theodoulidis2012](@cite), equations (1)–(3), (5), and (14)–(16), printed pp. 807–810 (PDF pages 2–5), with validity discussion on pp. 812–813.

**Transcription source.** Original IEEE page images. The leading ``-k^2``, ``\sqrt{2/\pi}\,kH``, lower summation bound ``n=1``, factorial factor, powers, and Bessel order in (5) were checked visually. The equation-damaged Markdown conversion was not used to choose tokens.

## Source transcription

The source presents (5) as the second of three independent infinite-series evaluations of the same ``I_{\mathrm{Pollaczek}}``. It is derived by inserting the exact cosine identity

```math
\cos z=\sum_{n=0}^{\infty}\frac{(z/2)^n}{n!(2n-1)}J_n(z)
\qquad\text{(14)}
```

into (9a), followed by the integral identity used in (15), yielding source (16), which is identical to displayed (5). The source operation is recorded without changing the series.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z(j\omega)`` | unchanged | Earth-return impedance per unit length | ``\Omega/\mathrm m`` |
| ``J_{\mathrm{Pollaczek}}`` | unchanged | Parent interface integral | ``\mathrm m^{-1}`` |
| ``I_{\mathrm{Pollaczek}}`` | unchanged | Decomposition auxiliary | ``\mathrm m^{-2}`` |
| ``x,h_1,h_2,H,r,R`` | unchanged | Separation, depths, depth sum, direct/image distances | metres |
| ``\lambda,u,k`` | unchanged | Spectral variable, root and earth constant | ``\mathrm m^{-1}`` |
| ``K_\nu`` | unchanged | Modified Bessel function of second kind | order ``\nu`` |
| ``\sigma,\varepsilon_0\varepsilon_r,\mu_0`` | unchanged | Earth constitutive quantities | SI |
| ``\delta`` | unchanged | Skin depth | metres |

No notation was renamed.

## Evidence and approximation sources

- Parent model and geometry: p. 807, Fig. 1 and (1)–(2).
- Optional Sunde constitutive replacement and self prescription: p. 808.
- Exact series (5): p. 808; cosine identity and derivation: pp. 809–810, (14)–(16).
- Comparative convergence statement: p. 809.
- Parent-formula physical validity finding: pp. 812–813.

## Limitations and discrepancies

- The LCM identity year 2015 conflicts with the inspected 2012 publication; no identifier change is made.
- The copied existing bibliography entry lacks the printed DOI; metadata is preserved rather than repaired.
- The source recommends scaled Bessel functions for numerical overflow/underflow, but no scaled or rearranged expression is substituted here.
