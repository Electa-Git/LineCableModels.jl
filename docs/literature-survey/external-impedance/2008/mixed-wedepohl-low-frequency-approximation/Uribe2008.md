# Uribe's secondary transcription of the Wedepohl mixed approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead and buried line axes. |
| Calculated quantities | Wedepohl-attributed low-frequency mixed approximation and Uribe's normalized form |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Low-frequency retained-term approximation of the Wedepohl–Wilcox series; exact retained order and unambiguous frequency bound unresolved. |
| Main source | L. M. Wedepohl and D. J. Wilcox (1973), secondary mixed-form attribution through F. A. Uribe (2008) |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Secondary mixed-form transcription checked against page images against Uribe; exact original mixed locator unresolved |

**Description.** Low-frequency Wedepohl-attributed approximation, as printed by Uribe, for per-unit-length mutual ground impedance between an overhead and a buried power line.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed to the quasi-TEM Pollaczek parent used by Uribe. | Stated — introduction to (1a), p. 198, and §V.C. |
| Air propagation constant ``γ_air`` | Neglected. | Equation-implied — (6e). |
| Earth propagation constant ``γ_earth`` | Represented by complex depth ``p=1/\sqrt{j\omega\mu_0\sigma}``. | Stated — (1b), used in (6e). |
| Earth permittivity and displacement current | Neglected. | Equation-implied — ``p`` contains ``\sigma`` only. |
| Range of validity | Uribe says the retained terms apply “up to quite frequencies” as printed and then assesses numerical error over selected curves; no unambiguous quantitative bound is supplied by the restatement. | Stated — §V.C and §VI. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Equation-implied — (6e). |
| Arrangement | Mixed mutual overhead/buried interaction. | Stated — §V.C and following application sentence. |
| Earth structure | Homogeneous soil half-space below air. | Inherited/stated — Fig. 1. |
| Conductor and insulation geometry | Parallel filamentary line axes; ``d`` is conductor separation and ``h_1+h_2`` the summed height/depth. No insulation term. | Equation-implied — (6e). |
| Constitutive and field assumptions | Linear homogeneous nonmagnetic soil; low-frequency truncation of the Wedepohl–Wilcox series as described by Uribe. | Stated — §V.C. |
| Conventions | ``j=\sqrt{-1}``; positive heights/depths; ``\gamma`` in (6e) is Euler's constant; ``\underline Z_G=Z_G\pi/(\omega\mu_0)``. | Stated — immediately after (6e) and (5a). |

**Expression.** Uribe prints

```math
Z_{G-W}=\frac{j\omega\mu_0}{2\pi}
\left[
-\log\!\left(\frac{\gamma_E d}{2p}\right)
+0.5-\frac{4(h_1+h_2)}{3p}
\right],
\tag{6e}
```

where the source's ``\gamma`` is Euler's constant and ``p=1/\sqrt{j\omega\mu_0\sigma}``.

**Approximation.** Uribe states that Wedepohl–Wilcox's full series was cumbersome/slow in some ranges and that (6e) retains only the terms relevant up to the source's qualitative frequency condition. Exact retained series order and a remainder bound are not reproduced.

**Limitations.** Original mixed-form equation was not independently located in the credited 1973 source. The phrase “up to quite frequencies” is retained as an unclear source statement, not silently interpreted as a numerical range.

**Reference.** [Uribe2008](@cite), equations (6e)–(6f), printed p. 201, as the inspected secondary transcription.

**Transcription source.** Secondary for the mixed Wedepohl form; original publication by Uribe. Logarithm argument, Euler-constant statement, ``0.5``, ``4/3`` term, and normalized factors were image-verified.

## Source transcription

Uribe's normalized form is

```math
\underline Z_{G-W}=j\left[
-\log\!\left(
\gamma_E\sqrt j\,\xi\sqrt{\eta^2+(\zeta-1)^2}/2
\right)
+0.5-\frac{2}{3}\sqrt j\,\xi(\zeta+1)
\right].
\tag{6f}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{G-W}`` | unchanged | Wedepohl-attributed mixed approximation | ``\Omega/\mathrm m`` |
| ``p`` | unchanged | complex depth | m |
| ``d,h_1,h_2`` | unchanged | direct distance, overhead height, burial depth | m |
| source ``\gamma`` | ``\gamma_E`` | Euler's constant in (6e)/(6f) | dimensionless; renamed only to avoid collision with propagation constants |

The sole display renaming is ``\gamma_E`` for the source's Euler constant.

## Evidence and approximation sources

Uribe attributes the parent series to Wedepohl–Wilcox [7] and prints the mixed low-frequency restatement. The inspected 1973 paper supports the historical parent series but does not verify these exact mixed equations.

## Limitations and discrepancies

- Original mixed equation locator remains unresolved.
- The source's qualitative “quite frequencies” wording is ambiguous.
- No global validity/error bound is supplied.
