# Uribe's secondary transcription of the CCITT mixed approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead and buried line axes. |
| Calculated quantities | CCITT-recommended mixed mutual ground-impedance approximation and Uribe's normalized form |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Source-labelled recommended approximation; derivation, retained order, and universal range unresolved. |
| Main source | CCITT via G. Lucca and H. W. Dommel–J. Sawada, secondary attribution through F. A. Uribe (2008) |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Secondary transcription checked against page images against Uribe; credited recommendation/report not equation-verified |

**Description.** CCITT-recommended approximation, as printed by Uribe, for per-unit-length mutual ground impedance between an overhead and a buried line.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed to the quasi-TEM parent used by Uribe. | Stated — introduction to (1a), p. 198, and §V.B. |
| Air propagation constant ``γ_air`` | Neglected. | Equation-implied — (6c). |
| Earth propagation constant ``γ_earth`` | ``k_e^2=-j\omega\mu_0\sigma``. | Inherited from the definition below (6a) and used in (6c). |
| Earth permittivity and displacement current | Neglected. | Equation-implied — ``k_e`` uses conductivity only. |
| Range of validity | Uribe states that accuracy ranges of practical approximations were undetermined and evaluates this formula over selected parameter curves. No universal bound is supplied. | Stated — opening of §V and §VI. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Equation-implied — (6c). |
| Arrangement | Mixed mutual overhead/buried interaction. | Stated — §V.B. |
| Earth structure | Homogeneous soil half-space below air. | Inherited/stated — Fig. 1. |
| Conductor and insulation geometry | Parallel filamentary line axes with direct distance ``R_{12}``; no insulation term. | Stated/equation-implied — (6c). |
| Constitutive and field assumptions | Linear homogeneous nonmagnetic soil; approximate recommended formula. Derivation/order are not reproduced. | Stated — §V.B. |
| Conventions | ``j=\sqrt{-1}``; per-unit-length; ``\underline Z_G=Z_G\pi/(\omega\mu_0)``. | Stated — (5a), (6c)–(6d). |

**Expression.** Uribe prints

```math
Z_{G-C}=\frac{j\omega\mu_0}{2\pi}
\left\{
\ln\!\left(\frac{1.851}{jk_eR_{12}}\right)
+\frac{2jk_e(y_1+y_2)}{3}
\right\}.
\qquad\text{(6c)}
```

Here ``R_{12}=\sqrt{x^2+(y_1-y_2)^2}`` and ``k_e^2=-j\omega\mu_0\sigma`` as defined with the adjacent comparison formulas.

**Numerical interpretation.** Evaluation uses the dimensional (6c), with ``jk_e=\sqrt{j\omega\mu_0\sigma}`` on the positive-real-part branch. Coordinates above ground are positive and those below ground negative. Thus ``y_1+y_2=h_a-h_g`` and ``R_{12}=\sqrt{x^2+(h_a+h_g)^2}`` for positive height and depth magnitudes ``h_a,h_g``. Reversing the conductor order leaves the result unchanged. The normalized (6d) is retained as source text, not used as an independent rescaling. Only mixed mutual terms use this approximation; same-medium terms retain the existing integral.

**Approximation.** Source-labelled approximate recommendation. Uribe reproduces no expansion parameter, retained order, discarded terms, or analytical remainder for the CCITT derivation; those dependencies remain unresolved.

**Limitations.** The original CCITT wording and Dommel–Sawada report equation were not inspected. Coordinate symbols are preserved as printed. Uribe's tested error curves do not establish a general range.

**Reference.** [Uribe2008](@cite), equations (6c)–(6d), printed p. 201, as the inspected secondary transcription.

**Transcription source.** Secondary for CCITT/Dommel–Sawada; original publication by Uribe. The constant ``1.851``, leading logarithm, ``2/3`` correction, and signs were image-verified.

## Source transcription

Uribe's normalized form is

```math
\underline Z_{G-C}=j\left\{
\ln\!\left[
\frac{1.851\sqrt j}{\xi\sqrt{\eta^2+(\zeta-1)^2}}
\right]
-\frac{2j\sqrt j\,\xi(1+\zeta)}{3}
\right\}.
\qquad\text{(6d)}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{G-C}`` | unchanged | CCITT-attributed mixed approximation | ``\Omega/\mathrm m`` |
| ``R_{12}`` | unchanged | direct conductor distance | m |
| ``k_e`` | unchanged | conductive-soil wave number | ``\mathrm m^{-1}`` |
| ``\underline Z_{G-C}`` | unchanged | Uribe-normalized impedance | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Uribe describes (6c) as the CCITT formula and cites Lucca [2] plus the Dommel–Sawada report [6]. With those originals uninspected, the corpus does not assign a more precise priority or derivation.

## Limitations and discrepancies

- The derivation and formal applicability conditions are not reproduced by Uribe.
