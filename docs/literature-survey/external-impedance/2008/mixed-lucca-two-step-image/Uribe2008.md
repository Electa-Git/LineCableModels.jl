# Uribe's secondary transcription of the Lucca mixed approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead and buried line axes. |
| Calculated quantities | Lucca-attributed two-step image approximation and Uribe's normalized form |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Image-theory construction followed by suppression of the oscillatory exponential factor; the source gives no remainder bound. |
| Main source | G. Lucca (1994), secondary attribution through F. A. Uribe (2008) |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Secondary transcription checked against Uribe's page images; Lucca original not equation-verified. |

**Description.** Two-step image approximation, attributed by Uribe to Lucca, for the per-unit-length mutual ground impedance between one overhead line and one buried line.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed to the quasi-TEM Pollaczek parent used by Uribe. | Stated — introduction to (1a), p. 198, and §V, p. 201. |
| Air propagation constant ``γ_air`` | Neglected. | Equation-implied — (6a) contains only the soil wave number. |
| Earth propagation constant ``γ_earth`` | ``k_e^2=-j\omega\mu_0\sigma`` and ``\gamma=jk_e``. | Stated — definitions below (6a). |
| Earth permittivity and displacement current | Neglected. | Equation-implied — soil definition uses ``\sigma`` only. |
| Range of validity | Uribe says the commonly used approximations lacked established accuracy ranges and tests this formula over selected Table I/II domains. No universal bound is assigned. | Stated — opening of §V and §VI. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Equation-implied — (6a). |
| Arrangement | Mixed mutual overhead/buried interaction. | Stated — §V.A. |
| Earth structure | Homogeneous soil half-space below air. | Inherited/stated — Fig. 1 and parent model. |
| Conductor and insulation geometry | Parallel filamentary axes; no insulation contribution. ``x`` is horizontal spacing and ``y_1,y_2`` are source coordinates. | Stated/equation-implied — Fig. 1 and (6a). |
| Constitutive and field assumptions | Linear homogeneous nonmagnetic soil; first step uses image theory, second suppresses the oscillatory exponential factor. | Stated — §V.A. |
| Conventions | ``j=\sqrt{-1}``; per-unit-length impedance; Uribe normalizes by ``\underline Z_G=Z_G\pi/(\omega\mu_0)``. | Stated — (5a), (6a)–(6b). |

**Expression.** Uribe prints

```math
Z_{G-L}=\frac{j\omega\mu_0}{2\pi}
\left\{
\ln\!\left(\frac{\overline R_{12}}{R_{12}}\right)
-\frac{2\overline y}{3\gamma^3}
\left[\frac{\overline y^2-3x^2}{\overline R_{12}^{6}}\right]
\right\},
\qquad\text{(6a)}
```

where

```math
R_{12}=\sqrt{x^2+(y_1-y_2)^2},\quad
k_e^2=-j\omega\mu_0\sigma,\quad
\overline R_{12}=\sqrt{\overline y^2+x^2},\quad
\gamma=jk_e,\quad
\overline y=y_1-y_2+2/\gamma.
```

**Approximation.** Uribe describes Lucca's first step as image theory following Wait–Spies and the second as suppression of the integrand's oscillatory exponential factor. No remainder bound is reproduced.

**Limitations.** Original Lucca equations and derivation were not inspected. Accuracy is parameter-dependent; Uribe's tests are not universal bounds. Source coordinates are retained without conversion to Uribe's ``h_1,h_2`` notation.

**Reference.** [Uribe2008](@cite), equations (6a)–(6b), printed p. 201, as the inspected secondary transcription.

**Transcription source.** Secondary for Lucca; original publication by Uribe. Bars, powers, signs, and definitions were visually checked on p. 201.

## Source transcription

Uribe's dimensionless form is

```math
\underline Z_{G-L}=j\left\{
\ln\!\sqrt{\frac{\lambda^2+\eta^2}{\eta^2+(\zeta-1)^2}}
+\frac{2\lambda(\lambda^2-3\eta^2)}
{3\sqrt j\,\xi^3(\eta^2+\lambda^2)^3}
\right\},
\qquad
\lambda=\zeta-1+\frac{2\sqrt j}{\xi}.
\qquad\text{(6b)}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{G-L}`` | unchanged | Lucca-attributed mixed approximation | ``\Omega/\mathrm m`` |
| ``R_{12},\overline R_{12}`` | unchanged | direct and complex-image distances | m |
| ``\gamma`` | unchanged | ``jk_e`` in this formula | ``\mathrm m^{-1}`` |
| ``\underline Z_{G-L}`` | unchanged | Uribe-normalized impedance | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Uribe attributes the formula to Lucca reference [2], *Mutual impedance between an overhead and a buried line with earth-return*, 1994. This record verifies Uribe's transcription and comparison only, not original priority or token identity.

## Limitations and discrepancies

- The source gives no source-independent applicability range or analytical remainder.
