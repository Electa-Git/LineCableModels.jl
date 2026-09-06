# Knight continuous fitted solid-round internal-impedance approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Solid homogeneous circular cylinder of diameter ``d=2r``. |
| Calculated quantities | Continuous closed approximations for p.u.l. AC resistance and internal inductance of an isolated solid round conductor |
| Earth structure | None. |
| Model and approximation | TED-ML and PACAML are empirical modified-Lorentzian corrections of functions constrained to both exact asymptotes. Their reported errors refer to comparison with the report's Kelvin-function calculation, not independent measurements. |
| Main source | David W. Knight (version 2.08.1, 2016) |
| Citation key(s) | `:Knight2016` |
| Evidence status | Author-report page image verified |

**Description.** Bessel-free, doubly asymptotically correct fitted functions spanning the DC and strong-skin-effect limits for a uniform solid cylindrical conductor.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally invariant isolated-conductor internal problem. | Stated — §§1, 3. |
| Air propagation constant ``γ_air`` | Not applicable. | Scope. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Scope. |
| Earth permittivity and displacement current | Good-conductor approximation; metal displacement current excluded. | Stated — discussion after (3.2). |
| Range of validity | Report claims no upper-frequency restriction for the continuous fits; empirical maximum errors are ±0.09% for resistance and ±0.016% for inductance against its exact reference. | Stated — abstract, pp. 31, 48. |
| Earth permeability ``μ_earth`` | Not applicable; conductor absolute permeability ``\mu`` enters skin depth. | Stated — §3. |
| Arrangement | Isolated single conductor; proximity absent. | Stated — title/scope. |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Solid homogeneous circular cylinder of diameter ``d=2r``. | Stated — §§1, 3. |
| Constitutive and field assumptions | Uniform linear good conductor; fitted to the Kelvin/Bessel exact cylindrical skin-effect solution. | Stated — §§3, 5, 16. |
| Conventions | ``Z_i=R_{ac}+j2\pi fL_i``; ``\delta_i=\sqrt{\rho/(\pi f\mu)}``. | Stated — §§1, 3. |

**Expression.** For resistance,

```math
\Xi=\frac{R_{ac}}{R_{dc}}=
\frac{r^2}{(2r\delta_i'-\delta_i'^2)(1+y)},\quad
\delta_i'=\delta_i[1-e^{-r/\delta_i}],\quad
\delta_i=\sqrt{\frac{\rho}{\pi f\mu}},
```

```math
y=\frac{0.189774}{\{1+0.272481[z^{1.82938}-z^{-0.99457}]^2\}^{1.0941}},qquad
z=0.62006\frac r{\delta_i}.
\qquad\text{(Rac-TED-ML)}
```

For internal inductance, set ``q=d/(\delta_i\sqrt2)`` and

```math
\Theta_\infty=\frac4{q\sqrt2}\left[1+\frac{0.01209}{q+1}-\frac{0.63523}{q^2+1}+\frac{0.16476}{q^3+1}\right],
```

```math
\Theta_{da}=\Theta_\infty[1-e^{-\Theta_\infty^{-1.5819}}]^{1/1.5819},\quad
\Theta=\Theta_{da}(1-y),
```

```math
y=\frac{-0.198584}{\{1+0.25741[z^{1.2652}-z^{-0.39709}]^2\}^{2.62343}},\quad
z=0.38691q,\qquad \frac{L_i}{\ell}=\frac{\mu}{8\pi}\Theta.
\qquad\text{(Li-PACAML)}
```

**Approximation.** TED-ML and PACAML are empirical modified-Lorentzian corrections of functions constrained to both exact asymptotes. Their reported errors refer to comparison with the report's Kelvin-function calculation, not independent measurements.

**Limitations.** No proximity effect, hollow wall, magnetic nonlinearity, compositional layering, dielectric displacement or external return. The fitted decimal coefficients are part of the formulas and should not be re-rounded silently.

**Reference.** [Knight2016](@cite).  Knight, *Practical continuous functions for the internal impedance of solid cylindrical conductors*, version 2.08.1 (2016), DOI `10.13140/RG.2.1.3865.1284`, pp. 31 and 48.

**Transcription source.** Author-report page images. Every fitted coefficient, exponent, sign, skin-depth definition and claimed maximum error was visually verified.

## Source transcription

The report's exact Kelvin-function parent is retained by locator rather than re-attributed as Knight's original contribution. ``R_{dc}/\ell=\rho/(\pi r^2)`` converts ``\Xi`` to p.u.l. resistance.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\Xi`` | unchanged | AC/DC resistance factor | dimensionless |
| ``\Theta`` | unchanged | internal-inductance factor | dimensionless |
| ``\delta_i`` | unchanged | good-conductor skin depth | m |
| ``\rho,\mu`` | unchanged | conductor resistivity/permeability | ``\Omega\,m``, ``\mathrm{H/m}`` |

No notation was renamed.

## Evidence and approximation sources

The two boxed formulas and their empirical accuracy labels are treated as the report's contribution. The Bessel parent and classical asymptotes remain prior art.

## Limitations and discrepancies

- This is a versioned author technical report rather than a journal article; future versions may differ.
- The report calls its reference result “exact” while explicitly invoking a good-conductor approximation; that qualification is preserved.
