# Rallis DCIM approximation of the mixed Pollaczek integral

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filaments at overhead height ``h_1``, burial depth ``h_2`` and horizontal spacing ``x``. |
| Calculated quantities | Finite rational complex-image sum for mutual earth-return impedance between one overhead and one buried conductor |
| Earth structure | Homogeneous conductive half-space below air. |
| Model and approximation | GPOF fits the spectral factor over a finite sampling interval and (4.14), ``\int_0^\infty e^{-p\lambda}\cos(q\lambda)d\lambda=p/(p^2+q^2)``, integrates each term. The two-level sampling variant is an empirical refinement. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF representation; Pollaczek for the parent kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

**Description.** DCIM/GPOF approximation of the mixed overhead/buried Pollaczek integral, yielding a short sum of rational complex-image terms.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected by the Pollaczek parent. | Stated through (4.4). |
| Air propagation constant ``γ_air`` | Not retained. | Stated by parent reduction. |
| Earth propagation constant ``γ_earth`` | ``k^2=j\omega\mu_0\sigma``. | Stated — §4.2. |
| Earth permittivity and displacement current | Neglected. | Stated — definition of ``k``. |
| Range of validity | Demonstrated with 100 samples, 14 terms and ``T_0=28|k|``; two-level GPOF uses ``T_{02}=5|k|`` and ``T_{01}=100|k|``. No global bound is supplied. | Stated — pp. 75–76. |
| Earth permeability ``μ_earth`` | ``\mu_0``. | Stated — (4.3). |
| Arrangement | Mutual term between one overhead and one underground parallel conductor. | Stated — Fig. 4.2. |
| Earth structure | Homogeneous conductive half-space below air. | Stated — Fig. 4.2. |
| Conductor and insulation geometry | Filaments at overhead height ``h_1``, burial depth ``h_2`` and horizontal spacing ``x``. | Stated — Fig. 4.2. |
| Constitutive and field assumptions | Linear isotropic soil; quasi-static classical kernel; finite sum is an evaluator approximation. | Stated/inherited. |
| Conventions | ``j=\sqrt{-1}``, positive height/depth, per-unit-length mutual impedance. | Stated — (4.3)–(4.4). |

**Expression.** The mixed kernel is approximated by

```math
\frac{e^{-h_2\sqrt{\lambda^2+k^2}}}
{\lambda+\sqrt{\lambda^2+k^2}}
\approx\sum_{n=1}^{N}c_ne^{s_n\lambda},
\qquad\text{(4.11)}
```

```math
J_{ua}\approx2\sum_{n=1}^{N}c_n\frac{H_n}{H_n^2+x^2},
\qquad H_n=h_1-s_n,qquad n=1,2,\ldots,N,
\qquad\text{(4.15)}
```

and ``Z_{ua}=j\omega\mu_0J_{ua}/(2\pi)``.

**Approximation.** GPOF fits the spectral factor over a finite sampling interval and (4.14), ``\int_0^\infty e^{-p\lambda}\cos(q\lambda)d\lambda=p/(p^2+q^2)``, integrates each term. The two-level sampling variant is an empirical refinement.

**Limitations.** Mutual mixed geometry only. Coefficients must be regenerated from the fit; no error theorem or universal coefficient set is supplied. Displacement current, finite radii, layering and longitudinal propagation are absent.

**Reference.** [Rallis2013](@cite).  K. V. Rallis, doctoral thesis, 2012, DOI `10.12681/eadd/34633`, equations (4.3)–(4.4), (4.10)–(4.15), printed pp. 69, 75–76.

**Transcription source.** Original thesis page images. The fitted numerator, choice of ``\lambda`` rather than ``\sqrt{\lambda^2+k^2}``, rational denominator, factor two and ``H_n=h_1-s_n`` were visually verified.

## Source transcription

```math
J_{ua}=2\int_0^\infty
\frac{e^{-h_1\lambda-h_2\sqrt{\lambda^2+k^2}}}
{\lambda+\sqrt{\lambda^2+k^2}}\cos(x\lambda)\,d\lambda.
\qquad\text{(4.10)}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_{ua}`` | unchanged | mixed Pollaczek correction | dimensionless |
| ``h_1,h_2,x`` | unchanged | overhead height, burial depth, horizontal spacing | m |
| ``c_n,s_n`` | unchanged | fitted residues and complex poles | fit-dependent |
| ``H_n`` | unchanged | effective complex overhead height | m |

No notation was renamed.

## Evidence and approximation sources

Rallis contributes the DCIM fit and finite sum; the underlying kernel remains Pollaczek's. This record is an approximation formula, not a priority claim for the mixed physical model.

## Limitations and discrepancies

- The source's empirical GPOF settings are not validity guarantees.
