# Rallis DCIM approximation of Carson's integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent at heights ``h_1,h_2`` and horizontal separation ``x``. |
| Calculated quantities | Finite rational complex-image sum for the homogeneous-earth overhead Carson correction |
| Earth structure | Homogeneous conductive half-space. |
| Model and approximation | The denominator is fitted by complex exponentials with GPOF and each is integrated using the elementary Laplace–cosine identity. The exact Struve/Bessel form (4.18), separately reproduced by the thesis, is not attributed to this DCIM method and is already represented elsewhere in the corpus. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF evaluator; Carson for the parent kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

**Description.** DCIM/GPOF finite-sum approximation to Carson's overhead mutual correction, retained separately from the exact Struve/Bessel identity that the thesis also reproduces.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected by Carson's parent. | Stated through (4.17). |
| Air propagation constant ``γ_air`` | Not retained. | Stated by parent reduction. |
| Earth propagation constant ``γ_earth`` | ``k^2=j\omega\mu_0\sigma``. | Stated below (4.17). |
| Earth permittivity and displacement current | Neglected. | Stated — definition of ``k``. |
| Range of validity | Fit uses a finite, empirically chosen sampling interval and term count; no global error bound is provided. | Stated — §4.4. |
| Earth permeability ``μ_earth`` | ``\mu_0``. | Stated — (4.16). |
| Arrangement | Overhead; mutual parallel conductors. | Stated — Fig. 4.7. |
| Earth structure | Homogeneous conductive half-space. | Stated — §4.4. |
| Conductor and insulation geometry | Filament parent at heights ``h_1,h_2`` and horizontal separation ``x``. | Stated — Fig. 4.7. |
| Constitutive and field assumptions | Linear isotropic conduction-only soil; quasi-static Carson kernel. | Stated/inherited. |
| Conventions | ``H=h_1+h_2``; per-unit-length impedance. | Stated — (4.18). |

**Expression.** The DCIM fit and final Carson correction are

```math
\frac{1}{\lambda+\sqrt{\lambda^2+k^2}}
\approx\sum_{n=1}^{N}c_ne^{s_n\lambda},
\tag{4.20}
```

```math
J_c\approx2\sum_{n=1}^{N}c_n\frac{H_n}{H_n^2+x^2},
\qquad H_n=h_1+h_2-s_n,qquad n=1,2,\ldots,N,
\tag{4.22}
```

with ``Z=j\omega\mu_0[\ln(D/d)+J_c]/(2\pi)``.

**Approximation.** The denominator is fitted by complex exponentials with GPOF and each is integrated using the elementary Laplace–cosine identity. The exact Struve/Bessel form (4.18), separately reproduced by the thesis, is not attributed to this DCIM method and is already represented elsewhere in the corpus.

**Limitations.** Mutual overhead geometry is shown. Fit coefficients are numerical and setup-dependent. This representation inherits Carson's homogeneous, filamentary, conduction-only and quasi-static restrictions.

**Reference.** [Rallis2013](@cite).  K. V. Rallis, doctoral thesis, 2012, DOI `10.12681/eadd/34633`, equations (4.16)–(4.22), printed pp. 78–80.

**Transcription source.** Original thesis page images. The fit, factor two, rational complex-image term and ``H_n=h_1+h_2-s_n`` were visually verified.

## Source transcription

```math
J_c=\int_0^\infty
\frac{2e^{-(h_1+h_2)\lambda}}
{\lambda+\sqrt{\lambda^2+k^2}}\cos(x\lambda)\,d\lambda.
\tag{4.17}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_c`` | unchanged | Carson correction | dimensionless |
| ``d,D`` | unchanged | direct and image distances | m |
| ``c_n,s_n`` | unchanged | fitted residues and poles | fit-dependent |
| ``H_n`` | unchanged | complex image height | m |

No notation was renamed.

## Evidence and approximation sources

The exact special-function equation (4.18) is a previously known evaluation; (4.20)–(4.22) are the thesis's DCIM representation. They are kept distinct to avoid claiming new physics or duplicating the exact identity.

## Limitations and discrepancies

- The thesis heading says DCIM even though an exact special-function evaluation is also available; the record captures only the fitted DCIM formula.
