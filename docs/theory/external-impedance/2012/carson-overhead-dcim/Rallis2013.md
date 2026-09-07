# Rallis DCIM approximation of Carson's integral

## Identification and source

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
\qquad\text{(4.20)}
```

```math
\begin{aligned}
J_c&\approx2\sum_{n=1}^{N}c_n\frac{H_n}{H_n^2+x^2} \\
H_n&=h_1+h_2-s_n \\
n&=1,2,\ldots,N,
\end{aligned}\qquad\text{(4.22)}
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
\qquad\text{(4.17)}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_c`` | unchanged | Carson correction | dimensionless |
| ``d,D`` | unchanged | direct and image distances | m |
| ``c_n,s_n`` | unchanged | fitted residues and poles | fit-dependent |
| ``H_n`` | unchanged | complex image height | m |

No notation was renamed.

## Numerical interpretation

Select `EarthImpedance.Formula(:Pollaczek1926;evaluation=:dcim)` for
one-level GPOF, or `evaluation=:dcim_two_level` for two-level fitting.
These selections change all three placement evaluators, not the physical
Pollaczek assumptions. The default spectral evaluator is unchanged.

The fit uses 100 uniformly spaced samples and at most 14 retained singular
directions per level. Appendix A defines the shifted sample matrices
``Y_1,Y_2``; for ``Y_1=UDV^H``, the image poles follow from the
eigenvalues of ``D^{-1}U^HY_2V``. A least-squares solve determines the
residues. Numerically unresolved directions and growing exponentials are
excluded before the residue solve. `dcim_samples` and `dcim_terms`
expose the sample and rank limits.

The sampling variable is divided by ``|k|`` before fitting. The
two-level implementation first fits the far interval, then fits its
residual on the near interval. Each retained residue is stored with its
sampling origin; this avoids forming overflowing coefficients before
multiplication by a decaying Bessel function. Coefficients are generated
in double precision, including for higher-precision output arithmetic.
Neither higher output precision nor the sample count constitutes an
error bound for the finite fit.

The fitted spectral function is integrated to infinity as prescribed by
the image sum. No parent-integral fallback is hidden behind this selection.
Tests separately verify the matrix-pencil recovery, integration of the
fitted function, parent comparisons, and matrix assembly. In particular,
large lateral separation can expose substantial buried-conductor fit
error; the spectral or exact-series selections remain available.

The one-level interval is ``0\le\lambda/|k|\le28``, following
the construction of §4.3.2 used by §4.4. Two-level intervals are
``[5,100]`` and ``[0,5]``. The conductor radius enters only
the perfect-ground self logarithm; the self correction sets ``x=0``.


## Evidence and approximation sources

The exact special-function equation (4.18) is a previously known evaluation; (4.20)–(4.22) are the thesis's DCIM representation. They are kept distinct to avoid claiming new physics or duplicating the exact identity.

## Limitations and discrepancies

- The thesis heading says DCIM even though an exact special-function evaluation is also available; the record captures only the fitted DCIM formula.
