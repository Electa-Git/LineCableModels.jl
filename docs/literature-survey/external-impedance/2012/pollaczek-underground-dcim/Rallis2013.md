# Rallis DCIM approximation of the underground Pollaczek integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filament parent; depths ``h_1,h_2``, horizontal separation ``x``. |
| Calculated quantities | Finite discrete-complex-image sum for the homogeneous-earth buried/buried Pollaczek correction |
| Earth structure | Homogeneous conductive half-space below air. |
| Model and approximation | Complex poles ``s_n`` and residues ``c_n`` are fitted with GPOF after mapping ``\gamma=P_at+P_b`` over ``t\in[0,1]``; (4.8) then integrates each exponential analytically. This is a data-fitted finite-sum approximation, not an identity. |
| Main source | K. V. Rallis (2012) for the DCIM/GPOF representation; Pollaczek for the physical kernel |
| Citation key(s) | `:Rallis2013` |
| Evidence status | Original-thesis page image verified |

**Description.** Discrete Complex Image Method approximation that fits the Pollaczek spectral denominator by exponentials and integrates each resulting buried/buried term analytically into a finite ``K_1`` sum.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected by the adopted Pollaczek kernel. | Stated through (4.1)–(4.2). |
| Air propagation constant ``γ_air`` | Not retained. | Stated by the conduction-only parent. |
| Earth propagation constant ``γ_earth`` | ``k^2=j\omega\mu_0\sigma``. | Stated below (4.1), p. 68. |
| Earth permittivity and displacement current | Neglected. | Stated by ``k^2=j\omega\mu_0\sigma``. |
| Range of validity | Thesis uses 100 samples, endpoint ``T_0=10|k|`` and 14 terms; reports 8–10 terms normally adequate and up to 14 for higher accuracy. These are empirical settings. | Stated — pp. 70–72. |
| Earth permeability ``μ_earth`` | ``\mu_0``. | Stated — (4.1). |
| Arrangement | Underground; self/mutual parallel conductors. | Stated — Fig. 4.1. |
| Earth structure | Homogeneous conductive half-space below air. | Stated — Fig. 4.1 and (4.2). |
| Conductor and insulation geometry | Filament parent; depths ``h_1,h_2``, horizontal separation ``x``. | Stated — Fig. 4.1. |
| Constitutive and field assumptions | Linear isotropic soil, classical quasi-static Pollaczek assumptions; DCIM is a numerical approximation. | Stated/inherited — §4.2. |
| Conventions | ``\gamma=\sqrt{\lambda^2+k^2}``; positive burial depths; per-unit-length impedance. | Stated — (4.2), (4.5). |

**Expression.** The fitted denominator and resulting finite sum are

```math
\frac{1}{\lambda+\sqrt{\lambda^2+k^2}}
\approx\sum_{n=1}^{N}c_ne^{s_n\sqrt{\lambda^2+k^2}},
\tag{4.5}
```

```math
J_{uu}\approx\sum_{n=1}^{N}
\frac{2c_nkH_n}{\sqrt{H_n^2+x^2}}
K_1\!\left(k\sqrt{H_n^2+x^2}\right),
\qquad H_n=h_1+h_2-s_n.
\tag{4.9}
```

The source inserts this in ``Z_{uu}=j\omega\mu_0[K_0(kd)-K_0(kD)+J_{uu}]/(2\pi)``.

**Approximation.** Complex poles ``s_n`` and residues ``c_n`` are fitted with GPOF after mapping ``\gamma=P_at+P_b`` over ``t\in[0,1]``; (4.8) then integrates each exponential analytically. This is a data-fitted finite-sum approximation, not an identity.

**Limitations.** Coefficients are generated numerically and not tabulated as universal constants. Accuracy depends on sample path, term count, frequency and separation. The inherited physics is filamentary, homogeneous-earth and conduction-only.

**Reference.** [Rallis2013](@cite).  K. V. Rallis, doctoral thesis, 2012, DOI `10.12681/eadd/34633`, equations (4.1)–(4.9), printed pp. 68–70.

**Transcription source.** Original thesis page images. The fitted variable, exponent sign as printed, ``K_1`` argument, prefactor and complex-distance definition were visually verified.

## Source transcription

The identity used termwise is

```math
\int_0^\infty e^{-\beta\sqrt{\gamma^2+x^2}}\cos(bx)\,dx
=\frac{\beta\gamma}{\sqrt{\beta^2+b^2}}K_1\!\left(\gamma\sqrt{\beta^2+b^2}\right).
\tag{4.8}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_{uu}`` | unchanged | buried/buried Pollaczek interface integral | dimensionless |
| ``k`` | unchanged | conduction-only earth constant | ``\mathrm m^{-1}`` |
| ``c_n,s_n`` | unchanged | GPOF residues and poles | fit-dependent |
| ``H_n`` | unchanged | complex image depth | m |

No notation was renamed.

## Evidence and approximation sources

The thesis explicitly presents DCIM as its method; Pollaczek supplies the physical integral and Gradshteyn–Ryzhik the termwise identity. The record therefore classifies (4.9) as a new representation/evaluator, not new ground-return physics.
