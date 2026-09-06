# Ametani–Miyamoto–Mahseredjian displacement-current complex-depth approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel round conductors; direct distance ``d_{ij}``, image distance modified by ``2h'_e``. |
| Calculated quantities | Closed complex-image approximation for overhead self/mutual earth return with air-referenced earth permittivity |
| Earth structure | Homogeneous earth below air. |
| Model and approximation | It replaces the exact distributed return by a single complex image at depth ``2h'_e``, following Gary–Déri but changing the penetration constant to the air-referenced displacement-current form. |
| Main source | A. Ametani, Y. Miyamoto, and J. Mahseredjian (2014), modifying the Gary–Déri complex-depth form |
| Citation key(s) | `:Ametani2014` |
| Evidence status | Original publication page images checked |

**Description.** Closed overhead complex-image formula replacing the classical conduction-only complex depth by an air-referenced complex depth that retains ``\epsilon_e-\epsilon_0``.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Represented by subtraction of the air displacement term in the effective depth. | Stated — derivation from (12)–(14). |
| Air propagation constant ``γ_air`` | Air term ``j\omega\epsilon_0`` is subtracted through ``\epsilon_e-\epsilon_0``. | Stated — (18), p. 937. |
| Earth propagation constant ``γ_earth`` | Replaced by ``1/h'_e=\sqrt{j\omega\mu_0[\sigma_e+j\omega\epsilon_0(\epsilon_r-1)]}``. | Stated — (18). |
| Earth permittivity and displacement current | Retained relative to air. | Stated — (18). |
| Range of validity | Presented as approximate and tested over selected heights, resistivities, permittivities and high-frequency cases; the paper reports “reasonable accuracy,” not a formal bound. | Stated — §§2.6–4. |
| Earth permeability ``μ_earth`` | ``\mu_e=\mu_0`` in the proposed complex depth. | Stated — derivation around (12). |
| Arrangement | Overhead multiconductor; self and mutual. | Stated — Fig. 1 and §3.2. |
| Earth structure | Homogeneous earth below air. | Stated — §2.5 homogeneous limit. |
| Conductor and insulation geometry | Infinite parallel round conductors; direct distance ``d_{ij}``, image distance modified by ``2h'_e``. | Stated — (17) and Fig. 1. |
| Constitutive and field assumptions | Linear isotropic nonmagnetic media; quasi-TEM/complex-image engineering approximation. | Stated/inherited — §§2.5–2.7. |
| Conventions | ``j=\sqrt{-1}``; ``\epsilon_e=\epsilon_r\epsilon_0``; p.u.l. impedance. | Stated — (8), (17)–(18). |

**Expression.**

```math
Z_{ij}=j\omega\frac{\mu_0}{2\pi}\ln\left(\frac{S'_{ij}}{d_{ij}}\right),
\qquad\text{(17)}
```

```math
S'_{ij}=\sqrt{(h_i+h_j+2h'_e)^2+y^2},
\qquad
\frac1{h'_e}=\sqrt{j\omega\mu_0
[\sigma_e+j\omega\epsilon_0(\epsilon_r-1)]}.
\qquad\text{(18)}
```

**Approximation.** It replaces the exact distributed return by a single complex image at depth ``2h'_e``, following Gary–Déri but changing the penetration constant to the air-referenced displacement-current form.

**Limitations.** Homogeneous nonmagnetic planar earth and infinite parallel overhead conductors. Accuracy is empirical and degrades with parameter combinations shown in the source. This is not the source's exact integral and does not retain full TM/TE modal propagation.

**Reference.** [Ametani2014](@cite), equations (15), (17)–(18), printed p. 937 (PDF page 2).

**Transcription source.** Original IEEJ page image. The primed image distance, factor ``2h'_e``, logarithm and ``\epsilon_r-1`` term were visually verified.

## Source transcription

The parent Gary–Déri expression is printed with ``h_e^{-1}=\sqrt{j\omega\mu_0\sigma_e}``; (17)–(18) change only that depth prescription while keeping the image geometry.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``h'_e`` | unchanged | proposed air-referenced complex depth | m |
| ``S'_{ij}`` | unchanged | complex source-to-image distance | m |
| ``d_{ij}`` | unchanged | direct conductor distance/self radius | m |
| ``y,h_i,h_j`` | unchanged | horizontal spacing and heights | m |

No notation was renamed.

## Evidence and approximation sources

The paper explicitly labels (17) as the proposed approximate formula and compares it against the spectral result (12)–(14). This record is therefore separate from the companion integral record.

## Limitations and discrepancies

- The source's equation (14) has a missing ``j\omega`` conflict; equation (18) itself visibly retains the complete ``j\omega\epsilon_0(\epsilon_r-1)`` term and is copied exactly.
