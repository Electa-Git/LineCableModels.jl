# Pawlik–Woodhouse thin insulated wire external impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Single wire; scalar extracted parameters. Conductive core radius ``a``, insulation outer radius ``b``. |
| Calculated quantities | QTEM external series impedance of a thin insulated conductor above or below a lossy interface |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | This record uses the paper's analytical qTEM extraction, not its full modal equation (7), as the p.u.l. formula. |
| Main source | Brent Pawlik, Darren J. Woodhouse, and Terrence J. Summers (2020) |
| Citation key(s) | `:Pawlik2020` |
| Evidence status | Original publication page images checked |

**Description.** Full modal parent and qTEM extraction of line parameters for an insulated wire on either side of an interface between homogeneous half-spaces.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained in modal parent; qTEM restrictions applied for analytical parameters. | (7), §III. |
| Air propagation constant ``γ_air`` | Independent medium wave number. | Nomenclature. |
| Earth propagation constant ``γ_earth`` | Independent second-half-space wave number. | Nomenclature. |
| Earth permittivity and displacement current | Retained through ``σ+jωε``. | (13)–(18). |
| Range of validity | Infinite thin insulated wire above or below one plane interface. | Abstract/Fig. 1. |
| Earth permeability ``μ_earth`` | Differing half-space permeability allowed. | Abstract. |
| Arrangement | Single wire; scalar extracted parameters. | Derivation. |
| Earth structure | Two homogeneous half-spaces. | Figure 1. |
| Conductor and insulation geometry | Conductive core radius ``a``, insulation outer radius ``b``. | Figure 1. |
| Constitutive and field assumptions | Linear isotropic homogeneous regions; thin-wire filament external field. | §II. |
| Conventions | ``e^{j\omega t-\Gamma z}``; ``Z_s=Z_i+Z_e``. | Nomenclature, (10). |

**Expression.**

```math
Z_s=Z_i+Z_e,
\qquad
Z_e=\frac{j\omega\mu_1}{2\pi}\,[\Lambda+2(Q-jP)],
\qquad\text{(10,14)}
```

where ``\Lambda,Q,P`` are the logarithmic and spectral terms defined in (16)–(18). The separately printed internal/insulation contribution is

```math
Z_i=\frac{j\gamma_c}{2\pi(\sigma_c+j\omega\epsilon_c)a}
\frac{I_0(j\gamma_ca)}{I_1(j\gamma_ca)}
+\frac{j\omega\mu_i}{2\pi}\ln\frac ba.
\qquad\text{(12)}
```

**Approximation.** This record uses the paper's analytical qTEM extraction, not its full modal equation (7), as the p.u.l. formula.

**Limitations.** Single thin circular covered conductor and one interface; no arbitrary layer stack.

**Reference.** [Pawlik2020](@cite), equations (7), (10)–(18).

**Transcription source.** Original pages and title-page author order inspected.

## Source transcription

``Z_e`` remains distinct from ``Z_i``; the cable application does not move it into the internal taxonomy.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``P,Q`` | unchanged | interface spectral correction terms | source normalized |
| ``\Lambda`` | unchanged | logarithmic geometric term | dimensionless |
| ``Z_e`` | unchanged | external p.u.l. series impedance | ``Ω/m`` |

## Evidence and approximation sources

The paper gives the extracted line parameters explicitly.

## Limitations and discrepancies

The scalar result is not generalized here to a multiconductor matrix.
