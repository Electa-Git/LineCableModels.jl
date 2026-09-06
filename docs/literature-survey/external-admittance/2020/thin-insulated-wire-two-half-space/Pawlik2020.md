# Pawlik–Woodhouse thin insulated wire external admittance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Scalar single-wire extraction. Core radius ``a``, outer insulation radius ``b``; coating admittance separate. |
| Calculated quantities | QTEM external shunt admittance of an insulated wire above or below a lossy interface |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | qTEM analytical extraction from the paper's modal equation; the inverse is scalar for its one-wire configuration. |
| Main source | Brent Pawlik, Darren J. Woodhouse, and Terrence J. Summers (2020) |
| Citation key(s) | `:Pawlik2020` |
| Evidence status | Original publication page images checked |

**Description.** External shunt counterpart of the paper's insulated-wire qTEM line parameters.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Full modal parent reduced under stated qTEM conditions. | §III. |
| Air propagation constant ``γ_air`` | Retained via medium-1 constants. | Definitions. |
| Earth propagation constant ``γ_earth`` | Retained via medium-2 constants. | Definitions. |
| Earth permittivity and displacement current | Retained through ``σ_1+jωε_1``. | (15). |
| Range of validity | Single infinite thin insulated wire either side of interface. | Abstract. |
| Earth permeability ``μ_earth`` | Unequal half-space values allowed. | Abstract. |
| Arrangement | Scalar single-wire extraction. | §III. |
| Earth structure | Two homogeneous half-spaces. | Figure 1. |
| Conductor and insulation geometry | Core radius ``a``, outer insulation radius ``b``; coating admittance separate. | Figure 1. |
| Constitutive and field assumptions | Linear isotropic media. | §II. |
| Conventions | ``1/Y_s=1/Y_i+1/Y_e``. | (11). |

**Expression.**

```math
\begin{aligned}
\frac1{Y_s}&=\frac1{Y_i}+\frac1{Y_e} \\
Y_e&=2\pi(\sigma_1+j\omega\varepsilon_1)[\Lambda+2(N-jM)]^{-1},
\end{aligned}\qquad\text{(11,15)}
```

where ``\Lambda,N,M`` are defined in (16)–(18).

**Approximation.** qTEM analytical extraction from the paper's modal equation; the inverse is scalar for its one-wire configuration.

**Limitations.** Single wire, two half-spaces, no arbitrary-layer or multiconductor matrix generalization.

**Reference.** [Pawlik2020](@cite), equations (10)–(18).

**Transcription source.** Original paper images; complex constitutive factor and inverse grouping verified.

## Source transcription

The external ``Y_e`` and insulation ``Y_i`` are kept as different taxonomy records.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``M,N`` | unchanged | interface potential correction terms | source normalized |
| ``\Lambda`` | unchanged | geometric log term | dimensionless |
| ``Y_e`` | unchanged | external shunt admittance | ``S/m`` |

## Evidence and approximation sources

The source was inspected across the modal derivation and extracted-parameter section.

## Limitations and discrepancies

No entrywise inverse is inferred for multiconductor use.
