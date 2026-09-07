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

The external dependencies omitted from the earlier transcription are

```math
\Lambda=K_0(j\gamma_1b)-K_0(j\gamma_1B),
\qquad B=\sqrt{4h^2+b^2}.
\qquad\text{(16)}
```

```math
Q-jP=\int_0^\infty
\frac{e^{-2hu_1}\cos(\lambda b)}
{u_1+(\mu_1/\mu_2)u_2}\,d\lambda,
\qquad\text{(17)}
```

```math
N-jM=\int_0^\infty
\frac{[u_1+(\mu_2/\mu_1)u_2]\,e^{-2hu_1}\cos(\lambda b)}
{[u_1+(\mu_1/\mu_2)u_2][(k_2^2/k_1^2)u_1+(\mu_2/\mu_1)u_2]}
\,d\lambda.
\qquad\text{(18)}
```

The source definitions are

```math
\begin{aligned}
k_m&=\omega\sqrt{\mu_m\left(\varepsilon_m-j\frac{\sigma_m}{\omega}\right)},\\
\gamma_m&=\sqrt{\Gamma^2+k_m^2},\\
u_m&=\sqrt{\lambda^2-\gamma_m^2}.
\end{aligned}
```

The conductor is in source medium 1. For a buried conductor, the material
indices are exchanged so that medium 1 denotes the soil. The branch is
selected with ``\Im(\gamma_1)<0`` in the lossy case and continued to
the corresponding outgoing boundary in lossless air. The small-radius
separation of internal and external terms requires
``|\gamma_i b|\ll1`` and ``|\gamma_1 b|\ll1``; it does not set
``\Gamma=0``.

The source samples the insulation surface horizontally:
``x=b,\ y=h``. Thus the self correction retains
``\cos(\lambda b)`` and ``B=\sqrt{4h^2+b^2}``.
When the transverse constant vanishes, the Bessel difference tends to
``\ln(B/b)``. This is not a prescription to replace every self
correction by its zero-horizontal-separation limit.

For a prescribed source attenuation constant, the engine input is
``k_x=j\Gamma``. The external series term matches the same-medium
route of `EarthImpedance.Formula(:Pawlik2018)`. The potential coefficient
``j\omega/Y_e`` matches the same-medium route of
`EarthAdmittance.Formula(:MartinsBritto2024)`.
This comparison covers the scalar external term; it does not attribute a
multiconductor assembly or a solved modal spectrum to the 2020 paper.
The internal and insulation terms remain separate.

The additional analytical evaluations in section III impose further
material and propagation restrictions. In particular, (25) independently
retains the exponential factor missing from the earlier Theodoulidis
hypergeometric expression. That correction is documented in the
[series record](../../../external-impedance/2012/pollaczek-hypergeometric-series/Theodoulidis2012.md#identification-and-source).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``M,N`` | unchanged | interface potential correction terms | source normalized |
| ``\Lambda`` | unchanged | Direct-minus-image Bessel term; logarithmic only in its zero-transverse-constant limit | dimensionless |
| ``Y_e`` | unchanged | external shunt admittance | ``S/m`` |

## Evidence and approximation sources

The source was inspected across the modal derivation and extracted-parameter section.

## Limitations and discrepancies

No entrywise inverse is inferred for multiconductor use.
