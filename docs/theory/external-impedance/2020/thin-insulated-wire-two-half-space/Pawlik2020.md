# Pawlik–Woodhouse thin insulated wire external impedance

## Identification and source

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
\begin{aligned}
Z_s&=Z_i+Z_e \\
Z_e&=\frac{j\omega\mu_1}{2\pi}\,[\Lambda+2(Q-jP)],
\end{aligned}\qquad\text{(10,14)}
```

where ``\Lambda,Q,P`` are the logarithmic and spectral terms defined in (16)–(18). The separately printed internal/insulation contribution is

```math
Z_i=\frac{j\gamma_c}{2\pi(\sigma_c+j\omega\varepsilon_c)a}
\frac{I_0(j\gamma_ca)}{I_1(j\gamma_ca)}
+\frac{j\omega\mu_i}{2\pi}\ln\frac ba.
\qquad\text{(12)}
```

**Approximation.** This record uses the paper's analytical qTEM extraction, not its full modal equation (7), as the p.u.l. formula.

**Limitations.** Single thin circular covered conductor and one interface; no arbitrary layer stack.

**Reference.** [Pawlik2020](@cite), equations (7), (10)–(18).

**Transcription source.** Original pages and title-page author order inspected.

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
[series record](../../2012/pollaczek-hypergeometric-series/Theodoulidis2012.md#identification-and-source).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``P,Q`` | unchanged | interface spectral correction terms | source normalized |
| ``\Lambda`` | unchanged | Direct-minus-image Bessel term; logarithmic only in its zero-transverse-constant limit | dimensionless |
| ``Z_e`` | unchanged | external p.u.l. series impedance | ``Ω/m`` |

## Evidence and approximation sources

The paper gives the extracted line parameters explicitly.

## Limitations and discrepancies

The scalar result is not generalized here to a multiconductor matrix.
