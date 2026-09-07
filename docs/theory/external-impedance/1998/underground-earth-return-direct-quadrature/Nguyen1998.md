# Nguyen underground earth-return direct quadrature

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Self/mutual pair ``j,k``. Cable axes; conductor/insulation terms are separate. |
| Calculated quantities | Underground-cable earth-return self/mutual impedance and direct trapezoidal evaluator |
| Earth structure | Homogeneous half-space. |
| Model and approximation | The physical integral is inherited; the trapezoidal truncation introduces ``\Delta u`` and finite ``N`` numerical error. |
| Main source | T. T. Nguyen (1998) |
| Citation key(s) | `:Nguyen1998` |
| Evidence status | Original publication page images checked |

**Description.** A dimensionless direct numerical integration form of the Pollaczek underground earth-return impedance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in the underlying low-frequency earth-return reduction. | Parent derivation. |
| Air propagation constant ``γ_air`` | Not independently solved. | (38). |
| Earth propagation constant ``γ_earth`` | Represented through ``\alpha=\omega\mu_0/\rho_e``. | (39)–(41). |
| Earth permittivity and displacement current | Neglected. | Material definition. |
| Range of validity | Parallel buried cables in homogeneous earth. | Title and geometry. |
| Earth permeability ``μ_earth`` | ``\mu_0``. | (38). |
| Arrangement | Self/mutual pair ``j,k``. | (38). |
| Earth structure | Homogeneous half-space. | Derivation. |
| Conductor and insulation geometry | Cable axes; conductor/insulation terms are separate. | Assembly. |
| Constitutive and field assumptions | Linear conductive earth, harmonic steady state. | Definitions. |
| Conventions | Complex square roots and ``K_0`` follow source convention. | (38). |

**Expression.**

```math
z_{jk}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(d_{jk}\sqrt{j\alpha})-K_0(d'_{jk}\sqrt{j\alpha})\right]
+\frac{\sqrt{\omega\mu_0\rho_e}}{\pi h_m}J(j,k),
\qquad\text{(38)}
```

```math
\begin{aligned}
J(j,k)&=\int_0^\infty h_m\sqrt\alpha\,[\sqrt{u^2+j}-u]
e^{-2h_m\sqrt\alpha\sqrt{u^2+j}}
\cos(x_{jk}\sqrt\alpha\,u)\,du \\
h_m&=\frac{h_j+h_k}{2}.
\end{aligned}
```

The direct evaluator is

```math
J_N=\left[\tfrac12f(0)+\sum_{n=1}^{N-1}f(n\Delta u)+\tfrac12f(N\Delta u)\right]\Delta u.
\qquad\text{(47)}
```

**Approximation.** The physical integral is inherited; the trapezoidal truncation introduces ``\Delta u`` and finite ``N`` numerical error.

**Limitations.** Homogeneous earth and parallel underground cables; no displacement current.

**Reference.** [Nguyen1998](@cite), equations (38)–(47).

**Transcription source.** Original-page images and Markdown were cross-checked.

## Source transcription

Equations (43)–(46) give the recursive accumulation ``J_n``; (47) is the compact final quadrature.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``h_m`` | unchanged | mean burial depth | m |
| ``\alpha`` | unchanged | earth diffusion parameter | ``m^{-2}`` |
| ``J`` | unchanged | dimensionless residual integral | source normalized |

## Evidence and approximation sources

The source gives the quadrature formulation directly.

## Limitations and discrepancies

Step-size and stopping criteria remain evaluator inputs rather than universal source bounds.
