# Bormann–Tavakoli reluctance-network series impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Arbitrary conductor shapes connected by reluctance graph. Polygonal/arbitrary conductor regions and gaps. |
| Calculated quantities | Approximate series-impedance matrix for arbitrary multiconductor cross-sections |
| Earth structure | None. |
| Model and approximation | Flux is confined to constructed channels; the added/subtracted internal term corrects the network's conductor penetration approximation. |
| Main source | Dierk Bormann and Hanif Tavakoli (2013) |
| Citation key(s) | `:Bormann2013` |
| Evidence status | Original publication page images checked |

**Description.** Magnetic-reluctance network representation of conductor proximity and skin effects, corrected by an independent internal-impedance term.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in two-dimensional parameter extraction. | Model. |
| Air propagation constant ``γ_air`` | Not applicable to reluctance approximation. | Scope. |
| Earth propagation constant ``γ_earth`` | Not included. | Scope. |
| Earth permittivity and displacement current | Not included. | Scope. |
| Range of validity | Parallel multiconductor cross-sections representable by flux channels. | Method. |
| Earth permeability ``μ_earth`` | Not applicable. | Scope. |
| Arrangement | Arbitrary conductor shapes connected by reluctance graph. | Title/method. |
| Earth structure | None. | Internal model. |
| Conductor and insulation geometry | Polygonal/arbitrary conductor regions and gaps. | Figures. |
| Constitutive and field assumptions | Linear harmonic magnetic diffusion in source's channel approximation. | (2)–(8). |
| Conventions | Incidence matrix ``∂`` maps graph channels to conductors. | (1). |

**Expression.**

```math
\mathbf Z=j\omega(\partial^T\boldsymbol{\mathcal R}\partial)^{-1},
\tag{1}
```

and the corrected final model is

```math
\mathbf Z=\mathbf Z_{int}-\mathbf Z_{int}^{\mathcal R}
+j\omega(\partial^T\boldsymbol{\mathcal R}\partial)^{-1},
\tag{7}
```

```math
Z_{int,k}^{\mathcal R}=
\frac{j\omega\ell\mu_k}{\kappa_k\sum_{i:\partial_{ik}\ne0}w_i\theta_{k,i}},
\quad \theta_{\pm,i}=\tanh(\kappa_{\pm,i}d_{\pm,i}/2).
\tag{3,8}
```

**Approximation.** Flux is confined to constructed channels; the added/subtracted internal term corrects the network's conductor penetration approximation.

**Limitations.** Accuracy depends on reluctance-network topology and geometric discretization; no earth return.

**Reference.** [Bormann2013](@cite), equations (1)–(8).

**Transcription source.** Original PDF images; inverse placement and correction signs verified.

## Source transcription

The low- and high-frequency channel reluctances are source equations (2)–(6).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\mathcal R`` | unchanged | diagonal channel-reluctance matrix | source normalized |
| ``∂`` | unchanged | graph incidence matrix | dimensionless |
| ``\kappa`` | unchanged | conductor diffusion root | ``m^{-1}`` |

## Evidence and approximation sources

Found by arbitrary-geometry and proximity vocabulary in full text.

## Limitations and discrepancies

The source approximation ``\widetilde\theta=\tanh(\kappa d/4)`` must not be confused with the channel ``\theta`` above.
