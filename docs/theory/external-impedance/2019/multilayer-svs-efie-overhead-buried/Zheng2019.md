# Zheng–Shafieipour multilayer SVS–EFIE line impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Arbitrary conductor boundaries discretized by MoM. Finite conductor cross-sections; cable insulation enters material Green functions. |
| Calculated quantities | Numerical p.u.l. resistance and inductance of overhead or buried multiconductor lines in layered media |
| Earth structure | General multilayer medium; examples use air–soil. |
| Model and approximation | MoM boundary discretization and finite-difference/eigen expansion of the layered Green function; MQS excludes wave-radiation effects. |
| Main source | S. Zheng, M. Shafieipour, J. DeSilva, J. Nordstrom, and V. Okhmatovski (2019) |
| Citation key(s) | `:Zheng2019` |
| Evidence status | Original title/equation pages verified |

**Description.** Surface-volume-surface electric-field integral equation using a multilayer Green function to extract p.u.l. ``R`` and ``L`` for conductors above or below an air–soil interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Magneto-quasistatic longitudinally invariant extraction. | Abstract/§II. |
| Air propagation constant ``γ_air`` | Quasistatic layered Green function. | (6)–(7). |
| Earth propagation constant ``γ_earth`` | Encoded in the layered Green-function ODE eigenvalues ``S_j``. | (7). |
| Earth permittivity and displacement current | Magneto-quasistatic formulation; conductive soil effect retained. | Abstract. |
| Range of validity | Parallel multiconductor lines, overhead or buried. | Abstract. |
| Earth permeability ``μ_earth`` | Source examples use stated layer values. | Model. |
| Arrangement | Arbitrary conductor boundaries discretized by MoM. | §II. |
| Earth structure | General multilayer medium; examples use air–soil. | Introduction/(7). |
| Conductor and insulation geometry | Finite conductor cross-sections; cable insulation enters material Green functions. | §II. |
| Constitutive and field assumptions | Linear media, MQS, surface/volume equivalence. | Formulation. |
| Conventions | ``V_{p.u.l.}`` drives the boundary-current equation. | (1). |

**Expression.**

```math
-j\omega\mu_0\oint_{\partial S}G_\sigma(\rho,\rho')J_z(\rho')d\rho'
+\omega^2\mu_0\sigma\oint_{\partial S}\left[\iint_SG_\varepsilon(\rho,\rho')G_\sigma(\rho',\rho'')ds'\right]
J_z(\rho'')d\rho''=V_{p.u.l.}.
\qquad\text{(1)}
```

The spatial layered Green function is evaluated from

```math
G_\varepsilon(|x-x'|,y_n,y'_m)=
\sum_{j=0}^{N-1}T_{nj}R_{jm}d_m\frac{e^{-\sqrt{S_j}|x-x'|}}{2\sqrt{S_j}}.
\qquad\text{(7)}
```

After MoM solution the terminal admittance is inverted and decomposed as ``\mathbf Z=\mathbf R+j\omega\mathbf L``.

**Approximation.** MoM boundary discretization and finite-difference/eigen expansion of the layered Green function; MQS excludes wave-radiation effects.

**Limitations.** Numerical formulation, not a closed scalar kernel; truncation and mesh convergence are problem dependent.

**Reference.** [Zheng2019](@cite), equations (1), (6)–(7).

**Transcription source.** Original IPST title and equation pages visually checked.

## Source transcription

The source computes both conductor redistribution and interface effects; this record classifies its terminal line impedance by output, not by application label.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_z`` | unchanged | longitudinal auxiliary surface current | ``A/m`` |
| ``G_σ,G_ε`` | unchanged | layered conductive/electric Green functions | source normalized |
| ``S_j`` | unchanged | discrete Green-function eigenvalue | ``m^{-2}`` |

## Limitations and discrepancies

The source's finite-difference database construction is required for evaluation and is not reconstructed here.
