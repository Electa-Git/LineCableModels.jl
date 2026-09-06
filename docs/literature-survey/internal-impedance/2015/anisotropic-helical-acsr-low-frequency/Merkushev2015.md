# Merkushev–Elagin anisotropic helical ACSR internal impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | One central steel wire and six aluminum strands, homogenized as an anisotropic layer; strand radius ``R`` and helix pitch ``h``. |
| Calculated quantities | Low-frequency p.u.l. internal impedance of a single-layer aluminum conductor helically stranded over a magnetic steel core |
| Earth structure | None. |
| Model and approximation | Six discrete strands are homogenized into a continuous anisotropic surface layer; aluminum skin effect is ignored. The core solution retains cylindrical skin effect through Bessel functions. |
| Main source | A. G. Merkushev and I. A. Elagin (2015) |
| Citation key(s) | `:Merkushev2015` |
| Evidence status | Original publication page images checked |

**Description.** First-order commercial-frequency model replacing six helical aluminum strands by an anisotropic conducting layer around a solid magnetic steel core, so axial core magnetization contributes to the conductor internal impedance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally invariant cross-sectional field problem. | Stated — model construction, p. 402. |
| Air propagation constant ``γ_air`` | Not applicable; external return field is outside this internal model. | Scope — pp. 401–402. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Scope. |
| Earth permittivity and displacement current | Neglected under stationary conductor-field approximation. | Stated — p. 402. |
| Range of validity | First-order low-frequency condition ``\omega\ll1/(\mu_0\sigma_{Al}R^2)``; developed for commercial frequency and single-layer ACSR. | Stated — p. 402. |
| Earth permeability ``μ_earth`` | Not applicable; linear core relative permeability ``\mu_{St}`` is retained. | Stated — before (2). |
| Arrangement | One single-layer ACSR conductor. | Stated — abstract and Fig. 1. |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | One central steel wire and six aluminum strands, homogenized as an anisotropic layer; strand radius ``R`` and helix pitch ``h``. | Stated — pp. 401–402. |
| Constitutive and field assumptions | Linear steel magnetics; azimuthally averaged field; strand-boundary perturbations and aluminum strand skin effect neglected. | Stated — p. 402. |
| Conventions | Source uses ``k_{St}^2=-i\omega\mu_{St}\mu_0\sigma_{St}``; ``J_0,J_1`` are ordinary Bessel functions. | Stated — (2). |

**Expression.** The anisotropic layer parameters are

```math
I_w=\Sigma_zE_z|_{r=R}+\Sigma_\varphi E_\varphi|_{r=R},\quad
\Sigma_z=6\pi R^2\sigma_{Al}Q,\quad
\Sigma_\varphi=\theta\Sigma_z,\quad \theta=\frac{2\pi R}{h},
\tag{1}
```

```math
Q=\langle\cos^2\alpha\rangle
=\frac{2}{\pi}\int_1^3\frac{\arccos[(\rho^2+3)/(4\rho)]}{1+\theta^2\rho^2}\,\rho\,d\rho.
```

The p.u.l. internal impedance is

```math
Z_{int}=\frac{k_{St}}{\sigma_{St}h}\,
\frac{\gamma J_0(k_{St}R)}{1+\gamma\theta J_1(k_{St}R)},
\quad k_{St}^2=-i\omega\mu_{St}\mu_0\sigma_{St},\quad \gamma=\frac BA,
\tag{2}
```

```math
A=\frac{k_{St}RJ_0(k_{St}R)}{J_1(k_{St}R)},\qquad
B=\frac{2\Sigma_c}{\theta\Sigma_z}-\theta\frac{k_{St}RJ_1(k_{St}R)}{J_0(k_{St}R)},\qquad
\Sigma_c=\pi R^2\sigma_{St}.
```

**Approximation.** Six discrete strands are homogenized into a continuous anisotropic surface layer; aluminum skin effect is ignored. The core solution retains cylindrical skin effect through Bessel functions.

**Limitations.** Single reinforcement wire, six single-layer strands, linear core permeability and commercial/low frequency. Contact resistance, nonlinear saturation, individual strand fields, multilayer ACSR and external line inductance are excluded.

**Reference.** [Merkushev2015](@cite).  Merkushev and Elagin, *Technical Physics Letters* 41(4), 2015, DOI `10.1134/S1063785015040288`, equations (1)–(2), printed p. 402.

**Transcription source.** Original publication page image. The integration range, ``\rho`` factor, signs, Bessel-function orders, pitch factor and ``\Sigma_c`` definition were visually verified.

## Source transcription

The source states ``h=4\pi R/\tan\alpha_0`` for its director angle. It also warns that high core permeability can produce strong nonlinearity/saturation, outside (2).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\Sigma_z,\Sigma_\varphi`` | unchanged | axial/azimuthal layer conductivities per unit length | source convention |
| ``\sigma_{Al},\sigma_{St}`` | unchanged | aluminum/steel bulk conductivity | ``\mathrm{S/m}`` |
| ``Q`` | unchanged | strand-packing/director form factor | dimensionless |
| ``Z_{int}`` | unchanged | conductor p.u.l. internal impedance | ``\Omega/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

The paper describes its anisotropic implementation as original and (1) as the first-order low-frequency helical-layer approximation. Equation (2) is the derived impedance formula.

## Limitations and discrepancies

- The translated PDF calls ``R`` both the radius of conductors and the core/conductor interface radius; that idealized equality is preserved rather than geometrically repaired.
- The source's conclusions emphasize that a linear ``\mu_{St}`` model is not exhaustive at operating currents.

