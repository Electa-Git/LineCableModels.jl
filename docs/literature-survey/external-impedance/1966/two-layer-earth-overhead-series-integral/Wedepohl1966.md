# Wedepohl–Wasley two-layer overhead series-impedance integral

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Mutual pair at heights ``h_r,y`` and horizontal coordinates ``s_r,x``. Filamentary external conductors; no insulation term. |
| Calculated quantities | Mutual p.u.l. series-impedance contribution of a multilayer earth; conduction-only and displacement-current variants |
| Earth structure | Finite upper layer of thickness ``d`` above a lower half-space. |
| Model and approximation | Equation (8) neglects displacement current. Equation (9) restores it after prescribing free-space longitudinal propagation; it is not a solved full-wave modal result. |
| Main source | L. M. Wedepohl and R. G. Wasley (1966) |
| Citation key(s) | `:Wedepohl1966` |
| Evidence status | Original publication page images checked |

**Description.** Original two-layer-earth integral for multiconductor overhead-line series impedance, including the paper's later displacement-current correction.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in (8); (9) uses the free-space propagation approximation. | Stated in derivation. |
| Air propagation constant ``γ_air`` | Free-space value enters the correction leading to (9). | Stated before (9). |
| Earth propagation constant ``γ_earth`` | Layer transverse roots are built from ``m_2,m_3`` and, in (9), ``k_2,k_3``. | (8)–(9). |
| Earth permittivity and displacement current | Omitted in (8), retained approximately in (9). | Source distinction. |
| Range of validity | Infinite parallel overhead conductors over planar layers. | Problem statement. |
| Earth permeability ``μ_earth`` | Independent ``μ_2,μ_3`` retained. | (8). |
| Arrangement | Mutual pair at heights ``h_r,y`` and horizontal coordinates ``s_r,x``. | (8). |
| Earth structure | Finite upper layer of thickness ``d`` above a lower half-space. | Figure and (8). |
| Conductor and insulation geometry | Filamentary external conductors; no insulation term. | Derivation. |
| Constitutive and field assumptions | Linear isotropic layers and harmonic fields. | Definitions. |
| Conventions | Source phasor convention and principal decaying transverse roots are retained. | Printed kernels. |

**Expression.**

```math
\begin{aligned}
Z_{2rs}&=\frac{m_1^2}{\pi}\int_0^\infty
\frac{(\mu_2/\mu_1)\cos\{\alpha(x-s_r)\}e^{-\alpha(h_r+y)}}
{(\mu_2/\mu_1)\alpha+A}\,d\alpha \\
m_1&=(j\omega\mu_1)^{1/2},
\end{aligned}
```

```math
\begin{aligned}
A&=\frac{q_2}{\sinh(dq_2)}
\left[\cosh(dq_2)-
\frac{\mu_3q_2}{\mu_3q_2\cosh(dq_2)+\mu_2q_3\sinh(dq_2)}\right] \\
q_i&=(\alpha^2+m_i^2)^{1/2},
\end{aligned}
```

where ``m_2=(j\omega\mu_2/\rho_2)^{1/2}`` and ``m_3=(j\omega\mu_3/\rho_3)^{1/2}``. Equation (9) has the same outer integral with ``A'`` obtained by replacing ``q_i`` by ``(\alpha^2+m_i^2+k_i^2)^{1/2}``, where ``k_2^2=\omega^2(\mu_1\varepsilon_1-\mu_2\varepsilon_2)`` and likewise for layer 3.

**Approximation.** Equation (8) neglects displacement current. Equation (9) restores it after prescribing free-space longitudinal propagation; it is not a solved full-wave modal result.

**Limitations.** Two earth layers only; finite-radius, internal, and insulation terms are external to this formula.

**Reference.** [Wedepohl1966](@cite), equations (8)–(9).

**Transcription source.** Directly checked against enlarged PDF page images.

## Source transcription

The two printed variants are retained separately through ``A`` and ``A'`` rather than blended into an inferred hybrid.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``m_i,q_i`` | unchanged | layer diffusion and transverse roots | ``m^{-1}`` |
| ``d`` | unchanged | upper-earth thickness | m |
| ``Z_{2rs}`` | unchanged | earth series contribution | ``Ω/m`` |

## Evidence and approximation sources

The source gives the two-layer integral and series terms directly.

## Limitations and discrepancies

The source does not give a general arbitrary-layer recursion or an explicit modern square-root branch prescription.
