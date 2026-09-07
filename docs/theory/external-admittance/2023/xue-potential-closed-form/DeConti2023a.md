# De Conti–Duarte–Alipio closed-form underground potential coefficients

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Insulated circular cables; self uses ``h_m=h_n`` and outer insulated radius ``r=r_o``. |
| Calculated quantities | Self and mutual entries of ground-return potential-coefficient matrix ``P_g`` and source-prescribed assembled admittance ``Y_g=j\omega P_g^{-1}`` |
| Earth structure | Homogeneous earth below homogeneous air. |
| Model and approximation | The source applies the asymptote ``u_0/(u_0+\gamma_0^2\gamma_1^{-2}u_1)\approx\gamma_1^2/(\gamma_1^2+\gamma_0^2)`` for ``\|\gamma_1\|\gg\|\gamma_0\|`` to the parent integral ``\Theta_2``. The remaining Fourier–Bessel integral is evaluated as ``K_0(\gamma_1D)``, yielding (14) and (15). |
| Main source | A. De Conti, N. Duarte, and R. Alipio (2023) |
| Citation key(s) | `:DeConti2023a` |
| Evidence status | Original publication page images checked |

**Description.** Closed-form approximation of each self or mutual entry of the Maxwell ground-return potential-coefficient matrix for insulated cables buried in homogeneous earth, followed by inversion of the assembled matrix to obtain admittance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed zero in the compact Xue equations. | Stated — discussion below (11), p. 2892. |
| Air propagation constant ``γ_air`` | ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``, air conductivity zero. | Stated — p. 2892. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)}``. | Stated — p. 2892. |
| Earth permittivity and displacement current | Retained in ``\gamma_1`` and the potential prefactor ``j\omega/(\sigma_1+j\omega\varepsilon_1)``. | Stated — (9), p. 2892. |
| Range of validity | Approximation of ``u_0/(u_0+\gamma_0^2\gamma_1^{-2}u_1)`` requires ``|\gamma_1|\gg|\gamma_0|``. Tested over 100 Hz–10 MHz, resistivities 100–10,000 ``\Omega\,\mathrm m`` and typical 0.5–2 m geometries; reported errors are empirical. | Stated — below (14), p. 2893 and §IV. |
| Earth permeability ``μ_earth`` | ``\mu_1=\mu_0``. | Stated — p. 2892. |
| Arrangement | Underground; self and mutual entries assembled for multiple cables. | Stated — Fig. 1 and paragraph preceding (7), p. 2892. |
| Earth structure | Homogeneous earth below homogeneous air. | Stated — Fig. 1. |
| Conductor and insulation geometry | Insulated circular cables; self uses ``h_m=h_n`` and outer insulated radius ``r=r_o``. | Stated — p. 2892. |
| Constitutive and field assumptions | Linear homogeneous media, quasi-TEM Xue parent, nonmagnetic ground. | Stated — §II. |
| Conventions | ``j=\sqrt{-1}``; positive depths; ``P_g`` is assembled first, then inverted as a matrix; per-unit-length relation ``Y_g=j\omega P_g^{-1}``. | Stated — (7), p. 2892. |

**Expression.** Proposed potential-coefficient entry (15) and the source's assembly relation (7), printed pp. 2892–2893.

```math
P^{\mathrm{app}}_{g(m,n)}=
\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[
K_0(\gamma_1d)+
\frac{\gamma_1^2-\gamma_0^2}{\gamma_1^2+\gamma_0^2}K_0(\gamma_1D)
\right],
\qquad\text{(15)}
```

```math
\mathbf Y_g=j\omega\mathbf P_g^{-1}.
\qquad\text{(7)}
```

```math
\begin{aligned}
d&=\sqrt{(h_m-h_n)^2+r^2} \\
D&=\sqrt{(h_m+h_n)^2+r^2},
\end{aligned}
```

with the same ``\gamma_0,\gamma_1`` definitions as in the companion impedance record. Every ``P_{g(m,n)}`` is assembled before inversion; ``(P_g^{-1})_{mn}`` is not ``1/P_{g(m,n)}``.

**Approximation.** The source applies the asymptote ``u_0/(u_0+\gamma_0^2\gamma_1^{-2}u_1)\approx\gamma_1^2/(\gamma_1^2+\gamma_0^2)`` for ``|\gamma_1|\gg|\gamma_0|`` to the parent integral ``\Theta_2``. The remaining Fourier–Bessel integral is evaluated as ``K_0(\gamma_1D)``, yielding (14) and (15).

**Limitations.** Matrix assembly and inversion are indispensable. The approximation's stated bulk-constant ordering excludes arbitrary equality/contrast. It covers one earth half-space, nonmagnetic media, infinite parallel cables and quasi-TEM. It does not include insulation admittance itself. Root branches are unstated.

**Reference.** [DeConti2023a](@cite), equations (2)–(7), (9), (11), and (14)–(15), printed pp. 2892–2893 (PDF pages 2–3).

**Transcription source.** Original IEEE page images. The frequency factor in ``P``, reflection ratio, ``K_0`` image term, self geometry, and the assembled ``j\omega P^{-1}`` operation were visually checked. No scalar reciprocals were introduced.

## Source transcription

The compact parent entry is

```math
P_{g(m,n)}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}[\Lambda+\Theta_2],
\qquad\text{(9)}
```

```math
\Theta_2=2\int_0^\infty\frac{u_0}{u_1}
\frac{e^{-(h_m+h_n)u_1}}{u_0+\gamma_0^2\gamma_1^{-2}u_1}
\cos(r\lambda)\,d\lambda,
\qquad\text{(11)}
```

with ``\Lambda=K_0(\gamma_1d)-K_0(\gamma_1D)``. The source explains that its voltage reference differs from Magalhães et al.; this normalization is retained rather than merged with that record.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``P_{g(m,n)}`` | unchanged | self/mutual ground-return potential coefficient | source-normalized Maxwell coefficient |
| ``\mathbf P_g`` | bold added in prose only | assembled potential-coefficient matrix | inverted as a matrix |
| ``\mathbf Y_g`` | bold added in prose only | assembled ground-return admittance | ``\mathrm S/\mathrm m`` |
| ``\Lambda,\Theta_2`` | unchanged | direct/image Bessel term and interface correction | source normalized |
| ``d,D,h_m,h_n,r`` | unchanged | geometry | m |
| ``u_i,\gamma_i`` | unchanged | spectral roots and bulk constants | ``\mathrm m^{-1}`` |

The boldface in prose only makes the source's matrix assembly explicit; scalar source equations are unchanged.

## Evidence and approximation sources

The original Xue form (2), its decomposition (3), (5), (6), compact form (9), and integral (11) are printed on p. 2892. The approximation and final potential entry are (14)–(15), p. 2893. Equation (7) supplies the only admissible conversion to admittance.

## Limitations and discrepancies

- The source calls ``P_g`` ground-return potential coefficients and explicitly uses matrix inversion; entrywise reciprocal language would be incorrect.
