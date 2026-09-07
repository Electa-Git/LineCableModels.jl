# Di Lorenzo et al. three-medium seabed potential matrix

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Each cable has core radius ``r_1``, sheath radii ``r_2,r_3``, and outer insulation radius ``r_4``; inner and outer insulation coefficients are included in ``P'_L``. |
| Calculated quantities | Seabed-return self/mutual potential coefficients and assembled phase-domain per-unit-length admittance for submarine cables |
| Earth structure | Air over finite-depth seawater over semi-infinite seabed. |
| Model and approximation | The parent electromagnetic problem is reduced by the quasi-TEM assumption. Equations (33)–(35) remain infinite spectral integrals and require numerical quadrature; no fit or truncation order is specified. Equation (39) is the ideal coaxial lossless-insulation coefficient used by the source. |
| Main source | G. Di Lorenzo, E. Stracqualursi, M. Marzinotto, J. Brandao Faria, and R. Araneo (2023) |
| Citation key(s) | `:DiLorenzo2023` |
| Evidence status | Original publication page images checked; DOI retained |

**Description.** Three-medium sea/seabed return potential coefficient for buried submarine cables, combined with coaxial insulation coefficients in loop coordinates, transformed to phase coordinates, and inverted as a matrix to obtain admittance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Quasi-TEM: omitted from ``\alpha_n``. | Stated — opening of §III, p. 578. |
| Air propagation constant ``γ_air`` | ``\gamma_0=\sqrt{j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)}``, with ``\sigma_0=0`` in Table I. | Stated — below (26), p. 578. |
| Earth propagation constant ``γ_earth`` | Seabed ``\gamma_2=\sqrt{j\omega\mu_2(\sigma_2+j\omega\varepsilon_2)}``; seawater uses ``\gamma_1`` analogously. | Stated — below (26), p. 578. |
| Earth permittivity and displacement current | Retained through ``\kappa_2=\sigma_2+j\omega\varepsilon_2`` and all three ``\gamma_n``. | Stated — (33) and definition below (26). |
| Range of validity | Same shallow-water, seabed-buried configuration as the impedance formulation; numerical cases are demonstrations, not universal validity bounds. | Stated — §§III–IV. |
| Earth permeability ``μ_earth`` | General layer permeabilities retained; numerical cases set them equal to ``\mu_0``. | Stated — (34)–(35), Table I. |
| Arrangement | Underground/submarine; two coaxial single-core cables are explicitly assembled, with a cited extension to more complex layouts. | Stated — text before (36), p. 579. |
| Earth structure | Air over finite-depth seawater over semi-infinite seabed. | Stated — §II-C and Fig. 2(d). |
| Conductor and insulation geometry | Each cable has core radius ``r_1``, sheath radii ``r_2,r_3``, and outer insulation radius ``r_4``; inner and outer insulation coefficients are included in ``P'_L``. | Stated — Fig. 3 and (37)–(39), p. 579. |
| Constitutive and field assumptions | Linear homogeneous isotropic layers, infinite parallel coaxial cables, quasi-TEM; insulation is lossless in the printed ``\varepsilon_{cs},\varepsilon_{se}`` coefficients. | Stated — §III and (39). |
| Conventions | ``j=\sqrt{-1}``; form the full potential matrix, transform it, then invert it; ``t`` denotes transpose. | Stated — (36)–(40). |

**Expression.** Proposed potential coefficient and admittance assembly, equations (33), (36), and (40), printed p. 579.

```math
P'_{2,ij}=\frac{j\omega}{2\pi\kappa_2}\int_0^\infty
[F_3(\lambda)+G_b(\lambda)]\cos(q_{ij}\lambda)\,d\lambda,
\qquad\text{(33)}
```

```math
G_b(\lambda)=2\mu_1\mu_2\alpha_2e^{-\alpha_2(h_i+h_j-2h_s)}
\left\{
(\gamma_2^2-\gamma_1^2)(s_{10}A_{10}-d_{10}\Delta_{10}e^{-4\alpha_1h_s})
-2\mu_0\mu_1e^{-2\alpha_1h_s}
\right.
```

```math
\left.
\qquad\cdot
[\alpha_0^2\gamma_1^2(\gamma_2^2-\gamma_1^2)
+\alpha_1^2\gamma_0^2(\gamma_2^2+\gamma_1^2)
-2\alpha_1^2\gamma_1^2\gamma_2^2]
\right\},
\qquad\text{(34)}
```

```math
\begin{aligned}
\Delta_{10}&=\alpha_0\gamma_1^2\mu_0-\alpha_1\gamma_0^2\mu_1 \\
A_{10}&=\alpha_0\gamma_1^2\mu_0+\alpha_1\gamma_0^2\mu_1.
\end{aligned}\qquad\text{(35)}
```

Here ``F_3`` is equation (30), transcribed in the companion impedance record, and ``\kappa_2=\sigma_2+j\omega\varepsilon_2`` follows the paper's medium convention.

```math
\begin{aligned}
\mathbf P'&=(\mathbf T^t)^{-1}\mathbf P'_L(\mathbf T)^{-1} \\
\mathbf Y'&=j\omega(\mathbf P')^{-1}.
\end{aligned}\qquad\text{(36,40)}
```

```math
\mathbf P'_L=
\begin{bmatrix}
P'_{cs}&0&0&0\\
0&P'_{se}+P'_{m,ii}&0&P'_{m,ij}\\
0&0&P'_{cs}&0\\
0&P'_{m,ji}&0&P'_{se}+P'_{m,jj}
\end{bmatrix},
\quad
\mathbf T=
\begin{bmatrix}
1&0&0&0\\-1&1&0&0\\0&0&1&0\\0&0&-1&1
\end{bmatrix}.
\qquad\text{(37,38)}
```

```math
\begin{aligned}
P'_{cs}&=\frac{1}{2\pi\varepsilon_{cs}}\ln\left(\frac{r_2}{r_1}\right) \\
P'_{se}&=\frac{1}{2\pi\varepsilon_{se}}\ln\left(\frac{r_4}{r_3}\right).
\end{aligned}\qquad\text{(39)}
```

**Approximation.** The parent electromagnetic problem is reduced by the quasi-TEM assumption. Equations (33)–(35) remain infinite spectral integrals and require numerical quadrature; no fit or truncation order is specified. Equation (39) is the ideal coaxial lossless-insulation coefficient used by the source.

**Numerical interpretation.** The numerator printed in (34) is divided by the magnetic and electric boundary determinants obtained from (27). With ``E=e^{-2\alpha_1h_s}``, these are

```math
\begin{aligned}
D_m&=s_{10}s_{21}-d_{10}d_{21}E,\\
D_e&=A_{10}A_{21}+\Delta_{10}\Delta_{21}E,\\
A_{21}&=\alpha_1\gamma_2^2\mu_1+\alpha_2\gamma_1^2\mu_2,\\
\Delta_{21}&=\alpha_1\gamma_2^2\mu_1-\alpha_2\gamma_1^2\mu_2.
\end{aligned}
```

Thus numerical ``G_b`` equals the right-hand side printed in (34), divided by ``D_mD_e``. This gives it the same spectral units as ``F_3``. Solving the eight boundary conditions (27) independently gives the same electric correction; removing the water layer or making the water and seabed identical recovers the corresponding two-medium coefficient. The displayed source equation remains unchanged.

The direct part of ``F_3`` is integrated as ``K_0(\gamma_2d_{ij})``; the reflected magnetic term and the normalized electric correction are integrated numerically. Material and propagation factors are scaled within the quotients to avoid arithmetic underflow. Both conductors must lie below the water layer. The existing coaxial potential assembly supplies (36)–(39), followed by the full matrix inverse (40).

**Limitations.** The explicitly printed assembly is for two single-core coaxial cables. Matrix inversion is essential: admittance is not the entrywise reciprocal of ``P'_{2,ij}``. Semiconductive screens and lossy insulation are absent. Layer interfaces are planar and media homogeneous; root branches and global numerical-error bounds are unstated.

**Reference.** [DiLorenzo2023](@cite).  Di Lorenzo et al., *IEEE Transactions on Electromagnetic Compatibility* 65(2), 2023, DOI `10.1109/TEMC.2023.3241363`, equations (20), (30)–(40), printed p. 579.

**Transcription source.** Original IEEE PDF page image. The ``j\omega/(2\pi\kappa_2)`` prefactor, every exponent and sign in (34), the transpose/inverse ordering, loop matrix, incidence matrix, logarithmic insulation coefficients, and final matrix inversion were visually verified.

## Source transcription

Equation (32) identifies the two terms used in (33):

```math
E_{2y}=\frac{\partial}{\partial y}
\left[\frac{\partial\Pi'_{2x}}{\partial x}
+\frac{\partial\Pi'_{2z}}{\partial z}\right]
=\frac{\partial^2Q}{\partial x\partial y}.
\qquad\text{(32)}
```

The paper states that the first and second terms inside the square brackets on the left correspond respectively to ``F_3`` and ``G_b``. That source basis is retained rather than treating (33) as an algebraic conversion of (29).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``P'_{2,ij}`` | unchanged | seabed-return potential coefficient | source p.u.l. convention |
| ``P'_{m,ij}`` | unchanged | return-medium entry placed in loop matrix | equals the applicable self/mutual medium coefficient |
| ``P'_{cs},P'_{se}`` | unchanged | core–sheath and sheath–external insulation coefficients | p.u.l. potential coefficients |
| ``\mathbf P'_L,\mathbf P',\mathbf Y'`` | bold added | loop potential, phase potential and phase admittance matrices | matrix operations |
| ``F_3,G_b`` | unchanged | impedance-derived and electric-field correction kernels | spectral functions |
| ``\kappa_2`` | unchanged | seabed complex conductivity | ``\mathrm S/\mathrm m`` |

Boldface only makes source matrices explicit; the algebra and source symbols are unchanged.

## Evidence and approximation sources

The authors explicitly say that (33) is newly proposed for cables buried in the seabed and derive it through the same Hertzian-potential approach as earlier two-medium work. Equations (36)–(40) are the source-prescribed assembly, not a corpus-added scalar conversion.

## Limitations and discrepancies

The printed (34) is preserved. Numerical use includes the boundary determinants specified above; their normalization is checked against the source's Hertz-potential boundary system rather than inferred from the impedance kernel.
