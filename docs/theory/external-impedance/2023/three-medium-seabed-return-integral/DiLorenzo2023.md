# Di Lorenzo et al. three-medium seabed-return impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Source cable ``i`` and observation cable ``j`` are at seabed burial coordinates represented by ``h_i,h_j``; horizontal separation is ``q_{ij}``; exterior insulation is not part of the return kernel. |
| Calculated quantities | Self and mutual per-unit-length ground-return impedance for submarine cables buried in a seabed below a finite seawater layer and air |
| Earth structure | Three horizontal media: air ``#0`` for ``z>h_s``, seawater ``#1`` for ``0<z<h_s``, and semi-infinite seabed ``#2`` for ``z<0``. |
| Model and approximation | The derivation assumes quasi-TEM propagation. Within that model, (29)–(30) is presented as the final spectral integral; it is not a fitted or finite-order approximation. Numerical quadrature is still required. |
| Main source | G. Di Lorenzo, E. Stracqualursi, M. Marzinotto, J. Brandao Faria, and R. Araneo (2023) |
| Citation key(s) | `:DiLorenzo2023` |
| Evidence status | Original publication page images checked; DOI retained |

**Description.** Exact spectral-integral expression derived for the longitudinal self or mutual ground-return impedance of parallel submarine cables lying in the lowest of three horizontal media: air, a finite-depth seawater layer, and a semi-infinite seabed.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Quasi-TEM: longitudinal modal propagation is neglected in the transverse roots. | Stated — opening of §III, p. 578. |
| Air propagation constant ``γ_air`` | ``\gamma_0=\sqrt{j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)}``, with the paper's air data ``\sigma_0=0``. | Stated — below (26) and Table I, p. 578. |
| Earth propagation constant ``γ_earth`` | The lowest-medium seabed constant is ``\gamma_2=\sqrt{j\omega\mu_2(\sigma_2+j\omega\varepsilon_2)}``; seawater uses ``\gamma_1`` analogously. | Stated — below (26), p. 578. |
| Earth permittivity and displacement current | Retained in every ``\gamma_n`` through ``\sigma_n+j\omega\varepsilon_n``. | Stated — below (26), p. 578. |
| Range of validity | Derived for shallow-water configurations with cables buried in the seabed; numerical studies use 1 m water depth and the material values in Table I, but do not state universal frequency or geometric bounds. | Stated — §III and §IV, pp. 578–580. |
| Earth permeability ``μ_earth`` | General ``\mu_0,\mu_1,\mu_2`` are retained in the expression; the numerical study sets all three equal to ``\mu_0``. | Stated — (20), (30), Table I. |
| Arrangement | Underground/submarine; self and mutual terms between infinite parallel cables in medium 2. | Stated — text before (28), p. 579. |
| Earth structure | Three horizontal media: air ``#0`` for ``z>h_s``, seawater ``#1`` for ``0<z<h_s``, and semi-infinite seabed ``#2`` for ``z<0``. | Stated — §II-C and §III, pp. 577–578. |
| Conductor and insulation geometry | Source cable ``i`` and observation cable ``j`` are at seabed burial coordinates represented by ``h_i,h_j``; horizontal separation is ``q_{ij}``; exterior insulation is not part of the return kernel. | Stated — Fig. 2(d), (28)–(30). |
| Constitutive and field assumptions | Linear homogeneous isotropic layers; electric Hertzian potential; infinite longitudinal conductors; quasi-TEM. | Stated — §III and (25)–(27). |
| Conventions | ``j=\sqrt{-1}``; medium indices 0/1/2 mean air/sea/seabed; ``h_s`` is seawater-layer thickness in the paper's coordinates; per-unit-length impedance. | Stated — Figs. 2–3 and §III. |

**Expression.** Proposed seabed-return impedance, equations (29)–(30), printed p. 579.

```math
Z'_{2,ij}=\frac{j\omega\mu_2}{2\pi}\int_0^\infty
F_3(\lambda)\cos(\lambda q_{ij})\,d\lambda,
\qquad\text{(29)}
```

```math
F_3(\lambda)=\frac{1}{\alpha_2}\left[
e^{-\alpha_2|h_i-h_j|}
-\frac{s_{10}d_{21}-d_{10}s_{21}e^{-2\alpha_1h_s}}
{s_{10}s_{21}-d_{10}d_{21}e^{-2\alpha_1h_s}}
e^{-\alpha_2(h_i+h_j-2h_s)}
\right].
\qquad\text{(30)}
```

```math
\begin{aligned}
\alpha_n&=\sqrt{\lambda^2+\gamma_n^2} \\
\gamma_n&=\sqrt{j\omega\mu_n(\sigma_n+j\omega\varepsilon_n)},
\end{aligned}
```

```math
\begin{aligned}
s_{10}&=\mu_0\alpha_1+\mu_1\alpha_0 \\
d_{10}&=\mu_0\alpha_1-\mu_1\alpha_0,
\end{aligned}
```

```math
\begin{aligned}
s_{21}&=\mu_2\alpha_1+\mu_1\alpha_2 \\
d_{21}&=\mu_2\alpha_1-\mu_1\alpha_2.
\end{aligned}\qquad\text{(20)}
```

**Approximation.** The derivation assumes quasi-TEM propagation. Within that model, (29)–(30) is presented as the final spectral integral; it is not a fitted or finite-order approximation. Numerical quadrature is still required.

**Limitations.** It applies to three plane, homogeneous layers and infinite parallel conductors in the lowest medium. It does not include cable internal or insulation impedances. Spectral-root branch selection and numerical quadrature details are not stated alongside the equation. The paper evaluates selected sea/seabed cases rather than establishing global error bounds.

**Reference.** [DiLorenzo2023](@cite).  Di Lorenzo et al., *IEEE Transactions on Electromagnetic Compatibility* 65(2), 2023, DOI `10.1109/TEMC.2023.3241363`, equations (20), (25)–(30), and (35), printed pp. 578–579.

**Transcription source.** Original IEEE PDF page images. The factor ``1/(2\pi)``, leading direct term, two finite-layer exponentials, signs and ordering of ``s_{10},d_{10},s_{21},d_{21}`` were visually verified.

## Source transcription

The source obtains (29) after solving the eight Hertzian-potential amplitudes with four boundary conditions at each interface. Its auxiliary quantities in (35), also used by the companion admittance expression, are

```math
\begin{aligned}
\Delta_{10}&=\alpha_0\gamma_1^2\mu_0-\alpha_1\gamma_0^2\mu_1 \\
A_{10}&=\alpha_0\gamma_1^2\mu_0+\alpha_1\gamma_0^2\mu_1.
\end{aligned}\qquad\text{(35)}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z'_{2,ij}`` | unchanged | self/mutual seabed-return impedance | ``\Omega/\mathrm m`` |
| ``\lambda`` | unchanged | transverse spectral variable | ``\mathrm m^{-1}`` |
| ``\alpha_n`` | unchanged | transverse spectral root in medium ``n`` | ``\mathrm m^{-1}`` |
| ``\gamma_n`` | unchanged | bulk propagation constant of medium ``n`` | ``\mathrm m^{-1}`` |
| ``h_s,h_i,h_j,q_{ij}`` | unchanged | water depth, cable vertical coordinates and horizontal spacing | m |
| ``s_{mn},d_{mn}`` | unchanged | permeability-weighted interface sums/differences | source convention |

No source notation was renamed.

## Evidence and approximation sources

The paper attributes the layered Hertzian-potential method to earlier sources but explicitly presents (29)–(30) as its proposed formulation for cables buried in medium 2. Equations (25)–(28) expose the boundary-value derivation; this is an original formula record, not a conversion-only witness.

## Limitations and discrepancies

- The paper uses positive ``h_i,h_j`` in exponent combinations even though the geometric description places medium 2 at ``z<0``; the printed coordinate convention is preserved.
