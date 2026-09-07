# Mariscotti generalized potential in air and lossy earth

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Scalar two-conductor setup; self/mutual positions. Filamentary external field; no coating. |
| Calculated quantities | Self/mutual generalized-potential coefficient and source shunt-admittance extraction |
| Earth structure | Air over homogeneous lossy earth. |
| Model and approximation | Generalized potential within the source's thin-wire/modal model; no extra closed-form approximation is applied to the displayed integral. |
| Main source | Andrea Mariscotti (2019), building on D’Amore–Sarto |
| Citation key(s) | `:Mariscotti2019` |
| Evidence status | Original publication page images checked |

**Description.** Generalized-potential formulation for conductors in air or lossy earth, with application to railway conductors.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained as ``γ`` in generalized potentials. | (12)–(15). |
| Air propagation constant ``γ_air`` | ``k_1``/``χ_1``. | Definitions. |
| Earth propagation constant ``γ_earth`` | ``k_2``/``χ_2``. | Definitions. |
| Earth permittivity and displacement current | Retained using ``\bar\varepsilon=\varepsilon+\sigma/(j\omega)``. | Material definitions. |
| Range of validity | Infinite parallel thin conductors in either half-space. | Geometry. |
| Earth permeability ``μ_earth`` | Independent medium permeability retained. | Definitions. |
| Arrangement | Scalar two-conductor setup; self/mutual positions. | §3. |
| Earth structure | Air over homogeneous lossy earth. | Figure 1. |
| Conductor and insulation geometry | Filamentary external field; no coating. | Derivation. |
| Constitutive and field assumptions | Linear homogeneous isotropic media. | Model. |
| Conventions | ``E_z=fI+\gamma^2gI``, ``V=\gamma gI``. | (12)–(13). |

**Expression.** In medium 1,

```math
V_1=\frac{\gamma I e^{-\gamma z}}{j\omega2\pi\bar\varepsilon_1}
\left[K_0(\chi_1R)-K_1(\chi_1R')+
\int_{-\infty}^{\infty}
\frac{k_1^2e^{-u_1(x_1+h_1)}e^{-j\lambda(y_1-d_1)}}{k_1^2u_2+k_2^2u_1}\,d\lambda\right],
\qquad\text{(14)}
```

The corresponding expression in medium 2 is

```math
V_2=\frac{\gamma I e^{-\gamma z}}{j\omega2\pi\bar\varepsilon_2}
\left[K_0(\chi_2R)-K_1(\chi_2R')+
\int_{-\infty}^{\infty}
\frac{k_2^2e^{u_2(x_2+h_2)}e^{-j\lambda(y_2-d_2)}}{k_1^2u_2+k_2^2u_1}\,d\lambda\right],
\qquad\text{(15)}
```

The material and distance definitions are

```math
\begin{aligned}
\bar\varepsilon_m&=\varepsilon_m+\frac{\sigma_m}{j\omega},
& k_m^2&=\omega^2\mu_m\bar\varepsilon_m,\\
\chi_m^2&=-k_m^2-\gamma^2,
& u_m&=\sqrt{\lambda^2+\chi_m^2},\\
R&=\sqrt{(x-h_m)^2+(y-d_m)^2},
& R'&=\sqrt{(x+h_m)^2+(y-d_m)^2}.
\end{aligned}
```

The source gives the following separate propagation approximations:

```math
\begin{aligned}
\gamma&\simeq jk=j\sqrt{j\omega\mu(\sigma+j\omega\varepsilon)},
&&\text{(17a)}\\
\gamma&\simeq jk=j\sqrt{j\omega\mu(j\omega\varepsilon)}.
&&\text{(17b)}
\end{aligned}
```

The paper then defines

```math
\begin{aligned}
w_{ij}&=\frac{V}{I_t}=\gamma g \\
y_{ij}&=w_{ij}^{-1} \\
y_{ij}&=g_{ij}+j\omega c_{ij}.
\end{aligned}\qquad\text{(19,21)}
```

**Approximation.** Generalized potential within the source's thin-wire/modal model; no extra closed-form approximation is applied to the displayed integral.

**Limitations.** The scalar inverse in (19) is specific to the source setup and must not be read as elementwise inversion of a general potential matrix.

**Reference.** [Mariscotti2019](@cite), equations (12)–(21).

**Transcription source.** Original page images; the unusual printed ``K_1`` image term is retained rather than silently repaired.

## Source transcription

Both source potentials and the separate propagation approximations are displayed above.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\bar\varepsilon_i`` | unchanged | complex permittivity | ``F/m`` |
| ``u_i`` | unchanged | transverse spectral root | ``m^{-1}`` |
| ``w_{ij},y_{ij}`` | unchanged | potential coefficient and shunt quantity | source scalar normalization |

## Evidence and approximation sources

The broad sweep found this under capacitance/railway vocabulary, outside the earlier earth-admittance title set.

## Limitations and discrepancies

The printed image term and propagation radicals remain distinguishable from the numerical interpretation below.

## Numerical interpretation

The numerical coefficient is ``P=V/Q``, using charge continuity
``Q=\gamma I/(j\omega)``. With the engine longitudinal wavenumber
``k_x=j\gamma``, define

```math
\begin{aligned}
g_m&=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m),\\
\chi_m^2&=g_m+k_x^2,\qquad q_m=\sqrt{\lambda^2+\chi_m^2}.
\end{aligned}
```

The outgoing root has positive real part, with the lossless limit
continued from passive material. For a source and receiver in medium
``m``, the coefficient is

```math
P_{ij}=\frac{j\omega}{2\pi(\sigma_m+j\omega\varepsilon_m)}
\left[
K_0(\chi_m R)-K_0(\chi_m R')
+2\int_0^\infty
\frac{g_m e^{-Hq_m}\cos(x_{ij}\lambda)}
{g_1q_2+g_2q_1}\,d\lambda
\right].
```

Here ``H`` is the sum of the positive distances to the interface and
``x_{ij}`` is the lateral separation. The image uses the order-zero
two-dimensional potential Green function. Thus the printed ``K_1``
is not used as a distinct physical correction. The default air reference
follows ``\gamma\simeq jk_1`` together with (8), giving
``k_x^2=-g_1`` in lossless air. The additional imaginary factor in
the printed radical (17b) is not applied to that definition. A supplied
``k_x`` instead evaluates the prescribed-propagation potential;
no modal root search is implied.

The source's scalar reciprocal in (19) is not applied entry by entry.
The assembled matrix uses ``\mathbf Y=j\omega\mathbf P^{-1}``.
Self evaluation uses the wire boundary radius for ``R``; mutual
evaluation uses the centre separation. Only same-half-space pairs are
registered. The earth-side numerator in (15) is retained and is not
replaced by the different numerator of Xue's buried potential kernel.
