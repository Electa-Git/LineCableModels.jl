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

and (15) gives the corresponding medium-2 expression with ``\bar\varepsilon_2``, ``k_2^2e^{u_2(x_2+h_2)}``, and the same denominator. The paper then defines

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

Both air-side and earth-side potentials were inspected; one compact representative is displayed above.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\bar\varepsilon_i`` | unchanged | complex permittivity | ``F/m`` |
| ``u_i`` | unchanged | transverse spectral root | ``m^{-1}`` |
| ``w_{ij},y_{ij}`` | unchanged | potential coefficient and shunt quantity | source scalar normalization |

## Evidence and approximation sources

The broad sweep found this under capacitance/railway vocabulary, outside the earlier earth-admittance title set.

## Limitations and discrepancies

The printed ``K_1`` in the image term is dimensionally unusual but source-verbatim.
