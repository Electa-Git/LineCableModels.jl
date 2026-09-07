# D’Amore–Sarto multiconductor lossy-ground shunt matrix

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Arbitrary conductor count in matrix/modal form. Bare external conductor boundaries; finite coating is not a separate formula. |
| Calculated quantities | Multiconductor p.u.l. shunt-admittance matrix over lossy ground |
| Earth structure | Homogeneous half-space. |
| Model and approximation | Modal/numerical evaluation is required; no entrywise scalar inversion is introduced. |
| Main source | Marcello D’Amore and Maria Sabrina Sarto (1996), Part II |
| Citation key(s) | `:DAmore1996` |
| Evidence status | Original publication page images checked |

**Description.** Wide-frequency multiconductor shunt matrix formed from external-air and ground potential operators.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Determined through modal TL quantities. | Part II derivation. |
| Air propagation constant ``γ_air`` | Enters ``Y_e'``. | (38). |
| Earth propagation constant ``γ_earth`` | Enters ``Y_g'``. | (40). |
| Earth permittivity and displacement current | Retained. | Wide-frequency derivation. |
| Range of validity | Parallel multiconductor lines above a plane lossy interface. | Problem statement. |
| Earth permeability ``μ_earth`` | Source homogeneous-ground value. | Definitions. |
| Arrangement | Arbitrary conductor count in matrix/modal form. | (32)–(34). |
| Earth structure | Homogeneous half-space. | Geometry. |
| Conductor and insulation geometry | Bare external conductor boundaries; finite coating is not a separate formula. | Model. |
| Constitutive and field assumptions | Linear homogeneous isotropic media, harmonic fields. | Derivation. |
| Conventions | ``-d\mathbf I/dx=\mathbf Y'\mathbf V``. | (32). |

**Expression.**

```math
\begin{aligned}
-\frac{d\mathbf I}{dx}&=\mathbf Y'\mathbf V \\
\mathbf Y'&=\mathbf Y_g'(\mathbf Y_e'+\mathbf Y_g')^{-1}\mathbf Y_e',
\end{aligned}\qquad\text{(32–33)}
```

```math
\begin{aligned}
\mathbf Y_e'&=(\boldsymbol\zeta_e')^{-1} \\
\mathbf Y_g'&=(\boldsymbol\zeta_g')^{-1},
\end{aligned}\qquad\text{(34)}
```

with modal components

```math
\begin{aligned}
\widehat{\mathbf Y}_e'&=j\omega\varepsilon_0,2\pi\widehat{\boldsymbol\Lambda}^{-1} \\
\widehat{\mathbf Y}_g'&=j\omega\varepsilon_0\pi
(\widehat{\mathbf S}_{2g}^{h}-\widehat{\mathbf S}_{2g}^{0})^{-1}.
\end{aligned}\qquad\text{(38,40)}
```

**Approximation.** Modal/numerical evaluation is required; no entrywise scalar inversion is introduced.

**Limitations.** Homogeneous planar ground and the source's modal basis.

**Reference.** [DAmore1996](@cite), equations (32)–(40).

**Transcription source.** Original page images; product and inverse ordering verified.

## Source transcription

The matrix series combination in (33) is preserved literally.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``ζ_e',ζ_g'`` | unchanged | external and ground potential matrices | source normalized |
| ``Y_e',Y_g'`` | unchanged | partial shunt matrices | ``S/m`` |
| ``Y'`` | unchanged | total p.u.l. shunt matrix | ``S/m`` |

## Evidence and approximation sources

The source gives the potential and admittance equations directly.

## Limitations and discrepancies

The modal definitions from Part I/preceding sections remain dependencies.
