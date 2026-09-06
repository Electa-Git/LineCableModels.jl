# D’Amore–Sarto multiconductor lossy-ground series matrix

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Arbitrary conductor count represented by modal/physical matrices. External and ground formula is separated from ``Z_i'``; no finite dielectric shell here. |
| Calculated quantities | Multiconductor p.u.l. series matrix split into internal, external, and ground contributions |
| Earth structure | Homogeneous lossy half-space. |
| Model and approximation | The displayed matrices belong to the paper's wide-frequency modal formulation; numerical modal truncation/evaluation is required. |
| Main source | Marcello D’Amore and Maria Sabrina Sarto (1996), Part II |
| Citation key(s) | `:DAmore1996` |
| Evidence status | Original publication page images checked |

**Description.** Wide-frequency multiconductor transmission-line matrix over lossy ground, with modal quantities transformed back to physical conductor coordinates.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Solved through the modal TL system rather than prescribed in the final matrix. | Part II derivation. |
| Air propagation constant ``γ_air`` | Enters modal external field quantities. | (38). |
| Earth propagation constant ``γ_earth`` | Enters the ground modal operators. | (40). |
| Earth permittivity and displacement current | Retained through complex ground fields. | Wide-frequency model. |
| Range of validity | Parallel multiconductor lines above a flat interface. | Problem statement. |
| Earth permeability ``μ_earth`` | Source uses the stated homogeneous ground constants. | Definitions. |
| Arrangement | Arbitrary conductor count represented by modal/physical matrices. | (26). |
| Earth structure | Homogeneous lossy half-space. | Part II geometry. |
| Conductor and insulation geometry | External and ground formula is separated from ``Z_i'``; no finite dielectric shell here. | (25). |
| Constitutive and field assumptions | Linear homogeneous isotropic media, harmonic fields. | Derivation. |
| Conventions | ``-d\mathbf V/dx=\mathbf Z'\mathbf I``. | (24). |

**Expression.**

```math
\begin{aligned}
-\frac{d\mathbf V}{dx}&=\mathbf Z'\mathbf I \\
\mathbf Z'&=\mathbf Z_i'+\mathbf Z_e'+\mathbf Z_g'.
\end{aligned}\qquad\text{(24–25)}
```

In the modal basis the paper gives

```math
\begin{aligned}
\widehat{\mathbf Z}_e'&=\frac{j\omega\mu_0}{2\pi}\widehat{\boldsymbol\Lambda} \\
\widehat{\mathbf Z}_g'&=\frac{j\omega\mu_0}{\pi}
\left(\widehat{\mathbf S}_{1g}^{h}-k_0^{-2}\widehat{\mathbf S}_{2g}^{0}\widehat{\boldsymbol\Lambda}\right),
\end{aligned}\qquad\text{(38,40)}
```

with physical entries obtained by the transformation printed in (26).

**Approximation.** The displayed matrices belong to the paper's wide-frequency modal formulation; numerical modal truncation/evaluation is required.

**Limitations.** Homogeneous planar ground and parallel conductors; the record does not reinterpret modal entries as independent scalar wires.

**Reference.** [DAmore1996](@cite), equations (24)–(26), (38), (40).

**Transcription source.** Original Part II page images, including operator ordering and the ``k_0^{-2}`` factor.

## Source transcription

Part II defines the p.u.l. ``Z'`` and ``Y'`` matrices.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\widehat\Lambda`` | unchanged | modal geometric operator | source normalized |
| ``S_{1g},S_{2g}`` | unchanged | modal ground field operators | source normalized |
| ``Z'`` | unchanged | p.u.l. series matrix | ``Ω/m`` |

## Evidence and approximation sources

The source gives the required equations directly.

## Limitations and discrepancies

The compact record depends on the paper's preceding modal definitions; it is not standalone numerical code.
