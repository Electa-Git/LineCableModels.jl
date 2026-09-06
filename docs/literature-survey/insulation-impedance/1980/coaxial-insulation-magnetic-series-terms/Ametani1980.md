# Ametani coaxial-insulation series-impedance terms

## Identity and source

| Field | Value |
| --- | --- |
| Family | Insulation impedance |
| Geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. |
| Calculated quantities | Per-unit-length magnetic series contributions of the core–sheath, sheath–armor, and armor–exterior insulation regions |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation within the stated coaxial component model; no expansion or truncation is printed. |
| Main source | A. Ametani (1980) |
| Citation key(s) | `:Ametani1980` |
| Evidence status | PDF page images checked |

**Description.** Per-unit-length series-impedance contributions associated with the three concentric insulation regions of a single-core coaxial cable having core, sheath, and armor.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated; no longitudinal propagation constant occurs in these component expressions. | Unresolved — component list following (13). |
| Air propagation constant ``γ_air`` | Not applicable. | These are cable-internal insulation terms. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Earth-return impedance is the separate ``[Z_o]`` matrix. |
| Earth permittivity and displacement current | Not applicable. | Earth is absent from these terms. |
| Range of validity | No explicit frequency bound is supplied. | Not stated — printed pp. 903–904. |
| Earth permeability ``μ_earth`` | Not applicable. | Earth is absent. |
| Arrangement | Not applicable to overhead/underground placement; the terms enter ``[Z_i]`` independently of ``[Z_o]``. | Stated — decomposition (3), printed p. 902. |
| Earth structure | Not applicable. | Earth-return terms are separate. |
| Conductor and insulation geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. | Stated — Fig. 1(a) and component list, printed pp. 903–904. |
| Constitutive and field assumptions | Scalar relative permeabilities ``\mu_{i1}``, ``\mu_{i2}``, ``\mu_{i3}``; circular coaxial field. No dielectric-loss parameter appears. | Equation-implied — printed expressions. |
| Conventions | ``s=j\omega``; natural logarithm; per-unit-length components. | Stated — (4), printed p. 902, and component heading, p. 903. |

**Expression.** Insulation-region series impedances, component items (2), (6), and (10) following equation (13).

```math
z_{12}=\frac{s\mu_0\mu_{i1}}{2\pi}\ln\!\left(\frac{r_3}{r_2}\right),
\qquad
z_{23}=\frac{s\mu_0\mu_{i2}}{2\pi}\ln\!\left(\frac{r_5}{r_4}\right),

z_{34}=\frac{s\mu_0\mu_{i3}}{2\pi}\ln\!\left(\frac{r_7}{r_6}\right),
\qquad s=j\omega.
```

The source calls these the core outer, sheath outer, and armor outer *insulator impedances*. They enter ``z_{cs}=z_{11}+z_{12}+z_{2i}``, ``z_{sa}=z_{2o}+z_{23}+z_{3i}``, and ``z_{a4}=z_{3o}+z_{34}`` before cable-matrix assembly.

**Approximation.** Not an analytical approximation within the stated coaxial component model; no expansion or truncation is printed.

**Limitations.** Concentric circular insulation only. These are component contributions, not complete cable impedances. The paper gives no anisotropic, eccentric, or frequency-dispersive insulation treatment. Its prose about neglected displacement current is internally ambiguous; the printed expressions are preserved without repair.

**Reference.** [Ametani1980](@cite), component items (2), (6), and (10) following (13), printed pp. 903–904 (PDF pages 2–3), with assembly (9)–(10), p. 903.

**Transcription source.** Original publication; coefficients, radius ratios, layer subscripts, and ``s=j\omega`` were checked visually against the PDF page images. No Markdown conversion was used as authority.

## Source transcription

The expression above retains all three separately printed insulation layers. The associated assembly is

```math
z_{cs}=z_{11}+z_{12}+z_{2i},\qquad
z_{sa}=z_{2o}+z_{23}+z_{3i},\qquad
z_{a4}=z_{3o}+z_{34}.
\tag{10}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``z_{12},z_{23},z_{34}`` | unchanged | Three insulation-region series contributions | per unit length |
| ``r_2,r_3`` | unchanged | Core outer/sheath inner radii | length |
| ``r_4,r_5`` | unchanged | Sheath outer/armor inner radii | length |
| ``r_6,r_7`` | unchanged | Armor outer/outer-insulation radii | length |
| ``\mu_0`` | unchanged | Vacuum permeability | source constant |
| ``\mu_{i1},\mu_{i2},\mu_{i3}`` | unchanged | Insulation relative permeabilities | dimensionless |
| ``s`` | unchanged | Complex-frequency factor | ``s=j\omega`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1`` |

## Evidence and approximation sources

General decomposition and ``s=j\omega``: (3)–(4), printed p. 902. Geometry: Fig. 1(a), p. 903. Assembly and component equations: (8)–(13) and the following numbered list, pp. 903–904.

## Limitations and discrepancies

- The two inspected copies contain the same publication and do not define separate formulations.
- The opening displacement-current assumption and author discussion on p. 910 are internally inconsistent; this record makes no corrective inference.
- No equation-preserving Markdown conversion was located.
