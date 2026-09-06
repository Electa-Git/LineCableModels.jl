# Ametani lossless coaxial-insulation potential coefficients

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. |
| Calculated quantities | Core, sheath, and armor insulation potential coefficients and admittance-matrix assembly |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation within the concentric, lossless dielectric model. Omission of dielectric loss is a constitutive restriction, not a corpus modification. |
| Main source | A. Ametani (1980) |
| Citation key(s) | `:Ametani1980` |
| Evidence status | PDF page images checked |

**Description.** Lossless coaxial-insulation potential coefficients for the core–sheath, sheath–armor, and armor–exterior regions of a single-core cable, with the source's matrix conversion to per-unit-length shunt admittance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated; absent from the potential-coefficient expressions. | Unresolved — (4), (17)–(24). |
| Air propagation constant ``γ_air`` | Not applicable to ``[P_i]``. | The separate ``[P_o]`` handles air. |
| Earth propagation constant ``γ_earth`` | Not applicable. | These are cable-internal coefficients. |
| Earth permittivity and displacement current | Not applicable. | Earth is absent from ``P_{cj},P_{sj},P_{aj}``. |
| Range of validity | No explicit frequency bound; the printed insulation model contains permittivity but no dielectric-loss term. | Equation-implied — (22); dielectric losses declared negligible on p. 902. |
| Earth permeability ``μ_earth`` | Not applicable. | Earth is absent. |
| Arrangement | Not applicable to cable placement for ``[P_i]``; the outer-medium matrix changes with placement. | Stated — (17)–(20), p. 904. |
| Earth structure | Not applicable. | These are cable-internal terms. |
| Conductor and insulation geometry | Concentric annuli ``r_2<r<r_3``, ``r_4<r<r_5``, and ``r_6<r<r_7``. | Stated — Fig. 1(a), p. 903, and (22), p. 904. |
| Constitutive and field assumptions | Scalar relative permittivities ``\varepsilon_{i1}``, ``\varepsilon_{i2}``, ``\varepsilon_{i3}``; no conductivity/loss tangent, hence a lossless printed model. | Equation-implied — (22), p. 904. |
| Conventions | ``s=j\omega``; ``[Y]=s[P]^{-1}``; natural logarithm; per-unit-length matrices. | Stated — (4), p. 902. |

**Expression.** Internal potential coefficients and source-prescribed admittance conversion, equations (4), (21), and (22).

```math
[Y]=s[P]^{-1},\qquad s=j\omega,
\qquad\text{(4)}

P_{cj}=\frac{1}{2\pi\varepsilon_0\varepsilon_{i1}}\ln\!\left(\frac{r_3}{r_2}\right),
\quad
P_{sj}=\frac{1}{2\pi\varepsilon_0\varepsilon_{i2}}\ln\!\left(\frac{r_5}{r_4}\right),

P_{aj}=\frac{1}{2\pi\varepsilon_0\varepsilon_{i3}}\ln\!\left(\frac{r_7}{r_6}\right),
\qquad\text{(22)}

[P_{ij}]=
\begin{bmatrix}
P_{cj}+P_{sj}+P_{aj} & P_{sj}+P_{aj} & P_{aj}\\
P_{sj}+P_{aj} & P_{sj}+P_{aj} & P_{aj}\\
P_{aj} & P_{aj} & P_{aj}
\end{bmatrix}.
\qquad\text{(21)}
```

Rows/columns correspond to core, sheath, and armor. For core-and-sheath, (23) is the leading ``2\times2`` form; for core-only, (24) gives ``[P_{ij}]=P_{cj}``. ``[Y]`` is obtained only after assembling the applicable potential-coefficient matrices.

**Approximation.** Not an analytical approximation within the concentric, lossless dielectric model. Omission of dielectric loss is a constitutive restriction, not a corpus modification.

**Limitations.** Concentric circular, scalar, lossless insulation only. No conductivity, complex permittivity, loss tangent, dispersion, eccentricity, or anisotropy. ``[P_{ij}]`` is not admittance; the source-prescribed matrix inverse must not be replaced by elementwise inversion.

**Reference.** [Ametani1980](@cite), equations (4), (21)–(24), printed pp. 902 and 904 (PDF pages 1 and 3).

**Transcription source.** Original publication; radius ratios, permittivity subscripts, matrix entries, and the matrix inverse were visually checked against the PDF pages.

## Source transcription

The source additionally prints

```math
[P_{ij}]=\begin{bmatrix}P_{cj}+P_{sj}&P_{sj}\\P_{sj}&P_{sj}\end{bmatrix}
\qquad\text{(23)}
```

for a core-and-sheath cable and ``[P_{ij}]=P_{cj}`` in (24) for core-only. For an underground cable it states ``[P]=[P_i]`` in (18); that outer-medium choice is not used to redefine the insulation coefficients.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``[Y]`` | unchanged | Cable shunt-admittance matrix per unit length | ``s[P]^{-1}`` |
| ``[P]`` | unchanged | Total potential-coefficient matrix | matrix inverse required |
| ``[P_{ij}]`` | unchanged | Cable ``j`` internal submatrix | core/sheath/armor ordering |
| ``P_{cj},P_{sj},P_{aj}`` | unchanged | Three insulation potential coefficients | source does not annotate units beside (22) |
| ``r_2,r_3,r_4,r_5,r_6,r_7`` | unchanged | Successive conductor/insulation radii | length |
| ``\varepsilon_0`` | unchanged | Vacuum permittivity | source constant |
| ``\varepsilon_{i1},\varepsilon_{i2},\varepsilon_{i3}`` | unchanged | Insulation relative permittivities | dimensionless |
| ``s`` | unchanged | Complex-frequency factor | ``s=j\omega`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1`` |

## Evidence and approximation sources

Admittance conversion: (4), p. 902. Geometry: Fig. 1(a), p. 903. Placement separation: (17)–(20), p. 904. Three-layer and reduced matrices: (21)–(24), p. 904.

## Limitations and discrepancies

- The two PDFs are duplicate witnesses, not distinct formulations.
- The opening displacement-current statement and author reply on p. 910 are internally inconsistent. The printed insulation coefficients retain ``s\varepsilon`` behavior while omitting dielectric loss; no repair is made.
- No equation-preserving Markdown conversion was located.

