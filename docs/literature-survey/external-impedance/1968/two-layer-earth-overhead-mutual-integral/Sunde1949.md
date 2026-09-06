# Sunde two-layer-earth overhead mutual inductance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary external wire positions; heights ``h_1,h_2`` and separation ``y``; finite conductor radius and insulation do not enter. |
| Calculated quantities | Mutual external inductance/impedance kernel for infinite overhead wires above two-layer earth; large-separation two-term expansion |
| Earth structure | Finite upper layer of depth ``d`` over a semi-infinite lower layer. |
| Model and approximation | Equation (4.55) already neglects the direct-distance logarithmic term for large separation and replaces horizontal by radial separation within the stated practical range. Equation (4.56) then expands ``F(u)`` as in the uniform-earth case and retains only its first two terms. |
| Main source | E. D. Sunde; formula appears in a 1968 corrected Dover republication of a work first published in 1949, but first-edition token identity was not inspected |
| Citation key(s) | `:Sunde1949` |
| Evidence status | Corrected-republication page image verified; original-1949 priority and token identity unresolved |

**Description.** Mutual inductance per cycle for two infinitely long parallel wires above a two-layer conductive earth, with a finite upper layer and a semi-infinite lower layer, plus the source's large-separation two-term expansion.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed low-frequency earth-return reduction; no independent longitudinal constant remains in (4.47) or (4.55). | Equation-implied — §§4.7–4.8, pp. 114–119. |
| Air propagation constant ``γ_air`` | Omitted from the displayed low-frequency kernel. | Equation-implied — (4.47), (4.55). |
| Earth propagation constant ``γ_earth`` | Layer roots enter as ``\alpha_i``; adjacent chapter notation identifies intrinsic layer constants ``\gamma_i`` but the defining equation was not reverified in this pass. | Unresolved definition dependency — (4.47), p. 114 and (4.55), p. 119. |
| Earth permittivity and displacement current | Not present in the inspected chapter formulas; the specific omission statement is not repeated on pp. 114–119. | Equation-implied — (4.47)–(4.56); exact earlier definition unresolved. |
| Range of validity | Equation (4.55) is introduced for practical ranges using radial separation in place of horizontal separation; the logarithmic direct-distance term is neglected for large separations. Equation (4.56) retains the first two expansion terms. No numerical universal bound is printed beside (4.55). | Stated — §4.8, p. 119. |
| Earth permeability ``μ_earth`` | Not explicit in (4.47)–(4.56); layer-root definition needed. | Unresolved — inspected pages. |
| Arrangement | Overhead/overhead, mutual, infinite parallel wires. | Stated — heading and §4.8, p. 119. |
| Earth structure | Finite upper layer of depth ``d`` over a semi-infinite lower layer. | Stated — §§4.7–4.8 and exponentials in (4.47)/(4.55). |
| Conductor and insulation geometry | Filamentary external wire positions; heights ``h_1,h_2`` and separation ``y``; finite conductor radius and insulation do not enter. | Stated — (4.55) and surrounding text. |
| Constitutive and field assumptions | Horizontally uniform isotropic layers in the low-frequency earth-return model. | Stated/equation-implied — chapter context and (4.47). |
| Conventions | Source uses ``i`` for the imaginary unit elsewhere in the chapter and reports mutual inductance ``L`` per cycle with factor ``v/\pi``; ``v`` is the source electromagnetic unit-conversion constant. SI conversion is not imposed. | Stated — equations and chart units, pp. 114–119. |

**Expression.** Two-layer overhead-wire mutual inductance, equation (4.55), printed p. 119.

```math
L=\frac{v}{\pi}\int_0^\infty
F(u)e^{-(h_1+h_2)u}\cos(uy)\,du,
\qquad\text{(4.55)}
```

```math
F(u)=
\frac{\alpha_1+\alpha_2+(\alpha_1-\alpha_2)e^{-2\alpha_1d}}
{(\alpha_1+\alpha_2)(u+\alpha_1)+(\alpha_1-\alpha_2)(u-\alpha_1)e^{-2d\alpha_1}}.
```

For large spacing the source expands ``F`` and retains two terms:

```math
L=\frac{v}{\pi}\left[
F(0)\frac{\cos\theta}{r'_{12}}+F'(0)\frac{\cos2\theta}{(r'_{12})^2}
\right]
=\frac{v}{\pi}\left[
\frac{\cos\theta}{\gamma'_1r'_{12}}-
\frac{\cos2\theta}{(\gamma'_1r'_{12})^2}
\right].
\qquad\text{(4.56)}
```

**Approximation.** Equation (4.55) already neglects the direct-distance logarithmic term for large separation and replaces horizontal by radial separation within the stated practical range. Equation (4.56) then expands ``F(u)`` as in the uniform-earth case and retains only its first two terms.

**Limitations.** This record cannot establish whether the 1949 first edition printed identical formulas. Definitions of ``\alpha_i``, ``\gamma'_1``, ``v``, ``r'_{12}``, and ``\theta`` are chapter dependencies not repeated on p. 119; their exact earlier tokens remain unresolved and the formula is therefore verification-complete only at the displayed-equation level, not dependency-complete. Displacement current, permeability treatment and SI conversion are not guessed.

**Reference.** [Sunde1949](@cite), corrected 1968 Dover republication, §4.8, equations (4.55)–(4.56), printed p. 119 (PDF page 134); edition conflict applies.

**Transcription source.** Original page images of the accessible corrected republication. Numerator/denominator signs, both ``e^{-2\alpha_1d}`` factors, height exponential, cosine, expansion orders and prime marks were visually checked. This is not represented as verification against an unseen 1949 copy.

## Source transcription

For wires on the layer surface, the immediately preceding source prints

```math
L=\frac{v}{\pi}\int_0^\infty
\frac{(\alpha_1+\alpha_2)+(\alpha_1-\alpha_2)e^{-2d\alpha_1}}
{(\alpha_1+\alpha_2)(u+\alpha_1)+(\alpha_1-\alpha_2)(u-\alpha_1)e^{-2d\alpha_1}}
\cos(uy)\,du.
\qquad\text{(4.47)}
```

Equations (4.48)–(4.54) are restricted thin-layer and material-contrast limits; they are not merged into (4.55).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``L`` | unchanged | source mutual inductance per cycle | source electromagnetic units |
| ``u`` | unchanged | horizontal spectral variable | inverse length |
| ``\alpha_1,\alpha_2`` | unchanged | layer spectral roots | exact earlier definition unresolved |
| ``d`` | unchanged | upper-layer thickness | source length unit |
| ``h_1,h_2`` | unchanged | wire heights | source length unit |
| ``y`` | unchanged | source separation coordinate | source length unit |
| ``F(u)`` | unchanged | two-layer interface kernel | source normalized |
| ``r'_{12},\theta,\gamma'_1`` | unchanged | polar separation variables and normalized layer constant | earlier chapter definitions unresolved |
| ``v`` | unchanged | electromagnetic unit-conversion constant | not silently converted to SI |

No notation was renamed.

## Evidence and approximation sources

The surface kernel (4.47) supplies the interface denominator. Section 4.8 inserts the height exponential and then states the large-separation simplification and two-term expansion. Pages 115–118 contain empirical charts and thin-layer special cases that inform limits but do not replace the main formula.

## Limitations and discrepancies

- Sections 4.1–4.7 and the later cable chapters provide related field solutions and special cases. This record is restricted to (4.55)–(4.56).
