# Déri–Tevan–Semlyen–Castanheira homogeneous complex ground-return plane

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Ideal, infinitely long parallel thin wires represented by radius, heights, and separation; conductor internal impedance and insulation are outside the expression. |
| Calculated quantities | Per-unit-length self and mutual ideal-conductor/ground-return loop impedances for overhead wires |
| Earth structure | Homogeneous conductive half-space under a plane interface. |
| Model and approximation | The authors transform Carson's correction integral and introduce the explicit kernel approximation in their equation (37), then obtain (3) and (4). Thus the closed logarithms are approximations to Carson, not exact full-wave expressions. Their proof of the mutual result initially requires ``β<1``; the wider practical range is supported by numerical testing. |
| Main source | Equations originally proposed by Dubanton and published by Gary; heuristic justification, analytical relation to Carson, and numerical error evaluation by A. Déri, G. Tevan, A. Semlyen, and A. Castanheira (1981) |
| Citation key(s) | Primary attribution: `:Dubanton1969`; equation-complete English derivation: `:Deri1981`; auxiliary complex-image derivation: `:Wait1969` |
| Evidence status | PDF page images checked for the English derivation; DOI agrees with its primary-PDF metadata |

**Description.** A homogeneous conducting earth is replaced by a perfectly conducting image plane at the complex penetration depth ``p``. The resulting closed logarithmic self and mutual impedances approximate Carson's infinite-integral expressions over Carson's own range of validity.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not defined or retained. The formulation is explicitly related to Carson's line-return equations under the long-wavelength condition stated by the authors. | Stated — Introduction, printed p. 3686. |
| Air propagation constant ``γ_air`` | Not defined; the air region contributes only through geometry and ``μ_0``. | Equation-implied — (3)–(4), printed p. 3686. |
| Earth propagation constant ``γ_earth`` | No independently named propagation constant; conduction enters through complex depth ``p=1/\sqrt{j\omega\mu_0\sigma}``. | Stated — (18), printed p. 3688. |
| Earth permittivity and displacement current | Earth capacitive displacement current is neglected. | Stated — Introduction, printed p. 3686. |
| Range of validity | Restricted to Carson's regime: frequency low enough to neglect earth displacement current and wavelength large relative to transverse geometry. Numerical comparisons report total-impedance magnitude error below 3% for ``β\leq2`` and about 0.5% or less for typical ``β<0.5`` geometries; these are tested errors, not universal bounds. | Stated — printed pp. 3686, 3691–3692. |
| Earth permeability ``μ_earth`` | ``\mu=\mu_0`` for the combined expressions (3)–(4). The paper treats arbitrary ``\mu`` only through separate ideal-ground and correction inductance terms. | Stated — text before (3) and discussion around (28)–(32), printed pp. 3686, 3689. |
| Arrangement | Overhead; self term for height ``h`` and radius ``r``; mutual term for parallel wires at heights ``h_k,h_\ell`` and horizontal separation ``d_{k\ell}``. | Stated — (3)–(4), printed p. 3686. |
| Earth structure | Homogeneous conductive half-space under a plane interface. | Stated — Introduction and (18), printed pp. 3686, 3688. |
| Conductor and insulation geometry | Ideal, infinitely long parallel thin wires represented by radius, heights, and separation; conductor internal impedance and insulation are outside the expression. | Stated/equation-implied — transmission-line geometry and (3)–(4), printed p. 3686. |
| Constitutive and field assumptions | Linear harmonic plane-field penetration in homogeneous earth with scalar conductivity ``\sigma`` and permeability ``\mu_0``; complex plane retained heuristically as an image surface for thin conductors. | Stated — printed pp. 3687–3688. |
| Conventions | ``j`` is the imaginary unit; ``\omega=2\pi f``; depth is positive downward; (9′) uses ``dE/dx=-j\omega\mu_0H``, consistent with ``e^{j\omega t}``. Natural logarithms are used. | Stated/equation-implied — (3)–(4), (9′), (18)–(19), printed pp. 3686–3688. |

**Expression.** Source equations (3), (4), (18), (18′), and (19).

```math
Z_s=j\omega\frac{\mu_0}{2\pi}\ln\frac{2(h+p)}{r},
\qquad\text{(3)}

Z_m=j\omega\frac{\mu_0}{2\pi}\ln
\frac{\sqrt{(h_k+h_\ell+2p)^2+d_{k\ell}^2}}
     {\sqrt{(h_k-h_\ell)^2+d_{k\ell}^2}},
\qquad\text{(4)}

p=\frac{1}{\sqrt{j\omega\mu_0\sigma}},
\qquad\text{(18)}

\delta=\frac{1}{\sqrt{\pi f\mu_0\sigma}},
\qquad\text{(18')}

\frac{1}{p}=(1+j)\frac{1}{\delta}.
\qquad\text{(19)}
```

``Z_s`` is the self loop impedance and ``Z_m`` the mutual impedance. Equation (19) fixes the complex square-root choice used by the source. The image of a conductor is displaced by ``2p`` below the real surface.

**Approximation.** The authors transform Carson's correction integral and introduce the explicit kernel approximation in their equation (37), then obtain (3) and (4). Thus the closed logarithms are approximations to Carson, not exact full-wave expressions. Their proof of the mutual result initially requires ``β<1``; the wider practical range is supported by numerical testing.

**Limitations.** Homogeneous, nonmagnetic, plane earth and overhead parallel wires only. Earth displacement current and longitudinal full-wave effects are excluded. The expressions are ideal-conductor/ground-return loop terms and do not supply conductor internal impedance. Reported error curves depend on the tested dimensionless geometry and frequency parameters.

**Reference.** [Deri1981](@cite), equations (3)–(4), printed p. 3686; equations (18)–(19), printed p. 3688; error evaluation on printed pp. 3691–3692.

**Transcription source.** Original publication PDF. Every radical, height sum/difference, factor ``2p``, logarithmic ratio, and penetration-depth factor was checked against the rendered page images. The Docling Markdown was used only for navigation because it corrupts equations (1)–(24).

## Source transcription

For the numerical error comparison the source defines

```math
\beta=\frac{d_{k\ell}}{h_k+h_\ell},
\qquad\text{(44)}

\alpha=\frac{h_{\mathrm{ave}}}{\delta}
=\frac{h_k+h_\ell}{2\sqrt{2}}\sqrt{\omega\mu_0\sigma}.
\qquad\text{(48--48')}
```

The reported ``3%`` and ``0.5%`` figures apply only to the source's comparisons against Carson's evaluation in the tested cases.

## Notation map

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_s`` | unchanged | Self ideal-conductor/ground-return loop impedance | per unit length |
| ``Z_m`` | unchanged | Mutual ground-return impedance | per unit length |
| ``p`` | unchanged | Complex depth of the image/return plane | length; branch fixed by (19) |
| ``\delta`` | unchanged | Real penetration depth | length |
| ``h,h_k,h_\ell`` | unchanged | Conductor heights above earth | length |
| ``r`` | unchanged | Conductor radius | length |
| ``d_{k\ell}`` | unchanged | Horizontal separation of conductors ``k,\ell`` | length |
| ``\sigma`` | unchanged | Homogeneous-earth conductivity | ``\mathrm{S/m}`` |
| ``\mu_0`` | unchanged | Free-space permeability, also assigned to earth | ``\mathrm{H/m}`` |
| ``\omega,f`` | unchanged | Angular frequency and frequency | ``\omega=2\pi f`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1`` |
| ``\alpha,\beta`` | unchanged | Dimensionless frequency/height and separation ratios used in the error study | dimensionless |

## Evidence and approximation sources

- Carson-regime assumptions and original Dubanton/Gary attribution: printed p. 3686.
- Closed self/mutual expressions: (3)–(4), printed p. 3686.
- Plane-field equations and homogeneous penetration depth: (9), (16)–(19), printed pp. 3687–3688.
- Analytical approximation step: (35)–(39), printed pp. 3689–3690; mutual counterpart (43)–(47), printed p. 3691.
- Numerical-error variables and results: (44), (48)–(49), Fig. 7, printed p. 3691; Fig. 8 and Conclusions, printed p. 3692.

## Limitations and discrepancies

- The formula is often called “Deri–Semlyen,” but the source explicitly attributes equations (3)–(4) to Dubanton and Gary; the four 1981 authors supply the justification, derivation, testing, and layered extension.
- The source's phrase “whole range of frequencies” is bounded by the stated Carson regime; it is not a full-wave all-frequency claim.
