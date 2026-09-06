# Papadopoulos–Tsiamitros–Papagiannis homogeneous-earth underground-cable admittance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; the outermost radius enters the self substitution. Insulation admittance is separate and is combined through the potential-coefficient matrix in Appendix C. |
| Calculated quantities | Per-unit-length mutual earth-return potential coefficient and admittance correction; self term by the source-prescribed substitution |
| Earth structure | Homogeneous earth half-space ``z\geq0`` with air in ``z<0`` and planar interface ``z=0``. |
| Model and approximation | The final spectral expression follows the source's lossless longitudinal-propagation approximation ``\gamma_x\simeq j\omega\sqrt{\mu_1\varepsilon_1}`` and the transform identity (4). No analytical truncation of the semi-infinite integral is stated. ``G(\lambda)`` is retained; dropping it is a different approximation explicitly associated by the source with ``\gamma_x=\gamma_1``. |
| Main source | Theofilos A. Papadopoulos, Dimitrios A. Tsiamitros, and Grigoris K. Papagiannis (2010) |
| Citation key(s) | `:Papadopoulos2010b` |
| Evidence status | Verified visually against the original PDF page image |

**Description.** Per-unit-length mutual earth-return potential coefficient and the source-defined admittance correction between two parallel single-core cables buried in homogeneous conducting earth. The kernel retains radial displacement current through ``G(\lambda)``.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source symbol ``\gamma_x`` with longitudinal factor ``e^{-\gamma_x x}``; for (6) the source sets ``\gamma_x\simeq j k_x=j\omega\sqrt{\mu_1\varepsilon_1}``. | Stated — printed p. 962, text below (3b) and preceding (4). |
| Air propagation constant ``γ_air`` | ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)`` under the indexed definition; air is labelled only by ``\mu_0,\varepsilon_0``, implying ``\sigma_0=0``. ``\gamma_0`` is retained explicitly in ``G`` and through ``\alpha_0``. | Equation-implied — definition below (2b), Fig. 1, and (6c), printed pp. 962–963. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1^2=j\omega\mu_1(\sigma_1+j\omega\varepsilon_1)``; retained explicitly in ``G`` and through ``\alpha_1``. | Stated — definition below (2b) and (6c), printed pp. 962–963. |
| Earth permittivity and displacement current | ``\varepsilon_1`` and the radial-displacement-current term ``G(\lambda)`` are retained. The source states that ignoring ``G`` makes the propagation constant equal to ``\gamma_1``. | Stated — (6b)–(6c) and following paragraph, printed p. 963. |
| Range of validity | Quasi-TEM propagation; no universal frequency or dimensionless bound is supplied. The selected ``\gamma_x`` prescription is described as a better high-frequency approximation than ``\gamma_x=0``. | Stated — printed p. 962, text below (3b); printed p. 968, conclusion. |
| Earth permeability ``μ_earth`` | Arbitrary scalar ``\mu_1`` in (6c); relative permeability is unity only in the numerical examples. | Equation-implied — (6c), printed p. 963. |
| Arrangement | Underground, parallel conductors; mutual term for cables ``i,j`` and self term by replacing ``y_{ij}`` with the outermost cable radius and ``h_j`` with ``h_i``. | Stated — Fig. 2 and text following (6c), printed pp. 962–963. |
| Earth structure | Homogeneous earth half-space ``z\geq0`` with air in ``z<0`` and planar interface ``z=0``. | Stated — Fig. 1 and headings above (1)–(2), printed p. 962. |
| Conductor and insulation geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; the outermost radius enters the self substitution. Insulation admittance is separate and is combined through the potential-coefficient matrix in Appendix C. | Stated — Section III and Appendix C, printed pp. 962–963 and 967–968. |
| Constitutive and field assumptions | Homogeneous scalar ``\mu_k,\varepsilon_k,\sigma_k`` in each medium; quasi-TEM field propagation. | Equation-implied — definitions below (2b), printed p. 962; stated in conclusion, printed p. 968. |
| Conventions | ``j`` is the imaginary unit; ``z`` is positive downward; longitudinal dependence ``e^{-\gamma_x x}``; per-unit-length quantities carry a prime. The propagation definition is consistent with ``e^{j\omega t}``, but the inspected section does not explicitly state the time convention. | Stated/equation-implied — Fig. 1, definition below (2b), and (3), printed p. 962. |

**Expression.** Source-defined per-unit-length mutual earth admittance and its earth-return potential coefficient, equations (6a)–(6c).

```math
Y'_{e_{ij}}=j\omega P_{e_{ij}}^{-1},
\qquad\text{(6a)}

P_{e_{ij}}
=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
 \int_{0}^{+\infty}[F(\lambda)+G(\lambda)]\cos(y_{ij}\lambda)\,d\lambda,
\qquad\text{(6b)}

F(\lambda)
=\frac{e^{-\alpha_1|h_i-h_j|}-e^{-\alpha_1(h_i+h_j)}}{\alpha_1}
 +\frac{2\mu_0e^{-\alpha_1(h_i+h_j)}}{\alpha_1\mu_0+\alpha_0\mu_1},

G(\lambda)
=\frac{2\mu_0\mu_1\alpha_1(\gamma_1^2-\gamma_0^2)e^{-\alpha_1(h_i+h_j)}}
 {(\alpha_1\mu_0+\alpha_0\mu_1)
  (\alpha_1\gamma_0^2\mu_1+\alpha_0\gamma_1^2\mu_0)},
\qquad\text{(6c)}

\alpha_k=\sqrt{\lambda^2+\gamma_k^2+k_x^2},\qquad
\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k),\qquad k=0,1,

\gamma_x\simeq j k_x=j\omega\sqrt{\mu_1\varepsilon_1}.
```

``P_{e_{ij}}`` is the mutual earth-return potential coefficient and ``Y'_{e_{ij}}`` is the corresponding per-unit-length admittance as printed. Indices ``0`` and ``1`` denote air and earth. The source assembles the full cable shunt-admittance matrix from the internal-insulation and earth-return potential-coefficient matrices in Appendix C. For the self term of cable ``i``, it prescribes ``y_{ij}\mapsto r_{i,\mathrm{outer}}`` and ``h_j\mapsto h_i``.

**Approximation.** The final spectral expression follows the source's lossless longitudinal-propagation approximation ``\gamma_x\simeq j\omega\sqrt{\mu_1\varepsilon_1}`` and the transform identity (4). No analytical truncation of the semi-infinite integral is stated. ``G(\lambda)`` is retained; dropping it is a different approximation explicitly associated by the source with ``\gamma_x=\gamma_1``.

**Limitations.** The expression assumes infinitely long parallel cables, a homogeneous planar earth, and quasi-TEM propagation. It provides an earth correction/potential coefficient, not the insulation admittance or the already-assembled total cable admittance. Equation (6a) is preserved in the source's printed scalar notation; this record does not replace it with an inferred matrix operation. The Markdown conversion drops the equation contents.

**Reference.** [Papadopoulos2010b](@cite), equations (6a)–(6c), printed p. 963 (PDF page 3), with definitions and approximation on printed p. 962 (PDF page 2) and matrix assembly in Appendix C, printed pp. 967–968.

**Transcription source.** Original publication. Equations, signs, separate ``F`` and ``G`` kernels, limits, powers, and denominators were checked visually against the source file, printed pp. 962–963. The Docling Markdown was a discovery aid only.

## Source transcription

The formula section expression preserves the source notation and order of (6a), (6b), and (6c). ``F`` is repeated from (5b) so this record is independently readable. The source explicitly identifies ``G`` with radial displacement currents and gives the self substitutions after (6c). Appendix C places the earth-return potential coefficients into a matrix before total-admittance assembly; this record does not infer a replacement for the source-printed (6a).

## Notation map

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Y'_{e_{ij}}`` | unchanged | Mutual earth-return admittance correction | ``\mathrm{S}/\mathrm{m}``; prime denotes per-unit-length quantity |
| ``P_{e_{ij}}`` | unchanged | Mutual earth-return potential coefficient | Units not stated beside (6); source-defined through (6a)–(6b) |
| ``F(\lambda)`` | unchanged | Spectral kernel also used by the impedance correction | Source-defined by (5b) |
| ``G(\lambda)`` | unchanged | Radial-displacement-current spectral kernel | Source-defined by (6c) |
| ``\lambda`` | unchanged | Transformed spectral integration variable | ``\mathrm{m}^{-1}``; ``0`` to ``+\infty`` |
| ``\alpha_k`` | unchanged | Vertical spectral propagation factor after transformation | ``\mathrm{m}^{-1}``; square-root branch not stated explicitly |
| ``\gamma_k`` | unchanged | Bulk propagation constant of medium ``k`` | ``\mathrm{m}^{-1}``; ``k=0`` air, ``k=1`` earth |
| ``\gamma_x`` | unchanged | Imposed longitudinal propagation constant | ``\mathrm{m}^{-1}``; dependence ``e^{-\gamma_x x}`` |
| ``k_x`` | unchanged | Lossless longitudinal wavenumber | ``\mathrm{m}^{-1}`` |
| ``h_i,h_j`` | unchanged | Positive burial depths | ``\mathrm{m}`` |
| ``y_{ij}`` | unchanged | Horizontal cable-axis separation | ``\mathrm{m}`` |
| ``\mu_0,\mu_1`` | unchanged | Air and earth permeability | ``\mathrm{H}/\mathrm{m}`` |
| ``\varepsilon_0,\varepsilon_1`` | unchanged | Air and earth permittivity | ``\mathrm{F}/\mathrm{m}`` |
| ``\sigma_0,\sigma_1`` | unchanged | Air and earth conductivity in the indexed propagation definition | ``\mathrm{S}/\mathrm{m}``; air conductivity equation-implied zero |
| ``\omega`` | unchanged | Angular frequency | ``\mathrm{rad}/\mathrm{s}`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1`` |

## Evidence and approximation sources

- Media, coordinates, and bulk propagation definitions: Fig. 1 and paragraph below (2b), printed p. 962.
- Longitudinal dependence and approximation: (3b), text below it, and (4), printed p. 962.
- Mutual potential coefficient, admittance relation, and ``F/G`` kernels: (6a)–(6c), printed p. 963.
- Radial-displacement interpretation and self substitution: paragraph following (6c), printed p. 963.
- Matrix assembly with insulation and earth potential coefficients: Appendix C, printed pp. 967–968.

## Limitations and discrepancies

- Conversion defect: the Markdown conversion preserves the labels but loses the mathematical contents of (3)–(6).
- The source does not state the square-root branch for ``\alpha_k`` in the inspected definition.
- Equation (6a) is printed as an inverse of ``P_{e_{ij}}``; the later matrix assembly must be consulted before treating elementwise inversion as a general multiconductor matrix rule.
- The high-frequency longitudinal approximation has no author-stated universal bound.
