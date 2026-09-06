# Papadopoulos–Tsiamitros–Papagiannis homogeneous-earth underground-cable impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; depths ``h_i,h_j`` and horizontal separation ``y_{ij}``; the outermost cable radius is used in the stated self substitution. Internal conductor and insulation terms are assembled separately. |
| Calculated quantities | Per-unit-length mutual earth-return impedance correction; self term by the source-prescribed substitution |
| Earth structure | Homogeneous earth half-space ``z\geq0`` with air in ``z<0`` and a planar interface at ``z=0``. |
| Model and approximation | Equation (5) applies the lossless longitudinal-propagation approximation ``\gamma_x\simeq j\omega\sqrt{\mu_1\epsilon_1}`` to the dipole-field starting solution and uses ``u^2-k_x^2=\lambda^2`` with identity (4). No series truncation is applied; the semi-infinite integral is evaluated numerically. |
| Main source | Theofilos A. Papadopoulos, Dimitrios A. Tsiamitros, and Grigoris K. Papagiannis (2010) |
| Citation key(s) | `:Papadopoulos2010b` |
| Evidence status | Verified visually against the original PDF page image |

**Description.** Per-unit-length mutual earth-return series-impedance correction between two parallel single-core cables buried at depths ``h_i`` and ``h_j``, with horizontal separation ``y_{ij}``, in a homogeneous conducting earth below air. The source obtains the self correction by a stated geometric substitution.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source symbol ``\gamma_x`` with longitudinal factor ``e^{-\gamma_x x}``; for the displayed final expression the source sets ``\gamma_x \simeq j k_x = j\omega\sqrt{\mu_1\epsilon_1}``. | Stated — printed p. 962, text below (3b) and preceding (4). |
| Air propagation constant ``γ_air`` | ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\epsilon_0)`` under the source's indexed definition; air is labelled only by ``\mu_0,\epsilon_0``, implying ``\sigma_0=0``. It is retained through ``\alpha_0``. | Equation-implied — printed p. 962, definition below (2b), Fig. 1, and (5b). |
| Earth propagation constant ``γ_earth`` | ``\gamma_1^2=j\omega\mu_1(\sigma_1+j\omega\epsilon_1)``; retained through ``\alpha_1``. | Stated — printed p. 962, definition below (2b). |
| Earth permittivity and displacement current | ``\epsilon_1`` is retained in ``\gamma_1``. | Stated — printed pp. 962–963, propagation-constant definition and discussion following (6c). |
| Range of validity | Quasi-TEM propagation; the source gives no universal frequency or dimensionless bound for (5). The ``\gamma_x\simeq j k_x`` choice is presented as a better high-frequency approximation than ``\gamma_x=0``. | Stated — printed p. 962, text below (3b); printed p. 968, conclusion. |
| Earth permeability ``μ_earth`` | Arbitrary scalar ``\mu_1`` in the expression; relative permeability is set to unity only in the numerical examples. | Equation-implied — (5a)–(5b), printed p. 963; example restriction stated on the same page. |
| Arrangement | Underground, parallel conductors; mutual term for cables ``i,j`` and self term by replacing ``y_{ij}`` with the outermost cable radius and ``h_j`` with ``h_i``. | Stated — Fig. 2 and text following (6c), printed pp. 962–963. |
| Earth structure | Homogeneous earth half-space ``z\geq0`` with air in ``z<0`` and a planar interface at ``z=0``. | Stated — Fig. 1 and headings above (1)–(2), printed p. 962. |
| Conductor and insulation geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; depths ``h_i,h_j`` and horizontal separation ``y_{ij}``; the outermost cable radius is used in the stated self substitution. Internal conductor and insulation terms are assembled separately. | Stated — Section III, Fig. 2, and text following (6c), printed pp. 962–963. |
| Constitutive and field assumptions | Homogeneous scalar ``\mu_k,\epsilon_k,\sigma_k`` in each medium; quasi-TEM field propagation. The earth-return correction does not include conductor skin-effect or insulation contributions. | Equation-implied — definitions below (2b) and matrix-assembly paragraph following (6c), printed pp. 962–963. |
| Conventions | ``j`` is the imaginary unit; ``z`` is positive downward; per-unit-length output; longitudinal dependence ``e^{-\gamma_x x}``. The sign in ``\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\epsilon_k)`` is consistent with an ``e^{j\omega t}`` phasor convention, but the inspected section does not state that convention explicitly. | Stated/equation-implied — Fig. 1, definition below (2b), and (3), printed p. 962. |

**Expression.** Per-unit-length mutual earth-return impedance correction, source equations (5a)–(5b).

```math
Z'_{e_{ij}}
=\frac{j\omega\mu_1}{2\pi}
 \int_{0}^{+\infty} F(\lambda)\cos(y_{ij}\lambda)\,d\lambda,
\tag{5a}

F(\lambda)
=\frac{e^{-\alpha_1|h_i-h_j|}-e^{-\alpha_1(h_i+h_j)}}{\alpha_1}
 +\frac{2\mu_0e^{-\alpha_1(h_i+h_j)}}{\alpha_1\mu_0+\alpha_0\mu_1},
\tag{5b}

\alpha_k=\sqrt{\lambda^2+\gamma_k^2+k_x^2},\qquad
\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\epsilon_k),\qquad k=0,1,

\gamma_x\simeq j k_x=j\omega\sqrt{\mu_1\epsilon_1}.
```

``Z'_{e_{ij}}`` is an earth-return correction in ``\Omega/\mathrm{m}``. Indices ``0`` and ``1`` denote air and earth. For the self term of cable ``i``, the source prescribes ``y_{ij}\mapsto r_{i,\mathrm{outer}}`` and ``h_j\mapsto h_i``. This record retains the two addends of ``F`` separately, as printed.

**Approximation.** The dipole-field starting point is the full electromagnetic-field solution, but the displayed (5) follows the explicit lossless longitudinal-propagation approximation ``\gamma_x\simeq j\omega\sqrt{\mu_1\epsilon_1}`` and the substitution ``u^2-k_x^2=\lambda^2`` using identity (4). No series truncation is applied to (5); its semi-infinite integral is evaluated numerically.

**Limitations.** The source supplies no formal universal frequency bound for the longitudinal approximation. The conductor is infinitely long and parallel to a planar interface, the earth is homogeneous, and the expression is only the earth correction—not the total cable series impedance. The Markdown conversion omits the mathematical content of (3)–(6) and cannot verify this transcription.

**Reference.** [Papadopoulos2010b](@cite), equations (5a)–(5b), printed p. 963 (PDF page 3), with definitions and approximation on printed p. 962 (PDF page 2).

**Transcription source.** Original publication. Equations, absolute value, signs, separate addends, limits, and denominators were checked visually against the source file, printed pp. 962–963. The Docling Markdown was used only to locate the section and is not the mathematical authority.

## Source transcription

The formula section expression preserves the paper's notation and source order for (5a), (5b), and the definitions immediately preceding them. The source then states the self-term substitutions ``y_{ij}\mapsto r_{i,\mathrm{outer}}`` and ``h_j\mapsto h_i``. It further states that setting ``k_x=0`` reproduces Sunde's earth-impedance formula and additionally neglecting earth permittivity reduces it to Pollaczek's formula; those attributed limiting witnesses are not merged into this formulation.

## Notation map

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z'_{e_{ij}}`` | unchanged | Mutual earth-return series-impedance correction | ``\Omega/\mathrm{m}``; prime denotes per-unit-length quantity |
| ``F(\lambda)`` | unchanged | Spectral impedance kernel | Source-defined by (5b) |
| ``\lambda`` | unchanged | Transformed spectral integration variable | ``\mathrm{m}^{-1}``; ``0`` to ``+\infty`` |
| ``\alpha_k`` | unchanged | Vertical spectral propagation factor after the ``u^2-k_x^2`` substitution | ``\mathrm{m}^{-1}``; square-root branch not stated explicitly |
| ``\gamma_k`` | unchanged | Bulk propagation constant of medium ``k`` | ``\mathrm{m}^{-1}``; ``k=0`` air, ``k=1`` earth |
| ``\gamma_x`` | unchanged | Imposed longitudinal propagation constant | ``\mathrm{m}^{-1}``; dependence ``e^{-\gamma_x x}`` |
| ``k_x`` | unchanged | Lossless longitudinal wavenumber used by the approximation | ``\mathrm{m}^{-1}`` |
| ``h_i,h_j`` | unchanged | Positive burial depths | ``\mathrm{m}``; measured in positive ``z`` direction |
| ``y_{ij}`` | unchanged | Horizontal cable-axis separation | ``\mathrm{m}`` |
| ``\mu_0,\mu_1`` | unchanged | Air and earth permeability | ``\mathrm{H}/\mathrm{m}`` |
| ``\epsilon_0,\epsilon_1`` | unchanged | Air and earth permittivity | ``\mathrm{F}/\mathrm{m}`` |
| ``\sigma_0,\sigma_1`` | unchanged | Air and earth conductivity in the indexed propagation definition | ``\mathrm{S}/\mathrm{m}``; air conductivity is equation-implied zero |
| ``\omega`` | unchanged | Angular frequency | ``\mathrm{rad}/\mathrm{s}`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1`` |

## Evidence and approximation sources

- Geometry and media: Fig. 1 and Fig. 2, printed p. 962.
- Bulk and spectral propagation definitions: paragraph below (2b), printed p. 962, and paragraph introducing (5), printed p. 963.
- Longitudinal approximation and transform identity: text below (3b) and (4), printed p. 962.
- Final mutual impedance and kernel: (5a)–(5b), printed p. 963.
- Self substitution, total-matrix assembly, and Sunde/Pollaczek limits: text following (6c), printed p. 963.

## Limitations and discrepancies

- Conversion defect: the Markdown conversion contains the equation labels but drops nearly all symbols in (3)–(6); it is unsuitable for equation verification.
- The source does not state the square-root branch for ``\alpha_k`` in the inspected definition.
- The source calls ``\gamma_x\simeq j\omega\sqrt{\mu_1\epsilon_1}`` a high-frequency approximation but gives no formal bound; the numerical test range does not define its validity domain.
