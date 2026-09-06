# Martins-Britto–Papadopoulos–Chrysochos generalized mixed impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel thin conductors, one in air and one in homogeneous soil. |
| Calculated quantities | Reciprocal generalized mixed mutual per-unit-length earth impedance |
| Earth structure | Homogeneous soil half-space below air. |
| Model and approximation | Integral representation within the paper's quasi-TEM model retaining longitudinal propagation, both permeabilities, and displacement current; numerical examples separately set ``k_x=0``. |
| Main source | A. G. Martins-Britto, T. A. Papadopoulos, and A. I. Chrysochos (2024), with earlier equivalent witnesses explicitly identified by the source |
| Citation key(s) | `:MartinsBritto2024` |
| Evidence status | Original publication page images checked |

**Description.** Reciprocal per-unit-length mutual earth impedance for two infinite parallel thin conductors placed on opposite sides of a homogeneous air/soil interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained as ``\gamma_x=jk_x`` for the excitation current; numerical examples set ``k_x=0``. | Stated — p. 984 and §III, p. 985. |
| Air propagation constant ``γ_air`` | ``\gamma_0=\sqrt{j\omega\mu_0(\sigma_0+j\omega\epsilon_0)}``; retained in ``a_0``. | Stated — p. 984 and definition below (4a). |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_1(\sigma_1+j\omega\epsilon_1)}``; retained in ``a_1``. | Stated — p. 984 and definition below (4a). |
| Earth permittivity and displacement current | Retained through ``\sigma_1+j\omega\epsilon_1``; neglecting it is a stated reduction to Pollaczek. | Stated — p. 984 and §II.B, p. 985. |
| Range of validity | Quasi-TEM derivation; Appendix A cites accuracy of the dominant mode up to 10 MHz for transmission-line problems. The paper's numerical validation spans 1 kHz–1 MHz; this is not a universal bound. | Stated — Appendix A p. 990 and §III.A p. 985. |
| Earth permeability ``μ_earth`` | Independent ``\mu_1`` retained; ``\mu_1=\mu_0`` is only a later simplification. | Stated — (5) and §II.B. |
| Arrangement | Mixed mutual air/soil interaction; reciprocity ``Z_{eij}^{01}=Z_{eij}^{10}``. | Stated — (5). |
| Earth structure | Homogeneous soil half-space below homogeneous air. | Stated — Fig. 1 and §II.A. |
| Conductor and insulation geometry | Infinite parallel thin conductors separated laterally by ``y_{ij}``; ``h_i>0`` and ``h_j<0`` are signed vertical coordinates. No insulation term enters the mixed external coefficient. | Stated — Fig. 1 and pp. 984–985. |
| Constitutive and field assumptions | Linear homogeneous media; quasi-TEM Hertzian-vector solution under Lorenz gauge. Internal impedance is assembled separately. | Stated — text before (5) and Appendix A. |
| Conventions | Vertical ``z`` points upward; interface ``z=0``; ``j=\sqrt{-1}``; result is per-unit-length. The source uses signed ``h_j<0``. | Stated — Fig. 1/§II.A and (5). |

**Expression.** The reciprocal mixed mutual impedance is

```math
Z_{eij}^{01}=Z_{eij}^{10}
=\frac{j\omega\mu_0\mu_1}{\pi}
\int_0^\infty
\frac{e^{-a_0h_i+a_1h_j}}
{a_0\mu_1+a_1\mu_0}
\cos(\lambda y_{ij})\,d\lambda,
\tag{5}
```

where

```math
\gamma_k=\sqrt{j\omega\mu_k(\sigma_k+j\omega\epsilon_k)},
\qquad
a_k=\sqrt{\lambda^2+\gamma_k^2+k_x^2},
\qquad
\gamma_x=jk_x,
```

for ``k=0,1``, with ``h_i>0`` and ``h_j<0``.

**Approximation.** Not an analytical approximation after imposing the source's quasi-TEM infinite-thin-conductor and homogeneous-medium model. Setting ``k_x=0``, ``\mu_1=\mu_0``, or omitting displacement current are later stated reductions and are not applied to (5) here.

**Limitations.** Homogeneous two-medium geometry only; no self term or layered recursion is supplied by (5). The source's Appendix-A quasi-TEM accuracy statement is cited from prior work and does not prove a universal 10 MHz boundary for every geometry/material set.

**Reference.** [MartinsBritto2024](@cite), equation (5), printed p. 985; definitions p. 984; derivation pp. 990–991.

**Transcription source.** Original IEEE page images. Superscripts, signed exponent, permeability-weighted denominator, root definitions, reciprocity, and longitudinal prescription were visually verified.

## Source transcription

The source derives (5) from the Hertzian-vector fields. Its Appendix B defines the mixed impedance by longitudinal integration in (B.4), applies the transform (B.14), and states that reciprocity supplies the reversed placement. Section II.B says (5) is identical to Dawalibi–Southey and Pawlik for homogeneous soil, reduces to Pollaczek when displacement current is neglected, and specializes to other cited nonmagnetic forms when ``\mu_1=\mu_0``.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{eij}^{01},Z_{eij}^{10}`` | unchanged | reciprocal mixed mutual earth impedances | ``\Omega/\mathrm m`` |
| superscripts ``01,10`` | unchanged | source/observation media, air ``0`` and soil ``1`` | source convention |
| ``h_i,h_j`` | unchanged | signed vertical coordinates | m; positive air, negative soil |
| ``y_{ij}`` | unchanged | horizontal separation | m |
| ``a_k`` | unchanged | longitudinally modified transverse root | ``\mathrm m^{-1}`` |

No notation was renamed.

## Evidence and approximation sources

The source presents a derivation but explicitly identifies earlier equivalent impedance formulas. The record therefore preserves Martins-Britto et al.'s inspected generalized witness without using it to erase Dawalibi–Southey, Pawlik, Pollaczek, Tsiamitros, Ametani, or Rallis as distinct historical/model/evaluator records.

## Limitations and discrepancies

- Numerical examples set ``k_x=0``; equation (5) itself retains ``k_x``.
- The paper's equivalence statements are author statements, not a corpus algebraic proof across all conventions.
- No mathematical correction was made to reconcile signed-coordinate and positive-depth versions in earlier sources.
