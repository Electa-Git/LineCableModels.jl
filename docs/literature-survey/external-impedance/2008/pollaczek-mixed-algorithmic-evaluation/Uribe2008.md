# Uribe algorithmic evaluation of the mixed Pollaczek integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel filamentary axes at overhead height, burial depth, and horizontal separation. |
| Calculated quantities | Positive-domain algorithmic evaluation of Pollaczek's mixed mutual ground impedance, with tail and zero-crossing prescriptions |
| Earth structure | Homogeneous conductive soil half-space below air. |
| Model and approximation | Exact integral transformation of the inherited Pollaczek kernel; finite tail truncation and empirical ``\lambda_e=12`` are numerical approximations. |
| Main source | F. A. Uribe (2008) for the evaluator; F. Pollaczek for the parent kernel |
| Citation key(s) | `:Uribe2008` |
| Evidence status | Original publication page images checked |

**Description.** Algorithmic representation and finite-range quadrature prescription for Pollaczek's per-unit-length mutual ground impedance between one overhead and one buried conductor over homogeneous soil.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed to the quasi-TEM zero-longitudinal-propagation Pollaczek parent. | Stated — text introducing (1a), p. 198. |
| Air propagation constant ``γ_air`` | Neglected in the parent; the air side contributes ``|\beta|``. | Equation-implied — (1a), p. 198. |
| Earth propagation constant ``γ_earth`` | ``1/p=\sqrt{j\omega\mu_0\sigma}``. | Stated — (1b), p. 198. |
| Earth permittivity and displacement current | Neglected; ``p`` contains conductivity only. | Equation-implied — (1b). |
| Range of validity | The algorithm was tested over Table I: ``0.1\le h_1,h_2\le10^2\,\mathrm m``, ``10^{-2}\le x\le600\,\mathrm m``, ``2\pi\le\omega\le2\pi10^6\,\mathrm{rad/s}``, and ``10^{-4}\le\sigma\le1\,\mathrm{S/m}``. These are tested ranges, not a universal theorem. | Stated — Table I, p. 200. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — (1a)–(1b). |
| Arrangement | Mutual interaction between one overhead and one buried conductor. | Stated — Fig. 1 and (1a). |
| Earth structure | Homogeneous imperfectly conducting soil half-space below air. | Stated — Fig. 1 caption and text before (1a). |
| Conductor and insulation geometry | Parallel filamentary axes at overhead height ``h_1``, burial depth ``h_2``, and horizontal separation ``x``; no insulation contribution is included. | Stated — nomenclature, Fig. 1, and (1a). |
| Constitutive and field assumptions | Linear homogeneous soil; quasi-TEM Pollaczek physical kernel. The contribution is an evaluator and convergence partition, not a new interface kernel. | Stated — §§I–III. |
| Conventions | ``j=\sqrt{-1}``; ``h_1,h_2,x`` are positive; ``\beta=u/|p|``; output is per-unit-length mutual ground impedance. | Stated — nomenclature and (1a)–(2e). |

**Expression.** The exact parent and Uribe's dimensionless positive-domain evaluator are

```math
Z_G(\omega)=\frac{j\omega\mu_0}{2\pi}
\int_{-\infty}^{+\infty}
\frac{e^{-h_1|\beta|}e^{-h_2\sqrt{\beta^2+1/p^2}}}
{|\beta|+\sqrt{\beta^2+1/p^2}}
e^{j\beta x}\,d\beta,
\qquad
p=\frac{1}{\sqrt{j\omega\mu_0\sigma}},
\qquad\text{(1a--1b)}
```

```math
Z_G=\frac{\omega\mu_0}{\pi}J(\xi,\eta),
```

```math
J(\xi,\eta)=\int_0^{+\infty}
[F(u)-u+jG(u)]e^{-\xi[\zeta u+F(u)]}
e^{-j\xi G(u)}\cos(\xi\eta u)\,du,
\qquad\text{(2e)}
```

where

```math
F(u)=\frac{\sqrt{u^2+\sqrt{u^4+1}}}{\sqrt2},\qquad
G(u)=\frac{\sqrt{-u^2+\sqrt{u^4+1}}}{\sqrt2},
```

```math
\xi=\frac{h_2}{|p|},\qquad \eta=\frac{x}{h_2},\qquad
\zeta=\frac{h_1}{h_2}.
```

The source truncates at

```math
u_{\max}=\frac{\lambda_e}{\xi(\zeta+1)},\qquad
\lambda_e=-\log[\epsilon_r\xi(\zeta+1)],
\qquad\text{(4c--4d)}
```

and reports the empirical choice ``\lambda_e=12`` for its applications.

**Approximation.** Equations (1a)–(2e) are algebraic/integral transformations of the inherited Pollaczek kernel. Replacing the infinite range by ``[0,u_{\max}]`` is a numerical tail approximation controlled by (4b)–(4d); ``\lambda_e=12`` is empirical. The regular and irregular zero-crossing formulas in (4h), (4k), and (4l) partition the quadrature intervals rather than alter the physical kernel.

**Limitations.** Homogeneous soil, fixed ``\mu_0``, no displacement current, parallel infinite axes, and mutual mixed geometry only. Table I and the reported RMS errors describe the author's test campaign, not a proof outside it.

**Reference.** [Uribe2008](@cite), equations (1a)–(5b), printed pp. 198–200.

**Transcription source.** Original IEEE page images. Absolute values, both exponentials, the complex depth, the nested roots defining ``F,G``, dimensionless variables, and truncation logarithm were visually verified.

## Source transcription

The source also normalizes ``\underline Z_G=Z_G(\pi/(\omega\mu_0))=J(\xi,\eta)`` in (5a)–(5b). Its regular zero crossings are ``u_k=\pi(2k-1)/(2\xi\eta)``; irregular crossings are obtained by inverting ``G(u)``. These numerical partition identities are dependencies of the algorithm but do not define a second physical impedance.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_G`` | unchanged | mixed mutual ground impedance | ``\Omega/\mathrm m`` |
| ``p`` | unchanged | complex depth | m |
| ``h_1,h_2,x`` | unchanged | overhead height, burial depth, horizontal separation | m |
| ``\xi,\eta,\zeta`` | unchanged | normalized depth, separation, and height ratios | dimensionless |
| ``F,G`` | unchanged | real and imaginary parts of ``\sqrt{u^2+j}`` | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Pollaczek supplies the physical integral. Uribe extends his earlier buried-conductor integration strategy to the mixed case by positive-domain normalization, tail selection, and regular/irregular zero-crossing partition. The evaluator is therefore a distinct published numerical formulation without being a new earth-return kernel.

## Limitations and discrepancies

- The source calls the solution convergent and reports numerical errors, but supplies no global error theorem for arbitrary physical parameters.
- ``\lambda_e=12`` is explicitly empirical.
- The closed-form comparison formulas printed later in the paper are retained in a separate record and are not substituted into this evaluator.
