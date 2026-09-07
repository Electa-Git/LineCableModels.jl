# Pettersson two-half-space series-impedance image approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Filamentary line source representing a sufficiently thin round wire of radius ``a``; no coating layer. |
| Calculated quantities | Generalized p.u.l. self/mutual series impedance of a thin wire above, at, or below a lossy planar interface |
| Earth structure | Two homogeneous half-spaces separated by a plane. |
| Model and approximation | The exact spectral integrand in (6) is replaced by the asymptotically matched expression (9), which is then integrated using identity (7). The self case ``x=0,y=h+a`` gives ``P\simeq\ln[1+1/(\beta h)]`` in (11). Separate on-interface distances are printed in (14)–(15). |
| Main source | Pär Pettersson (1994 publication; 1993 conference manuscript) |
| Citation key(s) | `:Pettersson1994` |
| Evidence status | Original publication page images checked |

**Description.** Coupled electrodynamic potential formulation and a single image-type approximation for series impedance of a thin wire on either side of, or directly at, a lossy planar interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Unknown wire mode ``\gamma_w`` enters the exact transverse roots; the image approximation sets ``\gamma_w=\gamma_1`` for ``h>0`` and ``\gamma_w=\gamma_1[(n^2+1)/2]^{1/2}`` for ``h=0``. | Stated — (1), (4), text before (6) and (12). |
| Air propagation constant ``γ_air`` | A medium may be air; ``\gamma_m=\gamma_0n_m`` and ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``. | Stated — definitions following (1). |
| Earth propagation constant ``γ_earth`` | Included through ``n_m=[\varepsilon_{rm}+\sigma_m/(j\omega\varepsilon_0)]^{1/2}`` and ``\gamma_{tm}=(\gamma_m^2-\gamma_w^2)^{1/2}``. | Stated — definitions following (1). |
| Earth permittivity and displacement current | Retained through the complex refractive index of both media. | Stated — definitions following (1). |
| Range of validity | Thin-wire approximation; examples report close agreement through at least 1 MHz overhead, while buried image and quasi-TEM values underestimate part of the exact response. | Stated — p. 1052. |
| Earth permeability ``μ_earth`` | Both media are nonmagnetic, ``\mu=\mu_0``. | Stated — Basic Theory, p. 1049. |
| Arrangement | Self or mutual coupling between parallel infinite wires, each above/on/below the same planar interface. | Stated — abstract and Basic Theory. |
| Earth structure | Two homogeneous half-spaces separated by a plane. | Stated — Fig. 1. |
| Conductor and insulation geometry | Filamentary line source representing a sufficiently thin round wire of radius ``a``; no coating layer. | Stated — pp. 1049–1050. |
| Constitutive and field assumptions | Linear isotropic media; full two-medium Helmholtz potentials in the parent, followed by a quasi-TEM substitution and integrand approximation. | Stated — pp. 1050–1051. |
| Conventions | ``I=I_0e^{-\gamma_wz+j\omega t}``; principal square roots have positive real part. | Stated — p. 1049 and below (1). |

**Expression.** With medium-1 quantities understood, the exact generalized series impedance is

```math
\begin{aligned}
Z&=\frac{j\omega\mu_0}{2\pi}(\Lambda+P) \\
\gamma_w&=\gamma_1\left(\frac{\Lambda+P}{\Lambda+Q}\right)^{1/2}.
\end{aligned}\qquad\text{(3,4)}
```

For ``h>0`` the image approximation is

```math
\begin{aligned}
\Lambda&=\ln\frac{d''}{d'} \\
P&\simeq\ln\frac{d_P}{d''} \\
d_P&=\sqrt{(y+h+2/\beta)^2+x^2},
\end{aligned}\qquad\text{(10)}
```

so ``\Lambda+P\simeq\ln(d_P/d')``. The shared approximation has ``\beta=\gamma_1(n^2-1)^{1/2}`` and uses ``b=1`` for ``P``.

**Approximation.** The exact spectral integrand in (6) is replaced by the asymptotically matched expression (9), which is then integrated using identity (7). The self case ``x=0,y=h+a`` gives ``P\simeq\ln[1+1/(\beta h)]`` in (11). Separate on-interface distances are printed in (14)–(15).

**Limitations.** Thin wires, one flat interface, nonmagnetic media, and the selected classic structure mode only. Internal wire impedance must be added separately for a lossy conductor. The source does not claim uniform error bounds.

**Reference.** [Pettersson1994](@cite).  Pettersson, *IEEE Transactions on Power Delivery* 9(2), 1994, equations (1)–(17), printed pp. 1049–1052.

**Transcription source.** Original IEEE page images. The exact prefactor, coupled modal relation, logarithmic distance ratio, ``2/\beta`` image depth and the separate above/on/below prescriptions were visually checked.

## Source transcription

The exact series potential is ``A_{zm}=I\mu_0(\Lambda_m+P_m)/(2\pi)``. Its spectral ``P`` and the scalar-potential ``Q`` share transverse roots but have different interface denominators, so the impedance approximation cannot be reused as admittance without the printed ``Q`` factor.

## Numerical interpretation

The image evaluator uses (10)–(11) for two wires in the same half-space.
Medium 1 is the wire's medium, including for buried wires. Both conduction
and displacement current enter through
``\kappa_m=\sigma_m+j\omega\varepsilon_m`` and
``n^2=\kappa_2/\kappa_1``. The root
``\beta=\sqrt{\gamma_2^2-\gamma_1^2}`` is chosen with positive real part;
the product in the source definition does not override this branch rule.
The sign of the Q image follows (10), with conjugate values at negative
frequencies. Self entries use zero horizontal separation and the wire radius
as the direct distance in (11).

For two wires on the interface, the air-side reference in (14) gives

```math
\begin{aligned}
\Lambda&=0, \\
d_P&=\sqrt{\left[y+\frac{\sqrt{2}(1+j)}{\beta}\right]^2+x^2}, \\
d_Q&=\pm\sqrt{\left[y+\frac{j\sqrt{2}(n^2+1)}
{(n^2+j)\beta}\right]^2+x^2}.
\end{aligned}\qquad\text{(14)}
```

Mutual entries use ``y=0`` and the horizontal separation ``x``.
For self entries, (15) gives

```math
\begin{aligned}
P&\simeq\ln\frac{2\sqrt{j}}{\beta a}, \\
Q&\simeq\frac{2}{n^2+1}
\ln\frac{j\sqrt{2}(n^2+1)}{(n^2+j)\beta a}.
\end{aligned}\qquad\text{(15)}
```

This last reduction requires ``|\beta a|\ll1``; values with
``|\beta a|\geq1`` are rejected. Mixed air/ground and
interface/off-interface mutual pairs are not assigned these image distances:
they require a common longitudinal-mode prescription, which the two
separate image substitutions do not supply. The coupled mode equation (4)
is not solved by this evaluator.

Series entries are ``j\omega\mu_0(\Lambda+P)/(2\pi)``.
Potential entries are ``j\omega(\Lambda+Q)/(2\pi\kappa_1)``;
the engine forms ``\mathbf Y=j\omega\mathbf P^{-1}`` after assembling
all potential contributions, rather than inverting individual mutual entries.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``A`` | ``\Lambda`` | source direct/interface logarithmic term (renamed only to avoid collision with vector potential) | dimensionless |
| ``P`` | unchanged | magnetic-vector-potential interface correction | dimensionless |
| ``\gamma_w`` | unchanged | unknown axial wire-mode constant | ``\mathrm m^{-1}`` |
| ``d',d'',d_P`` | unchanged | physical, mirror and complex-image distances | m |

The source's scalar ``A`` is displayed as ``\Lambda``; all other symbols are retained.

## Evidence and approximation sources

The source distinguishes the parent integral, quasi-TEM numerical evaluation, and proposed image method. Equation (9) is the approximation-generating step.

## Limitations and discrepancies

- Optical text extraction corrupts several radicals and primes; the page image, not the conversion, controls this transcription.
- The publication says the ``P`` part of (11) was known from Kostenko/Sunde, while the single above/on/below scheme and corresponding ``Q`` approximation are the claimed contribution; priority is not transferred to the whole parent equation.
