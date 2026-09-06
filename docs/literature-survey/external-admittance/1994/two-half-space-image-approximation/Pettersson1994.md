# Pettersson two-half-space shunt-admittance image approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Bare round thin wire represented by a line source; no coating. |
| Calculated quantities | Generalized p.u.l. self/mutual shunt admittance of a thin wire above, at, or below a lossy planar interface |
| Earth structure | Two homogeneous half-spaces with a flat boundary. |
| Model and approximation | Uses the same integrand replacement (9) as the series record, but ``b=n^2`` for ``Q``. On the interface the source uses (12)–(15), not the ``h>0`` distances above. |
| Main source | Pär Pettersson (1994 publication; 1993 conference manuscript) |
| Citation key(s) | `:Pettersson1994` |
| Evidence status | Original publication page images checked |

**Description.** Scalar-potential companion to Pettersson's series formulation, retaining conductive/displacement-current effects in both media and approximating the interface correction with a complex image distance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | ``\gamma_w`` is coupled to impedance through (4); the image evaluation uses the source's quasi-TEM substitutions. | Stated — (1), (4), text before (6)/(12). |
| Air propagation constant ``γ_air`` | ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``; air may be medium 1 or 2. | Stated — Basic Theory. |
| Earth propagation constant ``γ_earth`` | ``\gamma_m=\gamma_0n_m`` with conductive complex ``n_m``. | Stated — definitions following (1). |
| Earth permittivity and displacement current | Retained in ``n_m=[\varepsilon_{rm}+\sigma_m/(j\omega\varepsilon_0)]^{1/2}``. | Stated — below (1). |
| Range of validity | Thin-wire image approximation; no uniform error bound is stated. | Stated — pp. 1049, 1052. |
| Earth permeability ``μ_earth`` | Nonmagnetic media, ``\mu_0``. | Stated — p. 1049. |
| Arrangement | Parallel infinite thin wires; self/mutual coupling; either medium and interface placement. | Stated — abstract and Fig. 1. |
| Earth structure | Two homogeneous half-spaces with a flat boundary. | Stated — Basic Theory. |
| Conductor and insulation geometry | Bare round thin wire represented by a line source; no coating. | Stated — p. 1050. |
| Constitutive and field assumptions | Linear isotropic Helmholtz media; telegrapher reduction and selected classic mode. | Stated — pp. 1049–1050. |
| Conventions | ``e^{j\omega t}``; principal transverse square roots have positive real part. | Stated — p. 1049 and below (1). |

**Expression.** The exact shunt coupling and coupled mode equation are

```math
\begin{aligned}
Y&=j\omega\varepsilon_0n_1^2\,\frac{2\pi}{\Lambda+Q} \\
\gamma_w&=\gamma_1\left(\frac{\Lambda+P}{\Lambda+Q}\right)^{1/2}.
\end{aligned}\qquad\text{(3,4)}
```

For ``h>0`` the image form printed in (10) is

```math
\begin{aligned}
\Lambda&=\ln\frac{d''}{d'} \\
Q&\simeq\frac{2}{n^2+1}\ln\frac{d_Q}{d''} \\
d_Q&=\pm\sqrt{\left[y+h+\frac{n^2+1}{\beta}\right]^2+x^2}.
\end{aligned}\qquad\text{(10)}
```

Here ``\beta=\gamma_1(n^2-1)^{1/2}``; the sign of ``d_Q`` is selected so its imaginary part is negative for a wire in air and positive for a wire in ground.

**Approximation.** Uses the same integrand replacement (9) as the series record, but ``b=n^2`` for ``Q``. On the interface the source uses (12)–(15), not the ``h>0`` distances above.

**Limitations.** One planar interface, thin bare wires, nonmagnetic media, and the selected classic mode. ``Y`` is a generalized scalar self/mutual coupling, not by itself a multiconductor nodal-admittance assembly.

**Reference.** [Pettersson1994](@cite).  Pettersson, *IEEE Transactions on Power Delivery* 9(2), 1994, equations (1)–(17), printed pp. 1049–1052.

**Transcription source.** Original IEEE page images. The ``n_1^2`` multiplier, reciprocal ``\Lambda+Q``, ``2/(n^2+1)`` factor, image distance and sign prescription were visually verified.

## Source transcription

The scalar potential parent is ``V_m=I\gamma_w(j\omega\varepsilon_0n_m^2)^{-1}(\Lambda_m+Q_m)/(2\pi)`` in (1). This record preserves the explicit coupling to ``P`` through (4) rather than treating the admittance as electrostatic.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``A`` | ``\Lambda`` | direct/interface logarithmic term | dimensionless |
| ``Q`` | unchanged | scalar-potential interface correction | dimensionless |
| ``n=n_2/n_1`` | unchanged | complex refractive-index ratio | dimensionless |
| ``d_Q`` | unchanged | scalar complex-image distance | m |

Only the source's scalar ``A`` is renamed ``\Lambda``.

## Evidence and approximation sources

The paper says the ``Q`` approximation “may be new” and warns that the simplification available for ``\Lambda+P`` has no corresponding form for ``\Lambda+Q``.

## Limitations and discrepancies

- The printed radical/sign layout in (10) is easy to corrupt in text extraction; the original image controls it.
- The paper does not provide a matrix assembly or a universal error envelope, so neither is inferred here.

