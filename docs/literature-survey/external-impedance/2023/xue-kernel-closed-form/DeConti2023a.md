# De Conti–Duarte–Alipio closed-form underground impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Insulated circular cables at depths ``h_m,h_n``; mutual horizontal separation ``r``; self sets ``h_m=h_n`` and ``r=r_o``, the external radius including insulation. |
| Calculated quantities | Self and mutual closed-form approximation of Xue's homogeneous-earth underground-cable ground-return impedance |
| Earth structure | Homogeneous earth below homogeneous air. |
| Model and approximation | The source begins with ``Z_{g(m,n)}=j\omega\mu_0[\Lambda+\Theta_1]/(2\pi)`` and approximates the square-root ratio in ``d\Theta_1/dH`` by a constant plus an exponentially decaying term. It then replaces ``e^{-H\sqrt{\lambda^2+\gamma_1^2}}`` by ``e^{-H\gamma_1}`` only in that residual and integrates using the Bessel identity (12), yielding (13). No series order or formal remainder bound is provided. |
| Main source | A. De Conti, N. Duarte, and R. Alipio (2023) |
| Citation key(s) | `:DeConti2023a` |
| Evidence status | Original publication page images checked |

**Description.** Closed-form approximation of the self or mutual per-unit-length ground-return impedance between insulated cables buried in a homogeneous earth, retaining nonzero air and earth bulk propagation constants.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed zero in the compact Xue equations (8)–(11); the text explicitly contrasts this with another nonzero prescription. | Stated — discussion below (11), p. 2892. |
| Air propagation constant ``γ_air`` | ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``, with air conductivity zero. | Stated — definitions below (6), p. 2892. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_1(\sigma_1+j\omega\varepsilon_1)}``. | Stated — p. 2892. |
| Earth permittivity and displacement current | Retained through ``\varepsilon_1=\varepsilon_{r1}\varepsilon_0``. | Stated — p. 2892. |
| Range of validity | Numerical study uses 100 Hz–10 MHz, ``100\le\rho\le10000\ \Omega\,\mathrm m``, ``\varepsilon_{r1}=10``, depths/separations 0.5–2 m; errors depend inversely on separation, frequency and resistivity. These are tested ranges, not universal bounds. | Stated — §§IV and VII, pp. 2893–2900. |
| Earth permeability ``μ_earth`` | ``\mu_1=\mu_0``. | Stated — p. 2892. |
| Arrangement | Underground; self and mutual cable terms. | Stated — Fig. 1 and matrix assembly, p. 2892. |
| Earth structure | Homogeneous earth below homogeneous air. | Stated — Fig. 1. |
| Conductor and insulation geometry | Insulated circular cables at depths ``h_m,h_n``; mutual horizontal separation ``r``; self sets ``h_m=h_n`` and ``r=r_o``, the external radius including insulation. | Stated — matrix-assembly paragraph, p. 2892. |
| Constitutive and field assumptions | Linear homogeneous media, quasi-TEM Xue parent, nonmagnetic ground, no internal/insulation fields in this contribution. | Stated — §II and definitions, p. 2892. |
| Conventions | ``j=\sqrt{-1}``; positive burial depths; per-unit-length impedance. | Stated — definitions, p. 2892. |

**Expression.** Proposed closed-form impedance, equation (13), printed p. 2893.

```math
Z^{\mathrm{app}}_{g(m,n)}=\frac{j\omega\mu_0}{2\pi}\left[
K_0(\gamma_1d)+
\frac{\gamma_1-\gamma_0}{\gamma_0+\gamma_1}
e^{-(h_m+h_n)\gamma_1}
\left(\frac{2}{4+\gamma_1^2r^2}\right)
\right],
\qquad\text{(13)}
```

```math
d=\sqrt{(h_m-h_n)^2+r^2},\quad
D=\sqrt{(h_m+h_n)^2+r^2},\quad
\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0},\quad
\gamma_1=\sqrt{j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)}.
```

**Approximation.** The source begins with ``Z_{g(m,n)}=j\omega\mu_0[\Lambda+\Theta_1]/(2\pi)`` and approximates the square-root ratio in ``d\Theta_1/dH`` by a constant plus an exponentially decaying term. It then replaces ``e^{-H\sqrt{\lambda^2+\gamma_1^2}}`` by ``e^{-H\gamma_1}`` only in that residual and integrates using the Bessel identity (12), yielding (13). No series order or remainder bound is provided.

**Limitations.** One homogeneous earth half-space, nonmagnetic media, parallel infinite cables, and quasi-TEM. The error grows for greater separation and frequency and lower resistivity. Square-root/Bessel branches are unstated. The approximation is not the full Xue integral even though it is validated against it.

**Reference.** [DeConti2023a](@cite), equations (1), (3)–(13), printed pp. 2892–2893 (PDF pages 2–3).

**Transcription source.** Original IEEE page images. Equation (13)'s reflection ratio, decay exponent, rational denominator, geometry definitions, self substitution and full bulk constants were visually verified.

## Source transcription

The compact parent is

```math
Z_{g(m,n)}=\frac{j\omega\mu_0}{2\pi}[\Lambda+\Theta_1],
\qquad
\Lambda=K_0(\gamma_1d)-K_0(\gamma_1D),
\qquad\text{(8,3)}
```

```math
\Theta_1=2\int_0^\infty
\frac{e^{-(h_m+h_n)u_1}}{u_1+u_0}\cos(r\lambda)\,d\lambda,
\quad u_i=\sqrt{\lambda^2+\gamma_i^2}.
\qquad\text{(10)}
```

This retains the air term that disappears in Sunde's limiting equation.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{g(m,n)}`` | unchanged | self/mutual ground-return impedance | ``\Omega/\mathrm m`` |
| ``h_m,h_n,r`` | unchanged | burial depths and horizontal separation/self radius | m |
| ``d,D`` | unchanged | direct and image distances | m |
| ``\gamma_0,\gamma_1`` | unchanged | air and earth bulk constants | ``\mathrm m^{-1}`` |
| ``u_0,u_1`` | unchanged | spectral roots | ``\mathrm m^{-1}`` |
| ``\Lambda,\Theta_1`` | unchanged | Bessel direct/image and interface terms | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Equations (1)–(11) are attributed by the paper to Xue and algebraically compacted; only the asymptotic steps in §III-A and final (13) are this paper's closed-form contribution. The source reports reduction to Saad–Gaba–Giroux when ``\gamma_0=0`` and ``\gamma_1=\sqrt{j\omega\mu_1\sigma_1}``.

## Limitations and discrepancies

- The formulation is an approximation of Xue's impedance, not a new full-wave parent.
