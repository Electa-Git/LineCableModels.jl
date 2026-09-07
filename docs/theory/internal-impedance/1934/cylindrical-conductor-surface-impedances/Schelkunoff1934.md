# Schelkunoff cylindrical-conductor surface impedances

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Homogeneous solid cylinder of radius ``b``; or homogeneous annulus with inner radius ``a`` and outer radius ``b``. No insulation region enters these terms. |
| Calculated quantities | Solid-wire surface impedance; hollow-shell inner-surface, outer-surface, and transfer (mutual-surface) impedances |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical Bessel-function approximation within the reduced model. The physical reduction neglects conductor displacement current and drops ``\Gamma^2`` against ``\sigma^2`` before solving the radial Bessel equation. The displayed relation excludes the separate large-argument approximations (66)–(67). |
| Main source | S. A. Schelkunoff (1934) |
| Citation key(s) | `:Schelkunoff1934` |
| Evidence status | Verified visually against the original PDF page images |

**Description.** Per-unit-length longitudinal surface impedances of a homogeneous solid circular wire and a homogeneous hollow cylindrical shell. For the shell, the formulation relates the longitudinal electric field on each surface to the portions of conductor current associated with internal and external coaxial return paths and includes the transfer impedance between surfaces.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fields have the source dependence ``e^{-\Gamma z}``, but ``\Gamma^2`` is neglected relative to the metal intrinsic term ``\sigma^2`` in deriving (55). | Stated — exponential convention on printed p. 537; reduction from (54) to (55), printed pp. 548–549. |
| Air propagation constant ``γ_air`` | Not applicable. | The extracted surface impedances are solved inside the metal and contain no exterior-medium propagation constant. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Earth is not part of the formulation. |
| Earth permittivity and displacement current | Not applicable. | Earth is not part of the formulation. |
| Range of validity | The source requires ``\Gamma^2\ll\sigma^2`` for the reduced radial equation. Equations (65) and (75) themselves are not large- or small-argument Bessel approximations. | Stated — printed p. 549, paragraph preceding (55). |
| Earth permeability ``μ_earth`` | Not applicable. | Earth is not part of the formulation. |
| Arrangement | Coaxial return geometry. Solid-wire impedance is defined on the wire surface relative to total wire current. Shell currents may have internal and external return portions, producing inner, outer, and transfer terms. | Stated — printed pp. 551 and 553–554. |
| Earth structure | Not applicable. | Earth is not part of the formulation. |
| Conductor and insulation geometry | Homogeneous solid cylinder of radius ``b``; or homogeneous annulus with inner radius ``a`` and outer radius ``b``. No insulation region enters these terms. | Stated — section headings and text, printed pp. 551 and 552–554. |
| Constitutive and field assumptions | Circular symmetry, axial invariance apart from ``e^{-\Gamma z}``, homogeneous metal with scalar conductivity ``g`` and permeability ``\mu``. Conductor displacement current is neglected by setting ``\varepsilon=0``; proximity and eccentricity are absent. | Stated — printed p. 548 and geometry discussion on pp. 551–554. |
| Conventions | Implied time factor ``e^{j\omega t}``; longitudinal factor ``e^{-\Gamma z}``; ``j`` is the normalized imaginary unit, with the practical cgs-derived system and lengths in centimetres retained. Inner-surface enclosed current is ``-I_a`` and outer-surface enclosed current is ``I_b``. | Stated — printed pp. 533, 537, and 553. |

**Expression.** Exact modified-Bessel surface impedances within the source's reduced good-conductor field equation, equations (65) and (73)–(75).

```math
\begin{aligned}
\sigma^2&=j\omega\mu g=j2\pi f\mu g, \\
\eta&=\frac{\sigma}{g}=\frac{j\omega\mu}{\sigma},
\qquad \Re(\sigma)>0.
\end{aligned}
```

```math
Z_b=\frac{E_z(b)}{I}
=\frac{\eta I_0(\sigma b)}{2\pi b I_1(\sigma b)},
\qquad\text{(65)}
```

```math
D=I_1(\sigma b)K_1(\sigma a)-I_1(\sigma a)K_1(\sigma b),
\qquad\text{(73)}
```

```math
\begin{aligned}
E_z(a)&=Z_{aa}I_a+Z_{ab}I_b \\
E_z(b)&=Z_{ba}I_a+Z_{bb}I_b,
\end{aligned}\qquad\text{(74)}
```

```math
\begin{aligned}
Z_{aa}&=\frac{\eta}{2\pi aD}
\left[I_0(\sigma a)K_1(\sigma b)+K_0(\sigma a)I_1(\sigma b)\right], \\
Z_{bb}&=\frac{\eta}{2\pi bD}
\left[I_0(\sigma b)K_1(\sigma a)+K_0(\sigma b)I_1(\sigma a)\right], \\
Z_{ab}&=Z_{ba}=\frac{1}{2\pi g a bD}.
\end{aligned}\qquad\text{(75)}
```

``I_n`` and ``K_n`` are modified Bessel functions of order ``n`` of the first and second kinds. ``Z_b`` is the solid-wire surface impedance; ``Z_{aa}`` and ``Z_{bb}`` are the shell surface impedances for internal and external return, respectively; ``Z_{ab}=Z_{ba}`` is the transfer impedance between the two shell surfaces. The source reports these impedances in ``\Omega/\mathrm{cm}``.

**Approximation.** Not an analytical Bessel-function approximation within the reduced model. The physical reduction neglects conductor displacement current and drops ``\Gamma^2`` against ``\sigma^2`` before solving the radial Bessel equation. The separate large-argument approximations (66)–(67) are not silently substituted here.

**Limitations.** Circularly symmetric homogeneous cylindrical conductors and coaxial return paths only. The terms are surface/transfer impedances, not automatically complete terminal impedances of an arbitrary cable. The source's practical cgs-derived units must be converted explicitly in any later SI implementation. No earth or insulation contribution is present.

**Reference.** [Schelkunoff1934](@cite), equations (65), (73)–(75), printed pp. 551 and 554 (PDF pages 20 and 23), with field reduction and definitions on printed pp. 537 and 548–550.

**Transcription source.** Original publication. Every sign, radius factor, Bessel order, current orientation, equality ``Z_{ab}=Z_{ba}``, and the positive-real-part branch statement for ``\sigma`` was checked visually against the page images of the original PDF. Plain-text PDF extraction was used only for navigation.

## Source transcription

The displayed equations use the normalized ``j``, ``g``, ``\mu``, ``\sigma``, and ``\eta`` notation. The solid-wire field preceding (65) is

```math
E_z(\rho)=\frac{\eta I_0(\sigma\rho)}{2\pi b I_1(\sigma b)}I.
\qquad\text{(64)}
```

For the shell the source imposes

```math
\begin{aligned}
A I_1(\sigma a)+B K_1(\sigma a)&=-\frac{I_a}{2\pi a} \\
A I_1(\sigma b)+B K_1(\sigma b)&=\frac{I_b}{2\pi b},
\end{aligned}\qquad\text{(71)}
```

then derives (73)–(75). The source explicitly calls ``Z_{ab}`` the *transfer impedance* because it is not necessarily the total mutual impedance between two transmission lines.

## Numerical interpretation

The transfer term in the original printed (75) is positive:
``Z_{ab}=1/(2\pi gabD)``. The earlier negative sign in this
record was a transcription error, not a source convention.
The numerical surface term uses ``Z_{ms}=Z_{ab}``.
With the source's ``H_\varphi(a)=-I_a/(2\pi a)`` and
``H_\varphi(b)=I_b/(2\pi b)``, direct solution of (71)
reproduces all four entries of (74). The cable's loop-to-terminal
transformation is responsible for the separate transfer subtractions.

The evaluator uses SI resistivity, permeability, radii, and impedance
per metre. At zero frequency every annular surface coefficient tends
to ``R_{dc}=\rho/[\pi(b^2-a^2)]``. A solid conductor has only
the outer surface coefficient. Scaled Bessel ratios avoid large
intermediate exponentials; tests cover boundary excitations, dc,
negative-frequency conjugacy, thin walls, and the high-frequency
surface-layer limits.

## Notation map

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_b`` | unchanged | Solid-wire surface impedance per unit length | ``\Omega/\mathrm{cm}`` in source units |
| ``Z_{aa}`` | unchanged | Hollow-shell inner-surface impedance with internal return | ``\Omega/\mathrm{cm}`` |
| ``Z_{bb}`` | unchanged | Hollow-shell outer-surface impedance with external return | ``\Omega/\mathrm{cm}`` |
| ``Z_{ab}=Z_{ba}`` | unchanged | Transfer impedance between annulus surfaces | ``\Omega/\mathrm{cm}`` |
| ``E_z(a),E_z(b)`` | unchanged | Longitudinal electric-field intensity on inner/outer surfaces | volts/cm in source unit system |
| ``I,I_a,I_b`` | unchanged | Solid-wire total current; shell internal-return and external-return current portions | amperes; inner enclosed current is ``-I_a`` |
| ``a,b`` | unchanged | Inner and outer radii | centimetres in source unit system |
| ``\rho`` | unchanged | Radial coordinate | centimetres |
| ``g`` | unchanged | Metal conductivity | mhos/cm |
| ``\mu`` | unchanged | Metal permeability | henries/cm |
| ``\varepsilon`` | unchanged | Metal permittivity term | set to zero in conductor equations |
| ``\sigma`` | unchanged | Intrinsic metal propagation constant | ``\mathrm{cm}^{-1}``; root chosen with positive real part |
| ``\eta`` | unchanged | Intrinsic impedance of the metal | source-defined by (61) |
| ``\Gamma`` | unchanged | Longitudinal propagation constant | dependence ``e^{-\Gamma z}`` |
| ``I_n,K_n`` | unchanged | Modified Bessel functions of first/second kind | order ``n`` |
| ``j`` | unchanged | Imaginary unit | ``j^2=-1``; time factor ``e^{j\omega t}`` |
| ``\omega,f`` | unchanged | Angular frequency and frequency | ``\omega=2\pi f`` |

## Evidence and approximation sources

- Units and field/time conventions: footnote 1, printed p. 533; exponential convention, printed p. 537.
- Neglected conductor displacement current and pre-reduction equation containing ``\sigma^2-\Gamma^2``: (54), printed p. 548.
- ``\Gamma^2\ll\sigma^2`` reduction, radial Bessel solution, and branch for ``\sigma``: (55)–(57) and footnote 17, printed pp. 549–550.
- Intrinsic impedance ``\eta``: (60)–(61), printed p. 550.
- Solid-wire surface impedance: (64)–(65), printed p. 551.
- Hollow-shell current orientation and meaning of transfer impedance: printed p. 553.
- Hollow-shell boundary equations and surface/transfer impedances: (71)–(75), printed p. 554.

## Limitations and discrepancies

- Verification used the original PDF page images.
- The source uses centimetre-based practical units. This record preserves them and performs no hidden SI rescaling.
- The author distinguishes transfer impedance from total mutual line impedance; treating ``Z_{ab}`` as the latter without the rest of the field is unsupported.
