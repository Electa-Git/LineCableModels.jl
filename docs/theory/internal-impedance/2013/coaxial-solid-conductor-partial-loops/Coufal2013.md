# Coufal coaxial solid-conductor partial-loop impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | One inner and one outer coaxial return conductor. Two coaxial annuli separated by insulation. |
| Calculated quantities | Resistance and inductance of two coaxial tubular/solid conductors from radial partial loops |
| Earth structure | None. |
| Model and approximation | Current density is constant on each annular cell; convergence is obtained by increasing ``n``. |
| Main source | Oldřich Coufal (2013) |
| Citation key(s) | `:Coufal2013` |
| Evidence status | Original publication page images checked |

**Description.** Radial annular discretization of a coaxial conductor pair into coupled partial loops, followed by a phasor solve for loss and inductance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source explicitly assumes no ``z`` dependence/infinite field velocity. | p. 7. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal coaxial problem. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Internal problem. |
| Earth permittivity and displacement current | Not included. | Magnetic circuit. |
| Range of validity | Slowly varying current in long coaxial tubular conductors. | p. 7. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal problem. |
| Arrangement | One inner and one outer coaxial return conductor. | Figure 1. |
| Earth structure | None. | Internal problem. |
| Conductor and insulation geometry | Two coaxial annuli separated by insulation. | (40). |
| Constitutive and field assumptions | Piecewise-constant radial current density; linear magnetic coupling. | (42)–(53). |
| Conventions | Inner and outer partial-loop currents are opposite by (47). | (46)–(47). |

**Expression.** The ``n`` complex current-density unknowns satisfy

```math
\begin{aligned}
\left(\rho_i+\frac{\rho_o}{q}\right)\underline J_{ik}
+j\omega\sum_{\ell=1}^{n}\phi_{k\ell}\underline J_{i\ell}&=\underline U \\
k&=1,\ldots,n.
\end{aligned}\qquad\text{(53)}
```

The terminal current and impedance are

```math
\begin{aligned}
\underline I&=a\sum_{k=1}^{n}\underline J_{ik}
=-qa\sum_{k=1}^{n}\underline J_{ok} \\
\underline Z&=\frac{\underline U}{\underline I}=R_s+j\omega L.
\end{aligned}\qquad\text{(59–60)}
```

**Approximation.** Current density is constant on each annular cell; convergence is obtained by increasing ``n``.

**Limitations.** Coaxial geometry and quasistatic magnetic field only; not a general proximity model for separate axes.

**Reference.** [Coufal2013](@cite), equations (40)–(60).

**Transcription source.** Original open-access page images; equation (53) and terminal assembly verified.

## Source transcription

The paper calls reciprocal resistance “resistance-1” and terminal loss “resistance-2”; the symbols are retained without renaming those concepts.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``q`` | unchanged | outer/inner cross-section ratio | dimensionless |
| ``\phi_{k\ell}`` | unchanged | partial-loop flux coefficient | source units |
| ``a`` | unchanged | equal partial-cell area | ``m²`` |

## Limitations and discrepancies

The source flags a phase error in an earlier cited work; no correction is back-propagated beyond this paper's equations.
