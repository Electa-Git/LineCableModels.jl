# Weeks–Diao semiconducting-screen coaxial admittance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Single coaxial core-to-sheath radial path. Inner semiconducting screen, main insulation, outer semiconducting screen. |
| Calculated quantities | Core-to-sheath p.u.l. admittance of inner screen, lossy insulation, and outer screen in series |
| Earth structure | None. |
| Model and approximation | The series equation is the coaxial layer model; the ``y\simeq y_2`` reduction assumes both screen admittances greatly exceed the insulation admittance. |
| Main source | W. Weeks and Yi Diao (1984) |
| Citation key(s) | `:Weeks1984` |
| Evidence status | Original publication page images checked |

**Description.** Early explicit lossy-coaxial cable shunt model showing the two semiconducting screens and dielectric insulation as a series radial network.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Derived later from the assembled scalar line parameters. | Paper context. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal cable dielectric. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Insulation formula. |
| Earth permittivity and displacement current | Earth absent; insulation displacement and loss retained. | (10). |
| Range of validity | Long coaxial underground power cable with thin semiconducting layers. | Geometry. |
| Earth permeability ``μ_earth`` | Not applicable. | Shunt dielectric model. |
| Arrangement | Single coaxial core-to-sheath radial path. | §Characteristic Impedance. |
| Earth structure | None. | Scope. |
| Conductor and insulation geometry | Inner semiconducting screen, main insulation, outer semiconducting screen. | Text before (10). |
| Constitutive and field assumptions | Scalar coaxial quasistatic radial admittances; complex dielectric loss factor. | (10). |
| Conventions | Screen admittances ``y_1,y_3`` are in series with insulation ``y_2``. | Text/equation. |

**Expression.**

```math
\begin{aligned}
\frac1y&=\frac1{y_1}+\frac1{y_2}+\frac1{y_3} \\
\frac1y&\simeq\frac1{y_2}\quad(y_1,y_3\gg y_2),
\end{aligned}
```

```math
y_2=\frac{j\omega\varepsilon_2,2\pi}{\ln(a_2/a_1)}(1-jD_f)
=G_2+j\omega C_2.
\qquad\text{(10)}
```

**Approximation.** The series equation is the coaxial layer model; the ``y\simeq y_2`` reduction assumes both screen admittances greatly exceed the insulation admittance.

**Limitations.** Thin concentric layers and scalar cable mode; no eccentric/multiconductor dielectric matrix.

**Reference.** [Weeks1984](@cite), equations (9)–(10).

**Transcription source.** Original printed page images; the loss factor and radius ratio were visually checked.

## Source transcription

The paper's propagation and attenuation plots apply the explicit series shunt relation.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``y_1,y_3`` | unchanged | semiconductor-screen radial admittances | ``S/m`` |
| ``y_2`` | unchanged | lossy insulation admittance | ``S/m`` |
| ``D_f`` | unchanged | insulation dissipation factor | dimensionless |

## Evidence and approximation sources

This older formula was found through “semiconducting layer” and “lossy dielectric,” not earth-admittance titles.

## Limitations and discrepancies

The source's material frequency-independence assumption is explicitly cautioned by its authors.
