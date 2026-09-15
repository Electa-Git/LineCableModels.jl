# Lossy cable-layer admittivity (Ametani 2004 application)

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance; semicon admittance |
| Geometry | Semiconductor occupies annulus ``b'<r<c``; main insulation occupies ``c<r<r_0``; concentric cylindrical interfaces. |
| Calculated quantities | Shunt admittance of a cylindrical semiconducting layer; series radial combination with the main insulation admittance |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation within the scalar, concentric, radial dielectric model. The constitutive representation treats static resistivity ``\rho_2`` as a frequency-independent conduction term and adds it to displacement current through complex permittivity. |
| Main source | Standard frequency-domain constitutive relation; Ametani et al. (2004) application |
| Formula identifier | `:lossy` |
| Citation key(s) | `Ametani2004` |
| Evidence status | PDF page images checked |

**Description.** Generic lossy complex-admittivity relation for cable layers,
with Ametani et al. (2004) retained as an application reference for the
semiconducting-screen specialization and radial series assembly.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not retained in equations (14)–(15); the shunt relation is radial and per unit length. | Equation-implied — (14)–(15), printed p. 1525. |
| Air propagation constant ``γ_air`` | Not applicable. | No air region appears in the expression. |
| Earth propagation constant ``γ_earth`` | Not applicable. | The expression is earth-independent. |
| Earth permittivity and displacement current | Not applicable. | The expression concerns cable layers, not earth. |
| Range of validity | No explicit frequency bound is attached to (14)–(15). The paper's later conclusion says the admittance dominates semiconductor effects when layer thickness is small and resistivity high; this is a comparative observation, not a validity bound for (14). | Stated — conclusion item 5, printed p. 1530. |
| Earth permeability ``μ_earth`` | Not applicable. | The expression is earth-independent. |
| Arrangement | Not applicable. | Local coaxial shunt relation. |
| Earth structure | Not applicable. | The expression is earth-independent. |
| Conductor and insulation geometry | Semiconductor occupies annulus ``b'<r<c``; main insulation occupies ``c<r<r_0``; concentric cylindrical interfaces. | Stated — Fig. 1 and definitions after (15), printed pp. 1523, 1525. |
| Constitutive and field assumptions | Linear homogeneous isotropic layers represented by scalar permittivity; semiconductor conduction is incorporated as ``1/(j\omega\rho_2)`` in complex permittivity. | Equation-implied — (14), printed p. 1525. |
| Conventions | ``j=\sqrt{-1}``, ``\omega=2\pi f``; source uses ``e^{j\omega t}``; admittances are per unit length. | Stated — appendix opening and (14), printed pp. 1525, 1530. |

**Expression.** Source equations (14)–(15).

```math
\begin{aligned}
y_s&=\frac{j\omega\,2\pi\varepsilon_s}{\ln(c/b')} \\
\varepsilon_s&=\varepsilon_s'+\frac{1}{j\omega\rho_2},
\end{aligned}\qquad\text{(14)}
```

```math
\begin{aligned}
\frac{1}{Y}&=\frac{1}{y_s}+\frac{1}{y_i} \\
y_i&=\frac{j\omega\,2\pi\varepsilon_i}{\ln(r_0/c)}.
\end{aligned}\qquad\text{(15)}
```

``y_s`` is the semiconductor-layer shunt admittance, ``y_i`` the main-insulation shunt admittance, and ``Y`` their radial series combination between the core outer and sheath inner surfaces. ``r_0`` is the outer insulation radius, equal to the sheath inner radius.

**Approximation.** This is the standard scalar frequency-domain constitutive
relation with static material properties. It treats resistivity as a
frequency-independent conduction term and optionally adds polarization loss
through ``\tan\delta``. It is not an Ametani-specific empirical law.

**Limitations.** The formulation does not model frequency dependence of
material properties, anisotropy, interfacial polarization, or nonconcentric
geometry. It is a per-unit-length shunt admittance, not a longitudinal layer
impedance. The optional polarization-loss input must exclude conduction already
represented by resistivity.

**Application reference.** [Ametani2004](@cite). A. Ametani, Y. Miyamoto,
and N. Nagaoka, “Semiconducting Layer Impedance and its Effect on Cable
Wave-Propagation and Transient Characteristics,” *IEEE Transactions on Power
Delivery* 19(4), 1523–1531 (2004), DOI `10.1109/TPWRD.2003.822502`; equations
(14)–(15), printed p. 1525. The paper applies the standard relation to a
semiconducting screen and its radial combination with the main insulation.

**Source note.** The Ametani equations were checked against the original
publication PDF. They provide application geometry and context for this
registration; the constitutive relation is not claimed as newly developed by
that paper.

## Source transcription

The expressions above retain the source notation and order. The source text identifies the outer semiconducting layer on the insulation by reversing the material roles in Fig. 1; it does not print a separate equation for that placement.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``y_s`` | unchanged | Per-unit-length admittance of the semiconducting annulus | ``\mathrm{S/m}`` |
| ``y_i`` | unchanged | Per-unit-length admittance of the main insulation annulus | ``\mathrm{S/m}`` |
| ``Y`` | unchanged | Series radial combination of ``y_s`` and ``y_i`` | ``\mathrm{S/m}`` |
| ``\varepsilon_s,\varepsilon_s'`` | unchanged | Complex and real-part semiconductor permittivity | ``\mathrm{F/m}`` |
| ``\varepsilon_i`` | unchanged | Main-insulation permittivity | ``\mathrm{F/m}`` |
| ``\rho_2`` | unchanged | Semiconductor resistivity | ``\Omega\,\mathrm m`` |
| ``b',c,r_0`` | unchanged | Inner semiconductor, outer semiconductor, and outer insulation radii | m |
| ``j,\omega`` | unchanged | Imaginary unit and angular frequency | ``e^{j\omega t}`` |

## Evidence and approximation sources

- Geometry and medium assignment: Fig. 1, printed p. 1523.
- Complex permittivity and semiconductor admittance: (14), printed p. 1525.
- Radial series assembly and insulation definition: (15), printed p. 1525.
- Comparative statement that admittance dominates for thin, high-resistivity screens: conclusion item 5, printed p. 1530.

## Limitations and discrepancies

- The Markdown conversion emits only equation-number placeholders for (14)–(15); the source transcription comes from the PDF page image.
- The source writes ``\varepsilon_s=\varepsilon_s'+1/(j\omega\rho_2)``. This record preserves that convention and does not change the sign for a different phasor convention.

## Numerical interpretation

The executable layer coefficient is 2π(σ+jωε)/ln(r_out/r_in), with conduction included once. Insulation and semiconducting-screen selections share this evaluator. Physical annuli are combined in radial series before the nodal matrix transformation. The zero-frequency value is a constitutive limit; the public frequency-domain matrix interface requires a positive frequency.
