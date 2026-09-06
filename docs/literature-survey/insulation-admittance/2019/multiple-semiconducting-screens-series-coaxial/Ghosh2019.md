# Ghosh–Ghosh–Das multiple-semiconducting-screen coaxial admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Geometry | Two concentric semiconducting annuli separated by the main insulation, between core and sheath. |
| Calculated quantities | Effective core-to-sheath shunt admittance of two semiconducting screens and the intervening cable insulation |
| Earth structure | None in the shunt expression. |
| Model and approximation | Exact series combination within the coaxial quasistatic layer model. Each annulus is reduced to its scalar radial p.u.l. admittance. |
| Main source | Swarnankur Ghosh, Mousam Ghosh, and Supriyo Das (2019 accessible multiple-screen assembly) |
| Citation key(s) | `:Ghosh2019` |
| Evidence status | Original publication page images checked |

**Description.** Coaxial radial shunt path that treats the inner semiconducting screen, main insulation and outer semiconducting screen as three series admittances with complex permittivity.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Transverse coaxial shunt field; longitudinal variation neglected in the p.u.l. reduction. | Stated — §II-E. |
| Air propagation constant ``γ_air`` | Not applicable. | Scope. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Scope. |
| Earth permittivity and displacement current | Not an earth model; displacement and conduction in semiconductors are combined by complex permittivity. | Stated — (36). |
| Range of validity | Harmonic, concentric, quasistatic coaxial shunt path; no explicit universal frequency limit. | Stated — §II-E. |
| Earth permeability ``μ_earth`` | Not applicable. | Scope. |
| Arrangement | Single-core coaxial cable. | Stated — Fig. 1. |
| Earth structure | None in the shunt expression. | Scope. |
| Conductor and insulation geometry | Two concentric semiconducting annuli separated by the main insulation, between core and sheath. | Stated — Fig. 1, §II-E. |
| Constitutive and field assumptions | Linear homogeneous isotropic dielectric/semiconductor layers; purely radial coaxial electric field. | Stated — (35)–(37). |
| Conventions | ``\omega=2\pi f``; complex relative permittivity uses the printed ``+1/(j\omega\rho)`` for semiconductor and ``\epsilon_r'-j\epsilon_r''`` for insulation. | Stated — (36)–(37). |

**Expression.**

```math
Y_{e12}=\left(\frac1{y_{sem1}}+\frac1{y_{ins12}}+\frac1{y_{sem2}}\right)^{-1}
=G_{e12}+j\omega C_{e12},
\tag{35}
```

```math
y_{semi}=\frac{j\omega\epsilon_0\epsilon_{semi}}{\ln(r_{out}/r_{in})},qquad
\epsilon_{semi}=\epsilon_{r,semi}+\frac1{j\omega\rho_{semi}},
\tag{36}
```

```math
y_{ins12}=\frac{j\omega\epsilon_0\epsilon_{r,ins12}}{\ln(r_{out}/r_{in})},qquad
\epsilon_{r,ins12}=\epsilon_r'-j\epsilon_r''.
\tag{37}
```

**Approximation.** Exact series combination within the coaxial quasistatic layer model. Each annulus is reduced to its scalar radial p.u.l. admittance.

**Limitations.** Concentric layers, scalar isotropic material properties, no eccentric/proximity electric field, no longitudinal dielectric modes, and no external ground-potential coefficient.

**Reference.** [Ghosh2019](@cite).  Ghosh, Ghosh, and Das, *IEEE Access* 7, 2019, DOI `10.1109/ACCESS.2019.2955026`, equations (35)–(37), printed p. 169376.

**Transcription source.** Original IEEE Access page image. The reciprocal-series structure, logarithmic radius ratio, plus sign before the conductive complex-permittivity term and insulation loss sign were visually verified.

## Source transcription

The source analyzes only the core-to-sheath propagation path, so ``Y_{e12}`` is a scalar p.u.l. branch rather than a general nodal admittance matrix.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``y_{sem1},y_{sem2}`` | unchanged | individual semiconducting-screen radial admittances | ``\mathrm{S/m}`` |
| ``y_{ins12}`` | unchanged | main-insulation radial admittance | ``\mathrm{S/m}`` |
| ``Y_{e12}`` | unchanged | combined core-to-sheath branch admittance | ``\mathrm{S/m}`` |
| ``\rho_{semi}`` | unchanged | semiconducting-screen resistivity | ``\Omega\,m`` |

No notation was renamed.

## Evidence and approximation sources

The individual complex-permittivity layer rule is conventional; the accessible contribution is the explicit two-screen series assembly used with the paper's new impedance model.

## Limitations and discrepancies

- The source's prose calls ``\rho_{semi}`` resistivity while its complex-permittivity convention embeds it through ``1/(j\omega\rho)``; this is retained exactly.
- Citation priority for the single-screen admittance rests partly on source reference [24] and is not reassigned here.

