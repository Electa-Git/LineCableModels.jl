# Mouhaidali–Chadebec multilayer HVDC cable FEM admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | ``k+1`` conductor system; voltage excitations per conductor. Finite cable and semiconductor regions. |
| Calculated quantities | Cable shunt-admittance matrix in multilayer earth and dielectric regions from FEM energy/power extraction |
| Earth structure | Arbitrary meshed earth/seawater/seabed regions. |
| Model and approximation | Finite-element mesh and outer-domain truncation; material loss mechanisms enter only through their combined complex constitutive parameters. |
| Main source | A. Mouhaidali, O. Chadebec, S. Silvant, D. Tromeur-Dervout, and J.-M. Guichon (2018) |
| Citation key(s) | `:Mouhaidali2018` |
| Evidence status | Original-publication page image/text checked |

**Description.** Electrodynamic/electric-potential FEM used to extract shunt conductance and capacitance in realistic multilayer cable surroundings.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally homogeneous cross-section. | §III-A. |
| Air propagation constant ``γ_air`` | Quasistatic cross-sectional solve. | Formulation. |
| Earth propagation constant ``γ_earth`` | Constitutive conductive/dielectric regions solved directly. | (10). |
| Earth permittivity and displacement current | Retained, including complex permittivity/loss. | (10)–(11). |
| Range of validity | Parallel cables and meshed material layers. | Abstract. |
| Earth permeability ``μ_earth`` | Not controlling electric solve. | §II-B. |
| Arrangement | ``k+1`` conductor system; voltage excitations per conductor. | §II-D. |
| Earth structure | Arbitrary meshed earth/seawater/seabed regions. | Case study. |
| Conductor and insulation geometry | Finite cable and semiconductor regions. | Model. |
| Constitutive and field assumptions | Linear harmonic scalar-potential FEM. | (10). |
| Conventions | ``Y=G+jωC`` from complex power under unit voltage. | §III-C. |

**Expression.** The electric solve is the source's equation (10), with loss tangent

```math
\tan\delta=\frac{\omega\epsilon''+\sigma}{\omega\epsilon'}.
\qquad\text{(11)}
```

For an applied voltage ``V``, the source extracts

```math
C=\frac1{2\pi f}\frac{\Im S(\omega)}{V^2},
\qquad G=\frac{\Re S(\omega)}{V^2},
\qquad Y=G+j\omega C.
```

**Approximation.** Finite-element mesh and outer-domain truncation; material loss mechanisms enter only through their combined complex constitutive parameters.

**Limitations.** The printed FEM method yields a total shunt matrix, not a closed external-earth coefficient.

**Reference.** [Mouhaidali2018](@cite), equations (10)–(11) and §III-C.

**Transcription source.** Original formula and extraction pages inspected.

## Source transcription

The result is classified by shunt output, separately from the paper's series solve.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``S(ω)`` | unchanged | complex power in dielectric region | VA per length |
| ``G,C`` | unchanged | shunt conductance/capacitance | ``S/m``, ``F/m`` |
| ``ε',ε''`` | unchanged | complex-permittivity components | ``F/m`` |

## Evidence and approximation sources

The formulation includes multilayer earth and sea installations.

## Limitations and discrepancies

The source notes that conduction and dielectric relaxation losses cannot be separated from the combined loss tangent.
