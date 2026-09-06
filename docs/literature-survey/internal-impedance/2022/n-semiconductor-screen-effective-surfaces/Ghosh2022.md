# Ghosh–Das N-screen effective cable impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | One concentric cable; multiconductor assembly follows source matrix rules. Bonded concentric tubular conductor/semiconductor layers. |
| Calculated quantities | Effective inner/outer/mutual impedance recursion for cables with ``N`` semiconducting screens |
| Earth structure | None in this formula. |
| Model and approximation | Algebraic recursion is exact within the concentric tubular surface-impedance model. |
| Main source | Swarnankur Ghosh and Supriyo Das (2022) |
| Citation key(s) | `:Ghosh2022` |
| Evidence status | Original publication page images checked |

**Description.** Recursive combination of tubular conductor and semiconductor surface/transfer impedances for any number of concentric screens.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in radial tubular diffusion. | Derivation. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal cable model. |
| Earth propagation constant ``γ_earth`` | External earth matrix added separately. | Assembly. |
| Earth permittivity and displacement current | Not part of impedance recursion. | Scope. |
| Range of validity | Concentric underground cable with arbitrary count of semiconductor screens. | Title/model. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal model. |
| Arrangement | One concentric cable; multiconductor assembly follows source matrix rules. | (12). |
| Earth structure | None in this formula. | Scope. |
| Conductor and insulation geometry | Bonded concentric tubular conductor/semiconductor layers. | Figures. |
| Constitutive and field assumptions | Each homogeneous annulus uses Bessel surface impedances. | (1)–(5). |
| Conventions | Subscript ``in,out,m`` denotes the two surfaces and transfer coupling. | Definitions. |

**Expression.** For a semiconductor layer combined with an adjacent conductor layer,

```math
Z_{con_p,in}^{sem_p}=Z_{sem_p,in}-
\frac{Z_{sem_p,m}^2}{Z_{sem_p,out}+Z_{c_p,in}},
```

```math
\begin{aligned}
Z_{con_p,out}^{sem_p}&=Z_{c_p,out}-
\frac{Z_{c_p,m}^2}{Z_{sem_p,out}+Z_{c_p,in}} \\
Z_{con_p,m}^{sem_p}&=\frac{Z_{sem_p,m}Z_{c_p,m}}{Z_{sem_p,out}+Z_{c_p,in}}.
\end{aligned}
```

Equation (12) uses these effective surfaces to assemble diagonal and off-diagonal loop impedances for ``N`` screens.

**Approximation.** Algebraic recursion is exact within the concentric tubular surface-impedance model.

**Limitations.** Concentric bonded layers; no eccentric proximity fields or earth-return derivation.

**Reference.** [Ghosh2022](@cite), equations (1)–(12).

**Transcription source.** Original equation images; denominators and squared transfer terms verified.

## Source transcription

The internal recursion and the radial shunt network are separated into different taxonomy records.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{in},Z_{out}`` | unchanged | annulus surface impedances | ``Ω/m`` |
| ``Z_m`` | unchanged | through-wall transfer impedance | ``Ω/m`` |
| ``con_p,sem_p`` | unchanged | conductor/semiconductor layer labels | indexed |

## Evidence and approximation sources

This is a true ``N``-screen generalization beyond the earlier two-screen record.

## Limitations and discrepancies

The source indexing is retained; it should be validated carefully when coding arbitrary layer counts.
