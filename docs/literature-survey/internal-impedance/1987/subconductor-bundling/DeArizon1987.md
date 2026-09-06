# De Arizon–Dommel conductor-subdivision impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Arbitrary cables/conductors after subdivision and bundling. Cross-sections approximated by parallel conductive cells. |
| Calculated quantities | Cable impedance matrix from conductive subconductor subdivision, retaining skin and proximity redistribution |
| Earth structure | External primitive model dependent. |
| Model and approximation | Piecewise-uniform current in each cell; accuracy increases with cross-section subdivision. |
| Main source | Paloma De Arizon and Hermann W. Dommel (1987) |
| Citation key(s) | `:DeArizon1987` |
| Evidence status | Original publication page images checked |

**Description.** Primitive filament/subconductor circuit whose solved current redistribution yields cable series impedances with skin and proximity effects.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Neglected in the cross-section parameter extraction. | Derivation. |
| Air propagation constant ``γ_air`` | Not applicable to primitive conductor matrix. | Scope. |
| Earth propagation constant ``γ_earth`` | May be supplied in primitive external terms, not derived here. | Assembly. |
| Earth permittivity and displacement current | Not part of the subdivision formula. | Scope. |
| Range of validity | Parallel conductors represented by sufficiently fine subconductors. | Method. |
| Earth permeability ``μ_earth`` | Not part of internal redistribution. | Scope. |
| Arrangement | Arbitrary cables/conductors after subdivision and bundling. | (9)–(10). |
| Earth structure | External primitive model dependent. | Assembly. |
| Conductor and insulation geometry | Cross-sections approximated by parallel conductive cells. | Method. |
| Constitutive and field assumptions | Linear harmonic partial-inductance circuit. | (9). |
| Conventions | Reference-conductor reduction is performed after bundling. | (12). |

**Expression.**

```math
-\frac{d\mathbf V}{dx}=[\mathop{\mathrm{diag}}(R_i)+j\omega\mathbf L]\mathbf I,
\qquad\text{(9)}
```

followed by the paper's bundle reduction (10). With reference conductor ``k``,

```math
Z_{ij}^{*}=Z_{ij}+Z_{kk}-Z_{ik}-Z_{kj}.
\qquad\text{(12)}
```

**Approximation.** Piecewise-uniform current in each cell; accuracy increases with cross-section subdivision.

**Limitations.** Numerical primitive matrix; the record does not claim a closed Bessel solution.

**Reference.** [DeArizon1987](@cite), equations (9)–(12).

**Transcription source.** Original equation pages inspected.

## Source transcription

Skin and proximity arise from the solved unequal cell currents, not an empirical key or geometry label.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``R_i`` | unchanged | cell dc resistance per length | ``Ω/m`` |
| ``L`` | unchanged | primitive partial-inductance matrix | ``H/m`` |
| ``Z^*`` | unchanged | reference-reduced conductor matrix | ``Ω/m`` |

## Limitations and discrepancies

Mesh convergence remains geometry dependent.
