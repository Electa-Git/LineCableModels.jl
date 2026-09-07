# Morgan steel-cored aluminium conductor characteristics

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | One composite conductor; external proximity field may be included. Steel core and helically stranded aluminium layers. |
| Calculated quantities | AC resistance and internal inductance of steel-cored aluminium stranded conductors |
| Earth structure | None. |
| Model and approximation | Composite-strand field and material effects are represented by the source's layer/current approximations; this is not the exact field of every helical strand. |
| Main source | V. T. Morgan (1965) |
| Citation key(s) | `:Morgan1965` |
| Evidence status | Original publication page images checked |

**Description.** Historical ACSR model combining steel-core current, magnetic/hysteretic terms, skin effect, proximity effect, and stranding geometry.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not a wave-propagation solution; longitudinal conductor parameter. | Scope. |
| Air propagation constant ``γ_air`` | Not applicable. | Internal model. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Internal model. |
| Earth permittivity and displacement current | Not applicable. | Internal model. |
| Range of validity | Steel-cored aluminium stranded overhead conductor under source frequency/material restrictions. | Title/derivation. |
| Earth permeability ``μ_earth`` | Not applicable. | Internal model. |
| Arrangement | One composite conductor; external proximity field may be included. | §§4–5. |
| Earth structure | None. | Internal model. |
| Conductor and insulation geometry | Steel core and helically stranded aluminium layers. | Geometry. |
| Constitutive and field assumptions | Effective steel permeability/loss and strand-current distributions as developed by source. | (9)–(21). |
| Conventions | ``p`` is the proportion of total current in the steel core. | Definitions before (27). |

**Expression.** The resulting internal inductance is assembled as

```math
L_i=2\ln\frac{d}{2d_{11}}+\Delta L+\frac{\mu p^2}{2}+\delta L,
\qquad\text{(27)}
```

where the paper defines the strand/self geometric term, magnetic steel-core increment ``\Delta L``, current fraction ``p``, and skin/proximity correction ``\delta L`` through (6)–(26). The corresponding resistance is built from the dc sharing relations (6)–(8), magnetic/hysteretic losses (9)–(17), and skin/proximity terms (18)–(21).

**Approximation.** Composite-strand field and material effects are represented by the source's layer/current approximations; this is not the exact field of every helical strand.

**Limitations.** Historical unit convention and empirical magnetic-core properties require care before SI implementation.

**Reference.** [Morgan1965](@cite), equations (6)–(27).

**Transcription source.** Original page images; ``p^2`` in (27) was visually distinguished from a geometric symbol.

## Source transcription

The source's component equations remain dependencies because (27) is an assembly rather than an independent closed evaluator.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``p`` | unchanged | steel-core current fraction | dimensionless |
| ``d,d_{11}`` | unchanged | conductor/strand geometric distances | source units |
| ``\Delta L,\delta L`` | unchanged | magnetic and skin/proximity increments | source units |

## Evidence and approximation sources

This formula was found by the ACSR/skin/proximity branch of the full-text sweep.

## Limitations and discrepancies

No modern unit conversion or magnetic-loss refit is inferred.

## Coaxial assimilation scope — Deferred

The complete model requires operating current and the steel's
magnetizing-force-dependent permeability and loss angle. In source
centimetre units, (11) and (17) give

```math
\begin{aligned}
H&=\frac{4\pi I}{10q}\sum_m N_m,\\
\Delta R&=\frac{8\pi^2af}{q^2}
\left(\sum_m N_m\right)^2\mu(H)\tan\delta(H)\,10^{-9}.
\end{aligned}
```

The signed turns follow the alternating lays. Printed p. 332 instructs
the reader to calculate ``H`` from current and then obtain
``\mu`` and ``\tan\delta`` from the measured curves in Fig. 7.
The current-independent, real-permeability coaxial input does not
provide that constitutive calculation. Selecting arbitrary fixed curve
values would omit the model's primary current-dependent loss mechanism.
The record is therefore deferred, not rejected because it contains
helical strands.

The separate dc parallel-sharing construction is already represented by
the geometric conductor reduction. Equations (18)–(21) reproduce
Arnold's tubular skin/proximity approximation; their pair-spacing term
is not inserted into an isolated-conductor surface coefficient.
The source's inductance example on p. 334 supplies the
``10^{-9}\ \mathrm{H/cm}`` scale implicit in the bracket of (27).
