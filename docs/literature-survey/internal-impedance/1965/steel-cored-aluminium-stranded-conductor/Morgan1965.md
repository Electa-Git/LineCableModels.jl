# Morgan steel-cored aluminium conductor characteristics

## Identity and source

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
\tag{27}
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
