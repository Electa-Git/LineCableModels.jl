# Data entry validation

## Contents

```@contents
Pages = ["validation.md"]
Depth = 3
```

[`validate`](@ref) checks materialized inputs before they enter another
construction or numerical operation. It returns its argument unchanged when
the value is usable and throws a native Julia exception when it is not:

```julia
validated = validate(value)
@assert validated === value
```

Input validation never converts, fills defaults, rewrites nested values, or
mutates its input. Constructors normalize admitted scalar grammar before
checking the resulting object. A selected `Gridpoint` is checked when it
materializes, and computation entry points check their complete problem again.

## Validation failures

Construction and mutation check complete values before returning or changing
them. Rechecking at computation entry points detects changes to mutable vectors or
dictionaries retained by an otherwise immutable object.

The failure type communicates the category:

- `ArgumentError` for invalid choices or value kinds.
- `DomainError` for values outside a physical or mathematical domain.
- `DimensionMismatch` for incompatible shapes.
- `MethodError` when no supported operation exists for the supplied types.
- `KeyError` for missing library entries.

## Materialized constructors

Cable-part constructors accept resolved numeric geometry. For example:

```julia
part = Region(:core, Disk(radius), material)
```

Radius or thickness selection, repetition, and variation use explicit
`Grid` inputs. Scalar-complete calls invoke the action directly; varying calls
materialize the same action through `Gridspace`. Completed objects therefore
contain one resolved geometry and cannot drift from their declarations.

Mutable libraries validate a complete candidate before changing owned state.
Earth models, cable designs, and line systems are immutable descriptions;
rebuild the authoritative declaration when it changes.

Operating temperature is not a cable-part constructor input. Cable designs represent
the common material reference state and reject mixed material reference temperatures.
The line problem owns the finite operating temperature. The formulation selects
`temperature_dependence`; its default linear law validates the material-specific
range and positive correction factor. Identity and custom laws retain their own
applicability. [`compute`](@ref) applies the selected
correction without mutating the design.

## Uncertain geometry

A feasible nominal design is not a guarantee about its input distribution.
Internal wire overlaps and crossing layer boundaries remain construction
errors; sampling does not inflate bedding, shrink metal, change counts, or clip
inputs to make them fit. The ring diagnostic reports its count, radius, member
width, required chord, available chord and deficit in meters.

Declare dependent dimensions with a joint builder and justify feasibility over
its complete support. A few successful draws are a useful regression check,
not a support proof. Arbitrary callbacks and independent unbounded normal
inputs cannot be certified from marginal means and standard uncertainties.
`MonteCarlo(...; on_error=:retry, retain_details=true)` explicitly conditions on
successful construction and calculation; it does not preserve the original law.

The existing exterior cable-clearance policy is separate. System construction
may adjust practically touching cable placements with a warning, retaining its
propagated uncertainty reserve across Monte Carlo draws. That policy does not
authorize internal geometric or statistical repairs. Gmsh's tolerances are not
used as manufacturing clearances.

## Reference

```@docs
LineCableModels.validate
```

## Observable resolution and comparison

`observables(...; clip=true)` and `compare` share declared, quantity-aware
resolution in native physical units. No dimensionless epsilon is applied after
unit conversion. This policy is independent of matrix blocks, selected
formulations and visible plot ranges. It is a reporting cutoff, not a certified
bound on solver error or physical uncertainty.

Relative RMS requires both operands above resolution at every selected sample.
If not, it returns `missing`, with the operand counts and reason; it never removes
samples to obtain a favorable subset. Empty bands and unsupported quantities
remain distinct. Absolute RMS always measures the original arrays. A resolved
discrepancy is not evidence that the designated reference is physical truth.

`clip=false` and `observe` retain raw values. Detachment does not erase measurement
uncertainty, break measurement correlations, or clip Monte Carlo standard
deviations. Tests cover prefix/basis scaling, numerical types, threshold equality,
weak signals, complex phase, retained analysis revisions and zero solver calls
during explicit reanalysis.
