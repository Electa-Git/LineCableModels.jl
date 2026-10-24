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

Constructors normalize admitted scalar grammar and fill defaults before
validating the resulting object. A selected `Gridpoint` is checked when it
materializes. A computation does not check its problem again.

## Validation failures

Construction and mutation check complete values before returning or changing
them. Calling `validate` on a value checks it again, which detects changes to mutable
vectors or dictionaries retained by an otherwise immutable object.

The failure type communicates the category:

- `ArgumentError` for invalid choices or value kinds.
- `DomainError` for values outside a physical or mathematical domain.
- `DimensionMismatch` for incompatible shapes.
- `MethodError` when no supported operation exists for the supplied types.
- `KeyError` for missing library entries.

## Materialized constructors

Cable-part constructors accept resolved numeric geometry. This call constructs
a circular core from its radius and material:

```julia
part = Region(:core, Disk(radius), material)
```

Radius or thickness selection, repetition, and variation use explicit
`Grid` inputs. Scalar-complete calls invoke the action directly. Varying calls
materialize the same action through `Gridspace`. Completed objects
contain one resolved geometry and cannot drift from their declarations.

Mutable libraries validate a complete library replacement before changing owned state.
Earth models, cable designs, and line systems are immutable descriptions.
Rebuild the authoritative declaration when it changes.

Cable designs represent a common material reference state and reject mixed
material reference temperatures. The line problem defines the finite operating
temperature. The formulation selects
`temperature_dependence`. Its default linear law validates the material-specific
range and positive correction factor. Identity and custom laws retain their own
applicability. [`compute`](@ref) applies the selected
correction without mutating the design.

## Uncertain geometry

A feasible nominal design is not a guarantee about its input distribution.
Internal wire overlaps and crossing layer boundaries remain construction
errors. Sampling does not inflate bedding, shrink metal, change counts, or clip
inputs to make them fit. The ring diagnostic reports its count, radius, member
width, required chord, available chord and deficit in meters.

Declare dependent dimensions with a joint builder and justify feasibility over
its complete support. Successful draws establish feasibility for those
realizations only. Arbitrary callbacks and independent unbounded normal
inputs cannot be certified from marginal means and standard uncertainties.
`MonteCarlo` draws again a draw whose build or computation throws `DomainError`, so its
output distribution is conditional on successful construction and calculation.

Exterior cable clearance follows its own rule. A formation lays its cables at least
1 μm apart, using the outer radius of each realization of its design, so it never
overlaps. A system build raises a gap from zero up to 1 μm, between two cables or from a
cable to the air-earth interface, to 1 μm. A larger overlap throws
[`OverlapError`](@ref). The system checks its pair gaps, and a line problem checks the
interface gaps. Gaps use nominal values, so uncertainty never moves the nominal geometry.
`MonteCarlo` discards a draw that throws `OverlapError` and draws again. Internal
geometric and statistical repairs remain prohibited.

## Observable resolution and comparison

`observables(...; clip=true)` and `compare` share declared, quantity-aware
resolution in native physical units. No dimensionless epsilon is applied after
unit conversion. This reporting rule is independent of matrix blocks, selected
formulations and visible plot ranges. A certified bound on solver error or physical uncertainty cannot be inferred from this reporting cutoff.

Relative RMS requires both operands above resolution at every selected sample.
If not, it returns `missing`, with the operand counts and reason. It never removes
samples to obtain a favorable subset. Empty bands and unsupported quantities
remain distinct. Absolute RMS always measures the original arrays. A resolved
discrepancy is not evidence that the designated reference is physical truth.

`clip=false` and `observe` retain raw values. Detachment preserves measurement
uncertainty, measurement correlations, and Monte Carlo standard deviations.
