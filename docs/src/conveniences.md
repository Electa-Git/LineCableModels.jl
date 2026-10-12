# Construction vocabulary

The practical construction functions lower directly to `Region`, `Stack`,
`Group`, `Assembly`, and `Enclosure`. These stateless functions omit serialized wrapper types.

## Regions and ordered layers

```julia
phase = terminal(
    :phase,
    core(copper; r=8e-3),
    screen(semicon; t=0.8e-3),
    insulation(xlpe; t=6e-3),
)
```

Arguments to `terminal`, `layers`, and `build(CableDesign, ...)` are ordered
from the center outward. `solid` and `shell` are the general escape hatches for
physical roles not covered by `core`, `insulation`, `screen`, `sheath`,
`bedding`, `jacket`, and `filler`.

## Wires, wire strands, and ropes

Declare one fixed wire course with `wires`:

```julia
screen_wires = wires(
    copper;
    shape=Disk(0.5e-3),
    n=40,
    r=20e-3,
    gap_frac=0.02,
    lay=LayRatio(10),
)
```

With a `Disk` geometric boundary, `stranded` creates a center wire and infers the maximum
complete `6k` course inventory. Without compaction circular wires retain their
natural shape. `compact=true` requests area-preserving deformation of those
same complete courses. A `Sector` geometric boundary maps a center strand and complete
`6k` courses into the sector, then reconstructs each strand as an
area-preserving clipped disk. Its deformation is intrinsic. It has no
compaction selector or separately supplied center wire.

```julia
compacted_core = stranded(
    copper;
    shape=Disk(0.5e-3),
    lay=(LayRatio(13), Pitch(0.15), LayAngle(0.2)),
    compact=true,
    boundary=Disk(sqrt(37) * 0.5e-3),
)
```

A schedule sharing one definition lifts its raw course values directly:

```julia
round_core = stranded(
    copper;
    shape=Disk(0.5e-3),
    lay=LayRatio(13, 12, 11),
    boundary=Disk(3.5e-3),
)
```

The explicit tuple form above remains available when course definitions differ.

`rope(item; ...)` applies the same course grammar to an existing physical
item. Every child retains its internal path declarations. An outer course adds
its own path. `armor` uses the same ring, path, clearance, and compaction
operations while resolving its course radius from the preceding geometric boundary.

`@distribute` inserts only `n=capacity()`:

```julia
automatic = @distribute wires(
    copper;
    shape=Disk(0.5e-3),
    r=20e-3,
    gap_frac=0.03,
)
```

## Tape

One tape function covers conductive, semiconductive, and insulating systems:

```julia
insulating_tapes = @distribute tape(
    tape_insulation;
    section=Rectangle(1.4e-3, 0.5e-3),
    gap_frac=0.02,
    lay=LayRatio(10),
)
```

This group contributes no terminal. A conductive tape becomes a terminal only
under an explicit terminal scope:

```julia
screen_tape = @terminal :screen begin
    tape(
        copper;
        section=Rectangle(3.0e-3, 0.5e-3),
        n=4,
        lay=LayRatio(12),
    )
end
```

## Independent members and ducts

`cores(item; n, r, names)` retains one prototype and a repetition pattern.
`cores(member1, member2, ...)` and `assembly(member1, member2, ...)` preserve
explicit heterogeneous members and their poses.

```julia
cells = @assembly begin
    @at cell_a (-12e-3, 0.0)
    @at cell_b ( 12e-3, 2e-3) φ=0.1
end

bank = @duct shape=Rectangle(40e-3, 24e-3) fill=concrete begin
    cells
end
```

The corresponding pipe notation uses the same block grammar:

```julia
contained = @pipe shape=Disk(20e-3) fill=air begin
    cells
end
```

A duct bank is a concrete duct whose members are air-filled ducts. Each air-filled duct
contains one cable and is placed with `@at`:

```@example duct_bank
using LineCableModels
using CairoMakie

copper = Material(kind=:conductor, rho=1.7241e-8)
xlpe = Material(kind=:insulator, rho=1.0e14, eps_r=2.3)
air = Material(kind=:insulator, rho=Inf, eps_r=1.0)
concrete = Material(kind=:insulator, rho=100.0, eps_r=5.0)

cell(name) = @duct shape=Disk(5e-3) fill=air begin
    @terminal name begin
        core(copper; r=2e-3)
        insulation(xlpe; t=1.5e-3)
    end
end

bank = @duct shape=Rectangle(44e-3, 30e-3) fill=concrete begin
    @at cell(:a) (-12e-3, -6e-3)
    @at cell(:b) (  0.0,  -6e-3)
    @at cell(:c) ( 12e-3, -6e-3)
    @at cell(:d) (-12e-3,  6e-3)
    @at cell(:e) (  0.0,   6e-3)
    @at cell(:f) ( 12e-3,  6e-3)
end

duct_bank = @cable "duct-bank" begin
    bank
end

preview(duct_bank; backend=:cairo, display_plot=false, controls=false).figure
```

## Formations

Formation functions return ordinary placed-cable declarations:

```julia
placements = @trefoil design spacing=0.09 center=(0.0, -1.0) phase=(1, 2, 3) sheath=0
```

Nonzero one-based values identify active phases. `0` marks a conductor selected
for grounded or eliminated-conductor reduction. The values are identifiers, not
electrical polarity or phase-angle signs.

The function forms `trefoil`, `hflat`, and `vflat` accept the same data. A
`spacing` below the outer diameter plus 1 μm is raised to it for each realization of the
design. The cables then touch, 1 μm apart, and never overlap.
`at(placements, ...)` composes an outer translation and rotation without
copying the cable designs.

## Wire-pattern estimates

[`estimate_stranding`](@ref) and [`estimate_screen`](@ref) search deterministic wire
patterns and return a [`WireEstimate`](@ref). An infeasible estimate retains
ranked patterns and states which limits were not met.

```julia
estimate = estimate_stranding(1000.0)
closest = estimate[:closest_area]
fewest_layers = estimate[:fewest_layers]
```

## VDE designation parsing

The qualified [`LineCableModels.DataModel.vdeparse`](@ref) function decodes the
supported VDE/DIN 0271 and 0276 designation fields:

```julia
fields = LineCableModels.DataModel.vdeparse("N2XS(FL)2Y 1x630/35 76/132 kV RM")
```

The returned fields describe the designation. Unparsed compact-token text is
stored under `:unparsed_stub`.

## Reference

```@docs
LineCableModels.DataModel.vdeparse
Base.get
Base.delete!
```

```@autodocs
Modules = [
    LineCableModels.DataModel.BaseParams,
    LineCableModels.ParametricBuilder.WirePatterns,
]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = false
```
