# Construction vocabulary

The practical construction functions lower directly to `Region`, `Stack`,
`Group`, `Assembly`, and `Enclosure`. They carry no state and introduce no
serialized wrapper types.

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
from the centre outward. `solid` and `shell` are the general escape hatches for
physical roles not covered by `core`, `insulation`, `screen`, `sheath`,
`bedding`, `jacket`, and `filler`.

## Wires, strands, and ropes

Declare one fixed wire course with `wires`:

```julia
screen_wires = wires(
    copper;
    wire=Disk(0.5e-3),
    n=40,
    r=20e-3,
    gap_frac=0.02,
    lay=LayRatio(10),
)
```

`strand` creates a central wire and the requested outer courses. A scalar
count is the conventional base count (`k*n` on course `k`); a tuple is an exact
course schedule.

```julia
stranded = strand(
    copper;
    wire=Rectangle(0.35e-3, 0.8e-3),
    layers=3,
    n=(6, 11, 17),
    lay=(LayRatio(13), Pitch(0.15), LayAngle(0.2)),
)
```

`rope(item; ...)` applies the same course grammar to an existing physical
item. Every child retains its internal path declarations; an outer course adds
its own path. `armor` uses the same ring, path, clearance, and compaction
operations while resolving its course radius from the preceding boundary.

`@distribute` inserts only `n=capacity()`:

```julia
automatic = @distribute wires(
    copper;
    wire=Disk(0.5e-3),
    r=20e-3,
    gap_frac=0.03,
)
```

## Tape

One tape function covers conductive, semiconductive, and insulating systems:

```julia
insulating_tapes = @distribute tape(
    tape_insulation;
    section=Sector(8e-3, 8.5e-3, -0.08, 0.16),
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
        section=Sector(8e-3, 8.5e-3, -pi / 4, pi / 2),
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

Homogeneous repeated ducts use a pattern and retain one prototype:

```julia
bank = duct(
    cell;
    formation=Lattice(nx=3, ny=2, dx=12e-3, dy=12e-3),
    shape=Rectangle(44e-3, 30e-3),
    fill=concrete,
)
```

## Formations

Formation functions return ordinary placed-cable declarations:

```julia
placements = @trefoil design spacing=0.09 center=(0.0, -1.0) phase=(1, 2, 3) sheath=0
```

The function forms `trefoil`, `hflat`, and `vflat` accept the same data.
`at(placements, ...)` composes an outer translation and rotation without
copying the cable designs.

## Wire-pattern estimates

[`make_stranded`](@ref) and [`make_screened`](@ref) search deterministic wire
patterns and return a [`WireEstimate`](@ref). An infeasible estimate retains
ranked candidates and states which limits were not met.

```julia
estimate = make_stranded(1000.0)
closest = estimate[:match]
fewest_layers = estimate[:layers]
```

## VDE designation parsing

The qualified [`LineCableModels.DataModel.vdeparse`](@ref) function decodes the
supported VDE/DIN 0271 and 0276 designation fields:

```julia
fields = LineCableModels.DataModel.vdeparse("N2XS(FL)2Y 1x630/35 76/132 kV RM")
```

Unparsed compact-token text remains under `:unparsed_stub`. Parsing a
designation does not construct a cable.

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
