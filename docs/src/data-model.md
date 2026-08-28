# Cable data model

The stored physical grammar is:

```text
Material + Primitive
        ↓
      Region
        ↓
Stack / Group / Assembly / Enclosure
        ↓
      build
        ↓
completed CableDesign
        ↓
completed LineCableSystem
        ↓
LineParametersProblem → compute → LineParameters
```

`Material{T}` stores constitutive values and `kind::Symbol`. A primitive stores
intrinsic geometry. `Region` binds one primitive, one material, and one stable
physical tag.

`Stack` composes parts in outward order. `Group` repeats a part while assigning
one retained terminal. `Assembly` repeats a part under separate member names.
`Enclosure` owns a containing shape, fill, and optional wall. Placement patterns,
helical paths, and compaction values remain ordinary declarations within this
grammar.

## Complete construction

`build` is the public action that turns a complete physical declaration into a
completed domain object:

```julia
copper = Material(kind=:conductor, rho=1.7241e-8)
xlpe = Material(kind=:insulator, rho=1.0e14, eps_r=2.3)

root = Stack(
    Group(:phase, Conductor.Solid(:core, copper; r=10e-3)),
    Insulator.Shell(:insulation, xlpe; t=8e-3),
)

design = build(CableDesign, "example", root)
system = build(
    LineCableSystem,
    design,
    Pose2(0.0, -1.0);
    connections=(phase=1,),
    line_length=1000.0,
)
problem = LineParametersProblem(
    system;
    earth_props=Earth(rho=100.0),
    frequencies=[50.0],
)
parameters = compute(problem)
```

A completed `CableDesign` stores its cable identifier, nominal data,
authoritative physical root, `CableGeometry`, terminal order, and terminal map.
It does not store an analytical equivalent, reference frequency, engine input,
or solver matrices.

`build(CableDesign, ...)` validates the physical declaration, calls `resolve`
through the primitive and composition dispatch tree, constructs one
`CableGeometry`, assigns terminals, and freezes the derived ordering with the
root. `build(LineCableSystem, ...)` places completed designs and resolves global
terminal and connection order. Neither action runs a formulation.

## Intrinsic and resolved geometry

`AbstractPrimitive` and `AbstractShape` have separate responsibilities. A
`Disk`, `Rectangle`, or other primitive is intrinsic unresolved geometry.
`DiskShape`, `RectangleShape`, and the other shape types contain resolved
absolute geometry. `EmptyBoundary` is the explicit initial state for outward
stacking.

`CableGeometry` contains ordered `PlacedRegion` values and one outer resolved
shape. A `PlacedRegion` retains its source `Region`, resolved shape, terminal,
placement metadata, and longitudinal path declarations. It contains no
formulation result.

Adding a primitive requires local `resolve`, `boundary`, `area`, `centroid`, and
`support` methods. Construction does not use a central primitive switch.

## Stranded conductors

Circular and rectangular strands use the same convenience surface:

```julia
circular = Conductor.Stranded(
    :core, Disk(0.5e-3), copper;
    layers=4, n=6, lay=LayRatio(11),
)

rectangular = Conductor.Stranded(
    :core, Rectangle(0.35e-3, 0.8e-3), copper;
    layers=4, n=6, lay=LayRatio(11),
)
```

Both calls return an ordinary `Stack` containing one `Group` per physical
layer. Rectangular members retain their width, height, area, pose, and terminal
identity. Their layer loci follow the preserved annular-area radial relation
and enforce the tangential packing limit. No dedicated stranded-conductor
storage type is introduced.

## Parameter spaces

With no direct `Grid` or `Gridspace` argument, `build` returns one completed
object. With an explicit finite source, the same surface returns a
`Gridspace{Target}` whose callable invokes scalar `build` after resolving one
point:

```julia
identifiers = Grid(("design-a", "design-b"))
designs = build(CableDesign, identifiers, root)

eltype(designs) === CableDesign
first(designs) isa CableDesign
```

Raw tuples, vectors, polygons, stack members, and frequency vectors remain
atomic unless explicitly wrapped in `Grid`.

## Formulation boundary

A physical design may build successfully even when a formulation does not
support its geometry. The analytical adapter consumes `CableGeometry` and
performs its shape-, material-, path-, and terminal-specific preparation there.
Unsupported shapes fail before numerical work. Construction never substitutes
an equivalent circle or stores equivalent resistance or GMR as physical truth.

`CableConstants(design)` follows the ordinary problem path. It assigns the
selected core terminal to phase 1, assigns other terminals to phase 0, builds a
`LineParametersProblem`, calls `compute`, and observes `R`, `L`, and `C` from the
result.

## Persistence and tables

The JSON format records materials, physical roots, primitives, patterns, paths,
placements, connections, identifiers, nominal data, and explicit Grid
declarations. Decoding invokes `build`. Resolved geometry, terminal maps,
analytical inputs, and solver results are not serialized.

`DataFrame(design)` reports one row per physical placed region. Electrical
tables come from `DataFrame(CableConstants(design))` or observed
`LineParameters` results.
