# Cable data model

Cable construction follows one path:

```text
physical declaration
        ↓
      build
        ↓
completed CableDesign
        ↓
LineParametersProblem → compute → LineParameters
```

The notation layer describes physical order and electrical ownership. It does
not introduce another stored model.

## Constructing a cable

The following declaration describes one coaxial cable. Expressions inside
`@cable` run from the centre outward. `@terminal` assigns one electrical name
to every conductive descendant in its block.

```julia
copper = Material(kind=:conductor, rho=1.7241e-8)
xlpe = Material(kind=:insulator, rho=1.0e14, eps_r=2.3)

design = @cable "example" begin
    @terminal :phase begin
        core(copper; r=10e-3)
        insulation(xlpe; t=8e-3)
    end
end
```

`design` is a completed `CableDesign`; no second resolution call is required.
Its physical declaration remains authoritative while geometry and terminal
indices are derived once by `build`.

Place the cable and construct an ordinary calculation problem:

```julia
placed = @at design (0.0, -1.0) connections=(phase=1,)
problem = LineParametersProblem(
    [placed];
    earth_props=Earth(rho=100.0),
    frequencies=[50.0],
)
parameters = compute(problem)
```

System placement composes an outer transform with the design; it does not
rewrite the design's local geometry.

## Terminals and physical tags

Terminal names and physical tags have different jobs:

```julia
@terminal :a begin
    core(copper; r=4e-3, tag=:core)
    insulation(xlpe; t=2e-3, tag=:insulation)
end
```

`:a` is the retained electrical terminal. `:core` and `:insulation` identify
physical regions locally. The same local tags may be used under another
terminal because qualified identities are derived from the tree.

A nonconductive repeated group is valid and contributes no terminal. A
conductive screen becomes a separate terminal only when requested explicitly:

```julia
@terminal :screen begin
    sheath(copper; t=0.5e-3, tag=:metallic_screen)
end
```

## Repeated conductors

`strand` accepts circular or rectangular wire sections through the same
surface. `layers` counts outer courses; the central wire is additional.

```julia
core_part = strand(
    copper;
    wire=DiskDefinition(0.5e-3),
    layers=3,
    n=(6, 12, 18),
    lay=(LayRatio(13), LayRatio(12), LayRatio(11)),
)
```

Use `capacity()` when a course should contain the maximum count allowed by its
actual geometry and compaction law:

```julia
compact_core = @distribute strand(
    copper;
    wire=DiskDefinition(0.5e-3),
    layers=1,
    compact=FillFactor(0.9),
    lay=LayRatio(11),
)
```

Course schedules are ordinary physical tuples. Wrap a complete schedule in
`Grid` only when the schedule itself should vary.

## Independent cores and enclosures

`@assembly` preserves the terminals of independent members:

```julia
phase_a = @terminal :a begin
    core(copper, SectorDefinition(0, 4e-3, -pi / 6, pi / 3))
    insulation(xlpe; t=1e-3)
end

phase_b = @terminal :b begin
    core(copper, SectorDefinition(0, 4e-3, -pi / 6, pi / 3))
    insulation(xlpe; t=1e-3)
end

members = @assembly begin
    @at phase_a (-5e-3, 0.0)
    @at phase_b ( 5e-3, 0.0) φ=pi
end
```

`pipe` and `duct` describe containment. Nested ducts use the same operation;
there is no separate duct-bank object.

## Canonical grammar

The notation and practical functions lower immediately to the stored grammar:

```text
Material + primitive definition
              ↓
            Region
              ↓
Stack / Group / Assembly / Enclosure
              ↓
            build
              ↓
         CableDesign
```

`Stack` owns ordered outward composition. `Group` repeats and coalesces
conductive descendants into zero or one terminal. `Assembly` preserves child
terminal identities. `Enclosure` owns one containing domain, its fill, its
contents, and an optional wall. These types remain available for extension
code and unusual geometry; ordinary cable declarations use the vocabulary
shown above.

Primitive definitions describe intrinsic local geometry. Construction resolves
them against a boundary and records absolute geometry in `PlacedRegion` values
inside `CableGeometry`. `EmptyBoundary` represents the first stacking state.
Adding a primitive requires local `resolve`, `boundary`, `area`, `centroid`, and
`support` methods rather than a central type switch.

## Parameter spaces

Only an explicit `Grid` introduces variation:

```julia
designs = @cable "radius-study" begin
    @terminal :phase begin
        core(copper; r=Grid((8e-3, 10e-3, 12e-3)))
        insulation(xlpe; t=8e-3)
    end
end

eltype(designs) === CableDesign
```

Iteration and stochastic realization return completed ordinary designs. Raw
tuples, vectors, polygon points, course schedules, and frequency vectors remain
atomic unless wrapped in `Grid`.

## Formulation boundary

A physically valid design may build even when a formulation does not support
its geometry. The selected engine validates and adapts `CableGeometry` before
constructing its input. Construction does not substitute equivalent circles,
effective radii, or homogeneous conductors.

`equivalent(design)` is an explicit request for a new homogeneous
`CableDesign`; it is not part of ordinary problem construction.

`CableConstants(design)` follows the same problem and computation path as line
parameters, then observes `R`, `L`, and `C` from the result.

## Persistence and tables

JSON records authoritative materials, physical declarations, patterns, paths,
compaction laws, poses, terminal names, tags, and explicit Grid declarations.
Decoding invokes the same construction boundary. Resolved geometry, terminal
maps, engine workspaces, and solver results are not serialized.

`DataFrame(design)` reports resolved physical regions. Electrical tables come
from `CableConstants(design)` or observed `LineParameters` results.
