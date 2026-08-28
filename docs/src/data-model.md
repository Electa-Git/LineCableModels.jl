# Cable data model

The physical hierarchy is:

```text
Material → Primitive → Region → Stack / Group / Assembly / Enclosure
         → CableDesign → LineCableSystem → Engine input
```

Each name has one owner and one job. `Material` stores constitutive properties
and a physical class such as `:conductor`, `:insulator`, or `:semicon`.
Primitives describe intrinsic cross-section geometry. A `Region` associates one
primitive with one material and a stable tag.

`Stack` composes parts outward in order. `Group` repeats one part while retaining
one electrical terminal. `Assembly` repeats a part while retaining independent
member identities. `Enclosure` owns containing geometry, fill, and an optional
wall. Rings, lattices, helical paths, and compaction laws are concrete values
used by these compositions; they are not alternate cable records.

## Eager authorities

`CableDesign(root)` resolves geometry, terminal membership, terminal order, and
the current analytical equivalence in one constructor call. The supplied
physical root is authoritative. All other fields are derived and immutable.

`LineCableSystem(designs; positions, connections, ...)` places eager designs,
resolves global terminal order and connections, and freezes the derived system
geometry. It does not run a formulation. The Engine translates a complete
system into its existing analytical input in one direction.

```julia
copper = Material(kind=:conductor, rho=1.7241e-8)
xlpe = Material(kind=:insulator, rho=1.0e14, eps_r=2.3)

root = Stack(
    Group(:phase, Conductor.Solid(:core, copper; r=10e-3)),
    Insulator.Shell(:insulation, xlpe; t=8e-3),
)
design = CableDesign(root; cable_id="example")
system = LineCableSystem(
    design;
    position=Pose2(0.0, -1.0),
    connections=(phase=1,),
    line_length=1000.0,
)
```

## Parameter spaces

Ordinary scalar input constructs the eager target immediately. An explicit
`Grid` or nested `Gridspace` returns `Gridspace{Target}` and delays only finite
selection. Every selected point invokes the same target constructor:

```julia
radii = Grid((9e-3, 10e-3))
cores = Conductor.Solid(:core, copper; r=radii)
designs = CableDesign(Group(:phase, cores); cable_id="varying-core")

first(designs) isa CableDesign
```

There is no cable specification or builder object between a selected point and
the eager domain object.

## Geometry and analytical support

Intrinsic primitives and resolved geometry are distinct. `resolve` applies
context, placement, paths, and compaction without rewriting the declaration.
A physically valid geometry may still be unsupported by a formulation.
`CableDesign` and `LineCableSystem` then construct normally, while the Engine
rejects the unsupported geometry before analytical input construction. The
adapter never substitutes a circular shape, dummy fill, or approximate GMR
silently.

## Persistence

The versioned JSON format uses JSON Schema Draft 2020-12. It records materials,
physical declarations, placements, connections, and explicit parameter grids.
It does not record resolved geometry, terminal maps, analytical equivalence,
engine inputs, or solver state. Decoding invokes the same eager constructors,
so imported declarations cannot supply stale derived state.
