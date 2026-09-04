# Cable data model

Cable construction follows one path:

```text
physical declaration
        ↓
   @cable / build
        ↓
completed CableDesign
        ↓
  @system / build
        ↓
completed LineCableSystem ──────────┐
                                   ├→ LineParametersProblem → compute → LineParameters
@earth / build → EarthModel ───────┘
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

Place the cable, complete the system declaration, and construct an ordinary
calculation problem:

```julia
system = @system "example-system" line_length=1.0 begin
    @at design (0.0, -1.0) phase=1
end

problem = LineParametersProblem(
    system;
    earth_props=homogeneous(rho=100.0),
    frequencies=[50.0],
)
parameters = compute(problem)
```

System placement composes an outer transform with the design; it does not
rewrite the design's local geometry.

## Declaring earth

`homogeneous(...)` declares a single homogeneous earth. A layered model is
declared from the surface downward with `layer(...)`; the semi-infinite air
layer is implicit:

```julia
earth = @earth begin
    layer(rho=100.0, eps_r=10.0, thickness=5.0)
    layer(rho=500.0, eps_r=20.0)
end
```

Omitting `thickness` makes that earth layer semi-infinite. Consequently, every
finite layer must precede the bottom half-space. Use
`@earth vertical_layers=true` for vertical interfaces, or supply
`air_layer=layer(...)` when the air properties must be explicit. Layer
properties may contain `Grid` values; the result is then a
`Gridspace{EarthModel}` through the same construction path.

The completed `EarthModel` is immutable. Its read-only `layers` tuple begins
with the implicit or explicitly supplied air layer, followed by the declared
earth layers in block order. Programmatic code that already owns a complete
layer collection uses `build(EarthModel, layers; ...)`; earth models are never
extended incrementally with `add!`.

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

`stranded` accepts circular wire shapes. For a `Disk` boundary, `layers` counts
outer courses and the central strand is additional. With `compact=nothing`,
those wires remain circles on their natural touching course radii.

```julia
core_part = stranded(
    copper;
    shape=Disk(0.5e-3),
    lay=LayRatio(13, 12, 11),
    boundary=Disk(4.0e-3),
)
```

`boundary` is the authoritative finished core. The wire radius and boundary
infer the maximum centred triangular-lattice inventory. A lay schedule must
match the radial courses found by that inventory. No `Ring`, layer count, or
wire count is stored as a second description of the same geometry.

Omitting `compact` preserves the natural circular wires. `compact=true`
selects the maximum area-admissible inventory and deforms every strand while
preserving its source area:

```julia
compact_core = stranded(
    copper;
    shape=Disk(0.5e-3),
    compact=true,
    lay=LayRatio(11),
    boundary=Disk(sqrt(7) * 0.5e-3),
)
```

An uncompacted `Sector` boundary has no central wire and needs no artificial
course schedule or wire count. The circular-wire radius and exact boundary
determine the largest feasible polar-shell inventory. Successive shells grow
from the sector base toward its back; each shell follows the available angular
extent, remains symmetric about the bisector, and staggers against its neighbor.
The exact resolved filleted boundary remains authoritative for containment. A
sector phase core has no centre wire.

Course schedules are ordinary physical tuples. Wrap a complete schedule in
`Grid` only when the schedule itself should vary.

A stranded construction may prescribe the aggregate boundary separately from
the geometry of each strand:

```julia
sector_core = stranded(
    copper;
    shape=Disk(0.5e-3),
    boundary=Sector(
        span=deg2rad(119),
        r_base=1.10e-3,
        r_back=10.24e-3,
        fillet=1.02e-3,
    ),
)
```

`shape` remains the circular geometry of every resolved strand. The boundary
adds no homogeneous conductor and no fictitious strand at the cable origin.
Rectangular strands are admitted only in a circular core with an explicit
centre wire. They are always bent area-preservingly into annular courses.

## Independent cores and enclosures

`@assembly` preserves the terminals of independent members:

```julia
phase_a = @terminal :a begin
    core(copper, Sector(span=pi / 3, r_base=1e-3, r_back=4e-3))
    insulation(xlpe; t=1e-3)
end

phase_b = @terminal :b begin
    core(copper, Sector(span=pi / 3, r_base=1e-3, r_back=4e-3))
    insulation(xlpe; t=1e-3)
end

members = @assembly begin
    @at phase_a (-5e-3, 0.0)
    @at phase_b ( 5e-3, 0.0) φ=pi
end
```

`pipe` and `duct` describe containment. Both accept block notation when several
members make the physical nesting clearer:

```julia
contained = @pipe shape=Disk(15e-3) fill=air begin
    @at phase_a (-5e-3, 0.0)
    @at phase_b ( 5e-3, 0.0)
end
```

Nested ducts use the same operation; there is no separate duct-bank object.

When repeated sectors share an origin, their intrinsic span must equal the
angular pitch imposed by the repetition. This keeps neighboring sides parallel.
Bare sectors must retain positive physical side clearance. A conformal outer
insulation may close that clearance to zero, but it may not overlap its
neighbor.

## Canonical grammar

The notation and practical functions lower immediately to the stored grammar:

```text
Material + primitive
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

Primitives describe intrinsic local geometry. Construction resolves them
against a boundary and records absolute geometry in `PlacedRegion` values
inside `CableGeometry`. `EmptyBoundary` represents the first stacking state.

Adding a primitive requires local `resolve`, `boundary`, `area`, `centroid`, and
`support` methods rather than a central type switch.

## Sector cores

`Sector` describes one material-neutral cable-sector cross-section with an
optional fillet. Its symmetry axis is local `+x`; placement and terminal
identity remain outside the primitive.

```julia
sector = Sector(
    span=2pi / 3,
    r_base=1.10e-3,
    r_back=10.24e-3,
    fillet=1.02e-3,
)

phase = terminal(
    :phase,
    Region(:core, sector, copper),
    Region(:insulation, Shell(1.0e-3), xlpe),
)

sectorized = build(
    CableDesign,
    "three-core-sector",
    assembly(
        phase;
        pattern=Ring(3; r=0),
        names=(:a, :b, :c),
    ),
)
```

`resolve` derives exact arc contacts and composes each `Pose2`. `area`,
`perimeter`, `centroid`, and `support` use the exact boundary. `tessellate`
returns ordinary coordinate tuples solely for renderers and mesh adapters.
Resolving `Shell(t)` against one sector produces its exact parallel boundary;
a common layer around the disconnected assembly requires an explicit
`Enclosure`.

The package-owned formulation converts each sector conductor and conformal
dielectric to equivalent-area circles only while preparing its numerical
input. The exact sectors and their centroids remain authoritative in the
completed design and in persisted declarations.

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

Every block macro forwards `combine=:product` or `combine=:zip` to the same
construction used by its functional form. `@cable` also accepts
`nominal_data=(...)`; that descriptive data follows each completed design
without affecting its physical resolution.

Iteration and stochastic realization return completed ordinary designs. Raw
tuples, vectors, polygon points, course schedules, and frequency vectors remain
atomic unless wrapped in `Grid`.

## Formulation boundary

A physically valid design may build even when a formulation does not support
its geometry. `DataModel.flatten(design, frequency)` performs only local scalar circuit reductions:
parallel conductor resistance, recursive GMR, series dielectric admittance,
and the inverse conversions to effective material properties, returning a flat
component payload without changing the design.

`homogenize(design)` is an explicit request for a new homogeneous
`CableDesign`. Independent assembly members are flattened separately; the
operation calculates no mutual coupling, earth return, or line-parameter
matrix. The selected engine separately validates whether its formulation
supports the resulting cable topology.

`CableConstants(design)` consumes this canonical reduction through the
Engine-owned `CableConstantsProblem → CableConstantsFormulation → compute`
workflow. It evaluates one concentric assembly at a time at the requested
temperature and frequency. It contains no earth-return calculation; the
innermost terminal is active and every outward terminal is grounded. The
result stores aligned `R`, `L`, `C`, and `G` vectors with one entry per
assembly.

## Persistence and tables

JSON records authoritative materials, physical declarations, patterns, paths,
compaction laws, poses, terminal names, tags, and explicit Grid declarations.
Decoding invokes the same construction boundary. Resolved geometry, terminal
maps, engine workspaces, and solver results are not serialized.

Model declarations use their bounded `text/plain` displays. Electrical tables
come from `CableConstants(design)` or observed `LineParameters` results.

## Supported conductor formations

The following gallery exercises the public construction grammar directly. Each
example adds only one insulation layer so the conductor formation remains
visible.

### Solid circular and sector cores

A solid conductor may use either a circular disk or the exact filleted sector
primitive. Insulation follows the resolved boundary in both cases.

```@example supported_formations
using LineCableModels
using CairoMakie

gallery_copper = Material(kind=:conductor, rho=1.7241e-8, mu_r=1.0)
gallery_xlpe = Material(kind=:insulator, rho=1.0e14, eps_r=2.3)
gallery_air = Material(kind=:insulator, rho=Inf, eps_r=1.0)
gallery_jacket = Material(kind=:insulator, rho=1.0e13, eps_r=2.8)
gallery_interstitial = Material(kind=:insulator, rho=2.0e12, eps_r=3.5)
gallery_concrete = Material(kind=:insulator, rho=100.0, eps_r=5.0, mu_r=1.0)

solid_circular = build(
    CableDesign,
    "solid-circular",
    terminal(
        :core,
        core(gallery_copper; r=3.5e-3),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

solid_sector = build(
    CableDesign,
    "solid-sector",
    terminal(
        :core,
        Region(
            :core,
            Sector(
                span=deg2rad(118),
                r_base=1.5e-3,
                r_back=7.0e-3,
                fillet=0.6e-3,
            ),
            gallery_copper,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    [solid_circular, solid_sector];
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Circular stranded cores

Omitting `compact` preserves every strand as a circle on its centred
close-packed lattice. `compact=true` instead admits the maximum inventory by
area and deforms it while preserving every strand area and identity.

```@example supported_formations
circular_strand_radius = 0.45e-3
circular_strand_count = 1 + 6 + 12

circular_stranded = build(
    CableDesign,
    "circular-stranded",
    terminal(
        :core,
        stranded(
            gallery_copper;
            shape=Disk(circular_strand_radius),
            lay=LayRatio(14, 12),
            boundary=Disk(5circular_strand_radius),
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

circular_compacted = build(
    CableDesign,
    "circular-compacted",
    terminal(
        :core,
        stranded(
            gallery_copper;
            shape=Disk(circular_strand_radius),
            lay=LayRatio(14, 12),
            compact=true,
            boundary=Disk(sqrt(circular_strand_count) * circular_strand_radius),
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    [circular_stranded, circular_compacted];
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Rectangular stranded core

Rectangular source strands require a circular centre wire. They are always
bent into contiguous annular courses; the radial deformation preserves the
area of every source rectangle. Any residual outer annulus retains the stated
interstitial material.

```@example supported_formations
rectangular_stranded = build(
    CableDesign,
    "rectangular-stranded",
    terminal(
        :core,
        stranded(
            gallery_copper;
            center=Disk(0.55e-3),
            shape=Rectangle(0.75e-3, 0.35e-3),
            lay=LayRatio(13),
            boundary=Disk(3.5e-3),
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    rectangular_stranded;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Sectorized stranded cores

A sectorized core has no central strand. Without `compact`, the wire diameter
and sector boundary determine the largest close-packed inventory; the wires
remain exact circles. With `compact=true`, the maximum area-admissible
inventory is deformed while preserving every strand area and identity.

```@example supported_formations
sector_boundary = Sector(
    span=deg2rad(118),
    r_base=1.5e-3,
    r_back=8.0e-3,
    fillet=0.6e-3,
)
sector_stranded = build(
    CableDesign,
    "sector-stranded",
    terminal(
        :core,
        stranded(
            gallery_copper;
            shape=Disk(0.55e-3),
            boundary=sector_boundary,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

sector_compacted = build(
    CableDesign,
    "sector-compacted",
    terminal(
        :core,
        stranded(
            gallery_copper;
            shape=Disk(0.55e-3),
            lay=LayRatio(14),
            compact=true,
            boundary=sector_boundary,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    [sector_stranded, sector_compacted];
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Three- and four-core sector cables

A complete sector cable preserves the electrical identity of every core. Each
solid or stranded sector therefore owns its insulation before `cores` rotates
the members into a three- or four-phase assembly. The outer `pipe` fills the
interstices and adds one common tubular insulation layer. A sector's angular
span equals the assembly pitch, keeping every pair of facing wedge sides
parallel. The rows below show solid, naturally stranded, and compact-stranded
cores.

```@example supported_formations
three_core_sector = Sector(
    span=2pi / 3,
    r_base=0.60e-3,
    r_back=5.0e-3,
    fillet=0.20e-3,
)
four_core_sector = Sector(
    span=pi / 2,
    r_base=0.60e-3,
    r_back=5.0e-3,
    fillet=0.20e-3,
)

function enclosed_sector_cable(id, member, count)
    phases = cores(
        member;
        n=count,
        r=0,
        names=ntuple(index -> Symbol(:phase_, index), count),
    )
    return build(
        CableDesign,
        id,
        pipe(
            phases;
            shape=Disk(5.5e-3),
            fill=gallery_air,
            wall=insulation(
                gallery_xlpe;
                t=0.8e-3,
                tag=:tubular_insulation,
            ),
        ),
    )
end

solid_three_core_sector = enclosed_sector_cable(
    "solid-three-core-sector",
    terminal(
        :segment,
        solid(gallery_copper, three_core_sector; tag=:solid_sector),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    3,
)
solid_four_core_sector = enclosed_sector_cable(
    "solid-four-core-sector",
    terminal(
        :segment,
        solid(gallery_copper, four_core_sector; tag=:solid_sector),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    4,
)
stranded_three_core_sector = enclosed_sector_cable(
    "stranded-three-core-sector",
    terminal(
        :segment,
        stranded(
            gallery_copper;
            shape=Disk(0.35e-3),
            boundary=three_core_sector,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    3,
)
stranded_four_core_sector = enclosed_sector_cable(
    "stranded-four-core-sector",
    terminal(
        :segment,
        stranded(
            gallery_copper;
            shape=Disk(0.35e-3),
            boundary=four_core_sector,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    4,
)

compacted_three_core_sector = enclosed_sector_cable(
    "compacted-three-core-sector",
    terminal(
        :segment,
        stranded(
            gallery_copper;
            shape=Disk(0.35e-3),
            lay=LayRatio(14),
            compact=true,
            boundary=three_core_sector,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    3,
)
compacted_four_core_sector = enclosed_sector_cable(
    "compacted-four-core-sector",
    terminal(
        :segment,
        stranded(
            gallery_copper;
            shape=Disk(0.35e-3),
            lay=LayRatio(14),
            compact=true,
            boundary=four_core_sector,
            fill=gallery_air,
        ),
        insulation(gallery_xlpe; t=0.25e-3),
    ),
    4,
)

preview(
    [
        solid_three_core_sector,
        solid_four_core_sector,
        stranded_three_core_sector,
        stranded_four_core_sector,
        compacted_three_core_sector,
        compacted_four_core_sector,
    ];
    layout=(3, 2),
    panel_titles=(
        "solid — three cores",
        "solid — four cores",
        "stranded — three cores",
        "stranded — four cores",
        "compacted — three cores",
        "compacted — four cores",
    ),
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Milliken core

A Milliken core owns one centre wire and six separately bounded stranded
segments. `milliken` binds every conductor to one terminal. The segment
boundaries govern packing without becoming overlapping material domains; one
explicit fill owns every interstice around the centre and segment wires. Its
centre-wire radius is inferred after resolving one segment so that the centre
is tangent to the innermost strand in every rotated segment.

```@example supported_formations
milliken_segment_boundary = Sector(
    span=pi / 3,
    r_base=0.45e-3,
    r_back=5.0e-3,
    fillet=0.18e-3,
)
milliken_design = build(
    CableDesign,
    "six-segment-milliken",
    terminal(
        :phase,
        milliken(
            gallery_copper;
            shape=Disk(0.32e-3),
            segment=milliken_segment_boundary,
            segments=6,
            fill=gallery_interstitial,
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    milliken_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Rope formation

`rope` repeats a complete stranded core as one central member surrounded by six
peers. The enclosing terminal and insulation belong to the finished rope.

```@example supported_formations
rope_member_radius = 0.28e-3
rope_member = stranded(
    gallery_copper;
    shape=Disk(rope_member_radius),
    compact=true,
    boundary=Disk(sqrt(7) * rope_member_radius),
    fill=gallery_air,
)

rope_design = build(
    CableDesign,
    "stranded-rope",
    terminal(
        :core,
        rope(
            rope_member;
            layers=1,
            n=6,
            lay=LayRatio(18),
        ),
        insulation(gallery_xlpe; t=1.0e-3),
    ),
)

preview(
    rope_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Three coaxial cores in a pipe

Three independently named solid-core coaxials form a trefoil inside a circular
pipe. The pipe carries one outward insulating wall.

```@example supported_formations
coaxial_member = terminal(
    :phase,
    core(gallery_copper; r=1.6e-3),
    insulation(gallery_xlpe; t=1.2e-3),
)
coaxial_trefoil = cores(
    coaxial_member;
    n=3,
    r=3.5e-3,
    names=(:a, :b, :c),
    φ0=-pi / 2,
)
coaxial_pipe = pipe(
    coaxial_trefoil;
    shape=Disk(7.0e-3),
    fill=gallery_air,
    wall=jacket(gallery_jacket; t=1.0e-3, tag=:pipe_insulation),
)
coaxial_pipe_design = build(CableDesign, "coaxial-trefoil-pipe", coaxial_pipe)

preview(
    coaxial_pipe_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Mixed pipe formations in an elliptical duct

The composition grammar is not restricted to geometries currently admitted by
every numerical backend. Here a pipe containing three stranded coaxials sits
beside a pipe containing three insulated sector-stranded cores, all inside one
elliptical duct. The duct interior remains air; its concrete wall is a finite
normal offset of the ellipse rather than a zero-thickness outline.

```@example supported_formations
duct_round_member = terminal(
    :round,
    stranded(
        gallery_copper;
        shape=Disk(0.30e-3),
        compact=true,
        boundary=Disk(sqrt(7) * 0.30e-3),
        fill=gallery_air,
    ),
    insulation(gallery_xlpe; t=0.8e-3),
)
duct_round_trefoil = cores(
    duct_round_member;
    n=3,
    r=2.8e-3,
    names=(:r1, :r2, :r3),
    φ0=-pi / 2,
)
duct_round_pipe = pipe(
    duct_round_trefoil;
    shape=Disk(6.0e-3),
    fill=gallery_air,
    wall=jacket(gallery_jacket; t=0.8e-3, tag=:round_pipe_insulation),
)

duct_sector_boundary = Sector(
    span=2pi / 3,
    r_base=0.8e-3,
    r_back=3.8e-3,
    fillet=0.3e-3,
)
duct_sector_member = terminal(
    :sector,
    stranded(
        gallery_copper;
        shape=Disk(0.30e-3),
        compact=true,
        boundary=duct_sector_boundary,
        fill=gallery_air,
    ),
    insulation(gallery_xlpe; t=0.6e-3),
)
duct_sector_trefoil = cores(
    duct_sector_member;
    n=3,
    r=0.0,
    names=(:s1, :s2, :s3),
)
duct_sector_pipe = pipe(
    duct_sector_trefoil;
    shape=Disk(5.5e-3),
    fill=gallery_air,
    wall=jacket(gallery_jacket; t=0.8e-3, tag=:sector_pipe_insulation),
)

mixed_duct = duct(
    at(duct_round_pipe, -8.0e-3, 0.0),
    at(duct_sector_pipe, 8.0e-3, 0.0);
    shape=Ellipse(24.0e-3, 10.0e-3),
    fill=gallery_air,
    wall=jacket(gallery_concrete; t=1.5e-3, tag=:concrete_duct_wall),
)
mixed_duct_design = build(CableDesign, "mixed-elliptical-duct", mixed_duct)

preview(
    mixed_duct_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```
