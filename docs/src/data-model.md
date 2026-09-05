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

`stranded` accepts circular wire shapes. For a `Disk` boundary, it infers the
largest complete family of concentric courses: one centre wire and exactly
`6k` wires in course `k`. With `compact=nothing`, those wires remain circles on
their natural course radii.

```julia
core_part = stranded(
    copper;
    shape=Disk(0.5e-3),
    lay=LayRatio(13, 12, 11),
    boundary=Disk(4.0e-3),
)
```

`boundary` is the authoritative finished core. The wire radius and boundary
infer the maximum complete `6k` course inventory. A lay schedule must match the
radial courses found by that inventory. No `Ring`, layer count, or
wire count is stored as a second description of the same geometry.

Omitting `compact` preserves the natural circular wires. `compact=true`
selects the maximum complete `6k` inventory admitted by area and deforms every
strand while preserving its source area:

```julia
compact_core = stranded(
    copper;
    shape=Disk(0.5e-3),
    compact=true,
    lay=LayRatio(11),
    boundary=Disk(sqrt(7) * 0.5e-3),
)
```

A `Sector` boundary infers one bundle-centre strand and complete `6k` courses,
with total inventory `1 + 3L(L+1)`. The circular sites are mapped into the
resolved sector and retained while a prescribed-area power diagram allocates
space. Each strand is a disk grown inside its allocated cell until its clipped
area matches the source wire area. Free portions remain round; contact faces
follow cell boundaries. Sector deformation is intrinsic, so it has no separate
compaction choice. All remaining area belongs to the declared fill material.

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

`shape` defines each source strand and its conserved area; sector reconstruction
changes the resolved outline. The boundary adds no homogeneous conductor and
no fictitious strand at the cable origin.
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

phase = @terminal :phase begin
    solid(copper, sector; tag=:core)
    insulation(xlpe; t=0.5e-3)
end

sectorized = @cable "three-core-sector" begin
    cores(phase; n=3, r=0, names=(:a, :b, :c))
end
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

The following gallery uses the high-level block macros: `@cable` completes a
design, `@terminal` names its conductive descendants, and `@pipe` / `@duct`
express containment. Shapes and leaf operations such as `stranded` remain
ordinary calls inside those blocks. Each example keeps insulation simple so
the conductor formation remains visible.

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

solid_circular = @cable "solid-circular" begin
    @terminal :core begin
        core(gallery_copper; r=3.5e-3)
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

solid_sector = @cable "solid-sector" begin
    @terminal :core begin
        solid(
            gallery_copper,
            Sector(
                span=deg2rad(118),
                r_base=1.5e-3,
                r_back=7.0e-3,
                fillet=0.6e-3,
            );
            tag=:core,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

preview(
    [solid_circular, solid_sector];
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Circular stranded cores

Omitting `compact` preserves every strand as a circle on its concentric `6k`
courses. `compact=true` infers a complete `6k` inventory by area and deforms it
while preserving every strand area and identity.

```@example supported_formations
circular_strand_radius = 0.45e-3
circular_strand_count = 1 + 6 + 12

circular_stranded = @cable "circular-stranded" begin
    @terminal :core begin
        stranded(
            gallery_copper;
            shape=Disk(circular_strand_radius),
            lay=LayRatio(14, 12),
            boundary=Disk(5circular_strand_radius),
            fill=gallery_air,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

circular_compacted = @cable "circular-compacted" begin
    @terminal :core begin
        stranded(
            gallery_copper;
            shape=Disk(circular_strand_radius),
            lay=LayRatio(14, 12),
            compact=true,
            boundary=Disk(sqrt(circular_strand_count) * circular_strand_radius),
            fill=gallery_air,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

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
rectangular_stranded = @cable "rectangular-stranded" begin
    @terminal :core begin
        stranded(
            gallery_copper;
            center=Disk(0.55e-3),
            shape=Rectangle(0.75e-3, 0.35e-3),
            lay=LayRatio(13),
            boundary=Disk(3.5e-3),
            fill=gallery_air,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

preview(
    rectangular_stranded;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Sectorized stranded cores

A sectorized phase core contains its own bundle-centre strand and complete
`6k` courses. The source wire area and sector area determine the inventory.
Mapped circular sites retain these identities; prescribed-area power cells
allocate space, and a clipped-disk area solve determines each copper shape.
The polygons approximate circular arcs while preserving the source areas.

The supplied 0.55 mm wire radius below admits 37 wires and 74.23% copper fill.
The second design explicitly changes the wire radius to obtain 93% fill with
the same inventory and boundary. Compaction cannot remove the remaining void
while preserving both the source wire area and its count.

```@example supported_formations
sector_boundary = Sector(
    span=deg2rad(118),
    r_base=1.5e-3,
    r_back=8.0e-3,
    fillet=0.6e-3,
)
sector_stranded = @cable "sector-stranded" begin
    @terminal :core begin
        stranded(
            gallery_copper;
            shape=Disk(0.55e-3),
            boundary=sector_boundary,
            fill=gallery_air,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

dense_sector_wire_radius = sqrt(
    0.93 * area(LineCableModels.DataModel.resolve(
        LineCableModels.DataModel.EmptyBoundary(), sector_boundary
    )) / (37pi)
)
sector_dense = @cable "sector-stranded-93-percent" begin
    @terminal :core begin
        stranded(
            gallery_copper;
            shape=Disk(dense_sector_wire_radius),
            boundary=sector_boundary,
            fill=gallery_air,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

preview(
    [sector_stranded, sector_dense];
    panel_titles=("supplied wire — 74.23% fill", "resized wire — 93% fill"),
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Three- and four-core sector cables

A complete sector cable preserves the electrical identity of every core. Each
solid or stranded sector therefore owns its insulation before `cores` rotates
the members into a three- or four-phase assembly. The outer `@pipe` fills the
interstices and adds one common tubular insulation layer. A sector's angular
span equals the assembly pitch, keeping every pair of facing wedge sides
parallel. The rows below show solid and intrinsically compacted stranded cores.
At the supplied 0.35 mm strand radius, both stranded cores contain 37 wires;
the three-core sector has 66.62% copper fill and the four-core sector 88.53%.

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

sector_wall = insulation(gallery_xlpe; t=0.8e-3, tag=:tubular_insulation)

solid_three_member = @terminal :segment begin
    solid(gallery_copper, three_core_sector; tag=:solid_sector)
    insulation(gallery_xlpe; t=0.25e-3)
end
solid_three_core_sector = @cable "solid-three-core-sector" begin
    @pipe shape=Disk(5.5e-3) fill=gallery_air wall=sector_wall begin
        cores(solid_three_member; n=3, r=0, names=(:phase_1, :phase_2, :phase_3))
    end
end

solid_four_member = @terminal :segment begin
    solid(gallery_copper, four_core_sector; tag=:solid_sector)
    insulation(gallery_xlpe; t=0.25e-3)
end
solid_four_core_sector = @cable "solid-four-core-sector" begin
    @pipe shape=Disk(5.5e-3) fill=gallery_air wall=sector_wall begin
        cores(solid_four_member; n=4, r=0, names=(:phase_1, :phase_2, :phase_3, :phase_4))
    end
end

stranded_three_member = @terminal :segment begin
    stranded(
        gallery_copper;
        shape=Disk(0.35e-3),
        lay=LayRatio(14),
        boundary=three_core_sector,
        fill=gallery_air,
    )
    insulation(gallery_xlpe; t=0.25e-3)
end
stranded_three_core_sector = @cable "stranded-three-core-sector" begin
    @pipe shape=Disk(5.5e-3) fill=gallery_air wall=sector_wall begin
        cores(stranded_three_member; n=3, r=0, names=(:phase_1, :phase_2, :phase_3))
    end
end

stranded_four_member = @terminal :segment begin
    stranded(
        gallery_copper;
        shape=Disk(0.35e-3),
        lay=LayRatio(14),
        boundary=four_core_sector,
        fill=gallery_air,
    )
    insulation(gallery_xlpe; t=0.25e-3)
end
stranded_four_core_sector = @cable "stranded-four-core-sector" begin
    @pipe shape=Disk(5.5e-3) fill=gallery_air wall=sector_wall begin
        cores(stranded_four_member; n=4, r=0, names=(:phase_1, :phase_2, :phase_3, :phase_4))
    end
end

preview(
    [
        solid_three_core_sector,
        solid_four_core_sector,
        stranded_three_core_sector,
        stranded_four_core_sector,
    ];
    layout=(2, 2),
    panel_titles=(
        "solid — three cores",
        "solid — four cores",
        "stranded — three cores",
        "stranded — four cores",
    ),
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Milliken core

A Milliken core owns one cable-centre wire and six separately bounded stranded
segments, each with its own bundle-centre strand and `6k` courses.
`milliken` binds every conductor to one terminal. The segment
boundaries govern packing without becoming overlapping material domains; one
explicit fill owns every interstice around the centre and segment wires. Its
centre-wire radius is inferred after resolving one segment so that the centre
is tangent to the innermost strand in every rotated segment.

This example deliberately uses a densely filled segment: a source-wire radius
of approximately **0.238628 mm** gives **61 wires per segment**, with the
inventory `1 + 6 + 12 + 18 + 24`, and **93% segment copper fill**. There are
366 segment wires plus the shared central conductor. The latter has an inferred
radius of approximately 0.560 mm; including the inter-segment separators, the
complete core is approximately **84.6% copper by area**. These are illustrative
design dimensions, not measurements recovered from a cable photograph.

The radius is chosen from ``a=\sqrt{\eta A_{\mathrm{segment}}/(61\pi)}``,
with ``\eta=0.93``. Only the radius and boundary are passed to `milliken`:
the constructor still infers the complete-course inventory. Copper fill within
a segment is distinct from copper fill across the complete core, which also
includes the central conductor and the material between segments.

```@example supported_formations
milliken_segment_boundary = Sector(
    span=pi / 3,
    r_base=0.45e-3,
    r_back=5.0e-3,
    fillet=0.18e-3,
)
milliken_design = @cable "six-segment-milliken" begin
    @terminal :phase begin
        milliken(
            gallery_copper;
            shape=Disk(0.23862789195414474e-3),
            segment=milliken_segment_boundary,
            segments=6,
            fill=gallery_interstitial,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

preview(
    milliken_design;
    title="Milliken — six 61-wire segments and a central conductor",
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Four-pair umbilical

An umbilical can contain complete subcables rather than copies of a single
terminal. Here four screened two-core subcables surround a central polyethylene
(PE) filler, with four more PE fillers between them. Each subcable contains two
copper cores with high-density polyethylene (HDPE) insulation, two PE fillers,
an explicit nonwoven (TNT) filling region, a copper sheath and a PE jacket.
The complete assembly has thermoplastic elastomer (TPE) filling, an outer TNT
wrap and a finite galvanized-steel armour layer.

The dimensions and material constants below are **illustrative**, chosen to
reproduce the supplied cross-section's composition without claiming a
manufacturer specification. TNT is represented as a homogeneous insulating
material; the galvanized steel is one effective material, without a separate
zinc coating. This example demonstrates construction and preview, not support
for this topology in every numerical backend.

The ordinary Julia function describes one reusable subcable with the public
macros. Its name prefix gives each copy two independent core terminals and one
sheath terminal. Enclosing them in `@pipe` does **not** short them together.
All dimensions in the code are in metres.

```@example supported_formations
umbilical_pe = Material(kind=:insulator, rho=1.0e15, eps_r=2.3)
umbilical_hdpe = Material(kind=:insulator, rho=1.0e15, eps_r=2.4)
umbilical_tnt = Material(kind=:insulator, rho=1.0e13, eps_r=3.0)
umbilical_tpe = Material(kind=:insulator, rho=1.0e12, eps_r=3.5)
umbilical_steel = Material(kind=:conductor, rho=1.5e-7, mu_r=100.0)

function umbilical_pair(prefix)
    first_core = @terminal Symbol(prefix, :_a) begin
        core(gallery_copper; r=1.5e-3, tag=:copper_core)
        insulation(umbilical_hdpe; t=0.65e-3, tag=:hdpe_insulation)
    end
    second_core = @terminal Symbol(prefix, :_b) begin
        core(gallery_copper; r=1.5e-3, tag=:copper_core)
        insulation(umbilical_hdpe; t=0.65e-3, tag=:hdpe_insulation)
    end
    pair_wall = @terminal Symbol(prefix, :_sheath) begin
        sheath(gallery_copper; t=0.10e-3, tag=:copper_sheath)
        jacket(umbilical_pe; t=0.35e-3, tag=:pe_jacket)
    end
    pair_filler = filler(umbilical_pe, Disk(1.45e-3); tag=:pe_filler)

    return @pipe shape=Disk(4.5e-3) fill=umbilical_tnt wall=pair_wall begin
        @at first_core (-2.25e-3, 0.0)
        @at second_core (2.25e-3, 0.0)
        @at pair_filler (0.0, -2.85e-3)
        @at pair_filler (0.0, 2.85e-3)
    end
end

umbilical_pair_design = @cable "umbilical-subcable" begin
    umbilical_pair(:pair)
end

umbilical_legend_labels = Dict(
    :pipe_fill => "TNT fill", :tpe_fill => "TPE fill",
    :pe_filler => "PE filler", :pe_jacket => "PE jacket",
    :hdpe_insulation => "HDPE insulation", :tnt_wrap => "TNT wrap",
)

preview(
    umbilical_pair_design;
    title="Umbilical subcable — two cores and a separate copper sheath",
    legend_labels=umbilical_legend_labels,
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

The top and bottom subcables are rotated by 90° with `@at`; their complete
geometry rotates together. The four independent subcables retain **eight core
terminals and four sheath terminals**, plus the common armour terminal.
Every remaining region is assigned an insulating material, including the
interstices. No implicit air or unassigned gaps are introduced.

```@example supported_formations
umbilical_filler = filler(umbilical_pe, Disk(2.0e-3); tag=:pe_filler)
umbilical_diagonal = 9.85e-3 / sqrt(2)
umbilical_wall = @terminal :armour begin
    bedding(umbilical_tnt; t=1.5e-3, tag=:tnt_wrap)
    sheath(umbilical_steel; t=2.0e-3, tag=:galvanized_steel)
end

umbilical_design = @cable "four-pair-umbilical" begin
    @pipe shape=Disk(12.5e-3) fill=umbilical_tpe wall=umbilical_wall begin
        @at umbilical_pair(:east) (7.2e-3, 0.0)
        @at umbilical_pair(:north) (0.0, 7.2e-3, pi / 2)
        @at umbilical_pair(:west) (-7.2e-3, 0.0)
        @at umbilical_pair(:south) (0.0, -7.2e-3, pi / 2)
        umbilical_filler
        @at umbilical_filler (umbilical_diagonal, umbilical_diagonal)
        @at umbilical_filler (-umbilical_diagonal, umbilical_diagonal)
        @at umbilical_filler (-umbilical_diagonal, -umbilical_diagonal)
        @at umbilical_filler (umbilical_diagonal, -umbilical_diagonal)
    end
end

preview(
    umbilical_design;
    title="Umbilical — four screened pairs with steel armour",
    legend_group=region -> region.source.material == umbilical_tpe ?
                           :tpe_fill : region.source.tag,
    legend_labels=umbilical_legend_labels,
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

### Rope formation

`rope` repeats a complete stranded core as one central member surrounded by six
peers. Here each member's six course wires have `LayRatio(12)` about their own
member axis. The six outer members also have `LayRatio(18)` about the rope
axis. These are separate paths: the accepted approximation **multiplies their
overlength factors**, rather than replacing the inner lay with the outer one.
The enclosing terminal and insulation belong to the finished rope.

```@example supported_formations
rope_member_radius = 0.28e-3
strand_lay = LayRatio(12)
rope_lay = LayRatio(18)
rope_member = stranded(
    gallery_copper;
    shape=Disk(rope_member_radius),
    lay=strand_lay,
    compact=true,
    boundary=Disk(sqrt(7) * rope_member_radius),
    fill=gallery_air,
)

rope_design = @cable "stranded-rope" begin
    @terminal :core begin
        rope(
            rope_member;
            layers=1,
            n=6,
            lay=rope_lay,
        )
        insulation(gallery_xlpe; t=1.0e-3)
    end
end

preview(
    rope_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```

For local helix radius ``r`` and pitch ``p``, the length factor is
``\lambda=\sqrt{1+(2\pi r/p)^2}``. With `LayRatio(q)`, ``p=2rq``, so
``\lambda(q)=\sqrt{1+(\pi/q)^2}``, independent of the local radius.
For a wire following both paths,

```math
\lambda_{\mathrm{wire}}=\lambda_{\mathrm{strand}}\lambda_{\mathrm{rope}}.
```

The central wire of the central member has neither path and factor 1. Its six
neighbors have only the strand factor; the central wire of each outer member
has only the rope factor. The remaining 36 wires have both. The table below
evaluates these factors using [`overlength`](@ref). Its radius argument is
immaterial for these `LayRatio` laws; a `Pitch` law would require the actual
local radius of each path.

```@example supported_formations
strand_factor = overlength(Helix(strand_lay), rope_member_radius)
rope_factor = overlength(Helix(rope_lay), rope_member_radius)

using DataFrames
DataFrame(
    wires=["central wire", "central member: course", "outer members: centres",
           "outer members: courses"],
    count=[1, 6, 6, 36],
    strand=[1.0, strand_factor, 1.0, strand_factor],
    rope=[1.0, 1.0, rope_factor, rope_factor],
    total=[1.0, strand_factor, rope_factor, strand_factor * rope_factor],
)
```

These factors multiply each wire's resistance before the parallel reduction;
one blanket factor is not applied to all 49 wires. The transverse wire areas
and preview are unchanged by the lay choices.

### Three coaxial cores in a pipe

Three independently named solid-core coaxials form a trefoil inside a circular
pipe. The pipe carries one outward insulating wall.

```@example supported_formations
coaxial_member = @terminal :phase begin
    core(gallery_copper; r=1.6e-3)
    insulation(gallery_xlpe; t=1.2e-3)
end
coaxial_trefoil = cores(
    coaxial_member;
    n=3,
    r=3.5e-3,
    names=(:a, :b, :c),
    φ0=-pi / 2,
)
coaxial_wall = jacket(gallery_jacket; t=1.0e-3, tag=:pipe_insulation)
coaxial_pipe_design = @cable "coaxial-trefoil-pipe" begin
    @pipe shape=Disk(7.0e-3) fill=gallery_air wall=coaxial_wall begin
        coaxial_trefoil
    end
end

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
duct_round_member = @terminal :round begin
    stranded(
        gallery_copper;
        shape=Disk(0.30e-3),
        compact=true,
        boundary=Disk(sqrt(7) * 0.30e-3),
        fill=gallery_air,
    )
    insulation(gallery_xlpe; t=0.8e-3)
end
duct_round_trefoil = cores(
    duct_round_member;
    n=3,
    r=2.8e-3,
    names=(:r1, :r2, :r3),
    φ0=-pi / 2,
)
duct_round_wall = jacket(gallery_jacket; t=0.8e-3, tag=:round_pipe_insulation)
duct_round_pipe = @pipe shape=Disk(6.0e-3) fill=gallery_air wall=duct_round_wall begin
    duct_round_trefoil
end

duct_sector_boundary = Sector(
    span=2pi / 3,
    r_base=0.8e-3,
    r_back=3.8e-3,
    fillet=0.3e-3,
)
duct_sector_member = @terminal :sector begin
    stranded(
        gallery_copper;
        shape=Disk(0.30e-3),
        boundary=duct_sector_boundary,
        fill=gallery_air,
    )
    insulation(gallery_xlpe; t=0.6e-3)
end
duct_sector_trefoil = cores(
    duct_sector_member;
    n=3,
    r=0.0,
    names=(:s1, :s2, :s3),
)
duct_sector_wall = jacket(gallery_jacket; t=0.8e-3, tag=:sector_pipe_insulation)
duct_sector_pipe = @pipe shape=Disk(5.5e-3) fill=gallery_air wall=duct_sector_wall begin
    duct_sector_trefoil
end

duct_wall = jacket(gallery_concrete; t=1.5e-3, tag=:concrete_duct_wall)
mixed_duct_design = @cable "mixed-elliptical-duct" begin
    @duct shape=Ellipse(24.0e-3, 10.0e-3) fill=gallery_air wall=duct_wall begin
        @at duct_round_pipe (-8.0e-3, 0.0)
        @at duct_sector_pipe (8.0e-3, 0.0)
    end
end

preview(
    mixed_duct_design;
    backend=:cairo,
    display_plot=false,
    controls=false,
).figure
```
