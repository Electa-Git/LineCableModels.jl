# Minimal remaining cable fixtures

Defined on 2026-09-28 in [fixtures.jl](fixtures.jl). These three passive physical
constructors complement the existing bare-wire MVP. Loading the file starts no
meshing or computation. They use the current modeling API and normalized native
GetDP voltage extraction; they contain no solver, refinement or acceptance logic.

## Designs

All dimensions below are in millimetres. Metals use the default material-library
records at 20 degrees C; no datasheet resistance correction or equivalent-area
conductor is substituted.

| Constructor | Physical construction | Independent terminals | Purpose |
|---|---|---|---|
| `screened_cable()` | Solid Al core, radius 19.05; 49 explicit Cu wires, radius 0.475 on a 29.025 locus; Al foil, radii 29.9/30.05; PE fill and jacket to radius 32.5 | `core`, `screen`, `foil` | Small strands, one terminal with 49 disconnected metal regions, proximity, dissimilar metals and a 0.15 mm tubular wall |
| `tubular_cable()` | Solid Cu core, radius 23.15; PE insulation to radius 53.25; Pb sheath, radii 53.25/56.55; PE jacket to radius 59.55 | `core`, `sheath` | Thick 3.3 mm wall, both inner/outer surfaces and material-specific skin depth |
| `sector_cable()` | Three Al sectors rotated by 0/120/240 degrees; 1.1 insulation; 30 explicit Cu neutral wires, radius 0.79 on a 14.36 locus; PVC fill/jacket to radius 17.25 | `a`, `b`, `c`, `neutral` | Curved and straight sector boundaries, fillets, conformal insulation, multiple source/receiver pairs and shared neutral |

The first design is reduced from [Tutorial 2](../../../../examples/tutorial2.jl).
It retains the core envelope, complete 49-wire screen and thin Al barrier.
Insulation, screen fill and bedding are one PE region through radius 29.9;
there is no material discontinuity at the wire-envelope radii 28.55/29.5.
The stranded core becomes one solid conductor; semiconductor/tape construction
is replaced with PE over the same radial schedule. In particular, the copper
tape is omitted: this fixture does not claim wire/tape metal-contact coverage.
Wire geometry is explicit and straight; no helical lay correction is prescribed.

The second design retains the core envelope and lead-sheath dimensions from
[Tutorial 3](../../../../examples/tutorial3.jl). Its semiconductor and blocking
layers join the PE insulation. The stranded core becomes solid, and the steel
armor and associated outer layers are omitted. The thick Pb sheath and the
first fixture's thin Al foil provide two wall regimes without a thickness grid.

The third follows the supplied
[sector tutorial](https://github.com/MohamedNumair/LineCableModels.jl/blob/main/examples/tutorial2_sector.jl).
It retains the 119-degree sector opening, 10.24 back radius, 9.14 depth, 1.02
corner radius, insulation, neutral count/locus and cable radius. The current
`Sector` API expresses the base radius as `10.24 - 9.14 = 1.10`. Explicit
placements retain the tutorial's 119-degree conductors at 120-degree pitch;
they do not replace them with 120-degree sectors. PVC is lossless with relative
permittivity 8, as in the example. Remaining gaps receive explicit PVC fill.
The sector sleeves, bedding, neutral fill and jacket have identical PVC
properties. One fill around all 33 metal regions represents their union,
without overlapping offset sleeves or artificial boundaries tangent to the
neutral wires. The earlier separate annular fill failed native Gmsh boundary
layer recovery in a neutral wire at 10 kHz. All metal shapes, terminal
assignments, material admittivities and the outer boundary were verified
identical after removing these redundant dielectric seams. This fixture does
not qualify tangencies between different dielectric materials.
The nominal 95 mm2 designation is not an assertion that the geometric area is
exactly 95 mm2; use the resolved shape area for an independent DC check.

Retrieved tutorial SHA-256:
`348e21445d5fa5673a3ce45c720ed831e98267d3f1cb07394abbbfeeb56e3fb7`.

## Historical qualification selection

The selection below records the completed 36-case campaign; see the
[delivery assessment](delivery-assessment.md). The current
[Julia runner](README.md) defaults to quasi-fw and accepts explicit frequencies
and mesh controls. It does not automatically repeat the old native parity or
source-edit qualification work.

This is a bounded addition after the bare-wire MVP, not a replacement for its
air/earth/mixed placement and resistivity coverage. Use one placement per new
design: screen and tube at `(0,-1)` m, sector at `(0,+1)` m. Use 100 ohm m soil,
relative soil permittivity/permeability 1, temperature 20 degrees C and length
1 m. Keep the exterior/PML preset fixed to the qualified bare-wire preset;
do not coarsen it to conceal conductor cost.

- Frequencies: **0.1, 50, 10 000 and 1 000 000 Hz**. They sample low frequency,
  power frequency, the small-wire/wall transition range and strong skin effect.
  They are not a full sign-crossing search or every region's exact transition.
- Both `quasi_tem` and `quasi_fw`; no bundle, Kron or transposition reduction.
  Give each listed terminal its own connection index, retaining all columns and
  every self/mutual coefficient. Wire count does not become excitation count.
- One candidate sweep: **24 frequency/formulation/design cases, 72 source
  columns**. Reuse the completed local geometry/skin/proximity qualification.
- One explicitly refined high-frequency control per design and formulation:
  **6 cases, 18 columns**. Record componentwise changes without clipping or
  assuming symmetric Y. A poorly resolved component remains unresolved; this
  does not trigger automatic production refinement.
- One same-mesh native-export/Julia comparison at that high frequency for each
  design and formulation: **6 additional solves, 18 columns**. Reuse meshes and
  the candidate results; compare primitive Z/P before final Y. Test native
  material/frequency edits on the same fixtures without a new shape matrix.

Thus the initial bounded selection is **36 cases / 108 source columns**, not
108 independent full-field factorizations. Normal source-column factorization
reuse remains enabled. This count does not promise convergence or include a
new reference campaign. If an integration or accuracy issue remains, identify
the particular missing control before extending the selection.

Use the existing simple serial run/log/resume protocol. This document launches
nothing. No repetition of the 24/48/96-wire or sector-shape grids is required.

## Use and interpretation

```julia
using LineCableModels
include("test/manual/calculations/fem_conductor_mesh/fixtures.jl")
using .ConductorMeshFixtures

design = ConductorMeshFixtures.screened_cable()
connections = Dict(name => i for (i, name) in enumerate(design.terminal_order))
system = build(LineCableSystem, design, Pose2(0.0, -1.0);
    connections, line_length=1.0, system_id=design.cable_id)
problem = LineParametersProblem(system; frequencies=[0.1, 50.0, 1e4, 1e6],
    temperature=20.0, earth_props=homogeneous(rho=100.0, eps_r=1.0, mu_r=1.0))
```

Pass the problem to public `compute` or `export_data(:onelab, ...)` with the
selected formulation and prescribed mesh options. There is no old Python driver
or measurement-preparation stage. Retain source identities, raw Z/P/Y, terminal
ordering, conductor areas, mesh/solver costs and high-frequency current maps.

Use independent round/tube controls and integrated quantities for scientific
qualification. The analytical engine's equivalent-circle or screen models are
useful comparisons, not exact references for these explicit sector/wire
geometries. Keep acceptance criteria in the harness. These fixtures do not
establish general conductor-contact, armor, helical or arbitrary-shape validity.

## Construction checks and native geometry corrections

All three constructors and their `LineParametersProblem` objects were built with
the current Julia 1.12.7 worktree. Metal-region counts are 51, 2 and 33;
independent terminal counts are 3, 2 and 4. The resolved area of each Al sector
is 92.0730502159347 mm2; the entire neutral area is 58.82003925316172 mm2.

A geometry-only `export_data(:onelab, ...)` check used all four frequencies and
eight PML layers to check construction, not numerical accuracy. The screen and
tubular designs exported successfully. The initial sector export failed before meshing:

```text
Circle or ellipse arc 14 greater than Pi (angle=6.28319)
```

The reproduced cause was two contact angles differing by 1.36e-12 radians
whose physical points differed by only about 3e-15 m. The point registry merged
them correctly; the arc builder still requested a same-endpoint arc. It now
omits such zero-length segments, with a focused regression check.

The original fixture also contained intersecting individual PVC sleeves inside
identical PVC bedding. One PVC fill now represents this same material union.
The sector geometry exports successfully with unchanged conductor areas and
radii. This is a fixture representation correction, not a tolerance relaxation
or an engine overlap fallback.

The screen fixture initially split identical PE at each side of the wire layer.
Its tangent interfaces exposed a native boundary-layer meshing failure; changing
Gmsh's geometry tolerance from 1e-8 to 1e-12 or 1e-16 did not resolve it. The
fixture now uses one PE fill around the core and screen, retaining every metal
contour and all actual material interfaces. These fixtures do not qualify
arbitrary touching material constructions.

The sector mesh implementation uses the retained native transfinite construction
and the current conforming voltage paths. Normal/refined full-line numerical
comparisons and high-frequency field maps completed in the bounded batch;
see the delivery assessment for the measured changes and limitations.

Local export evidence: `/tmp/lcm-mesh-fixtures-ytbVCf/`; the bounded check script
is `/tmp/lcm-check-mesh-fixtures.jl`. The normal construction snippet above,
with the original separate-sleeve constructor at `Pose2(0,1)`, followed by the
public export call with `mesh_options=(pml_layers=8,)`, reproduced the old failure.
The current constructor uses the uniform PVC fill described above.
