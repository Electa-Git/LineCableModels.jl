# PML resolution and cost: measured diagnostics, 2026-09-29

The current 192-interval prescription has not been established as the minimum
needed for the accepted quasi-fw results. It is a successful retained preset,
not an error-estimator output. Its tensor-product corner construction explains
most of the mesh cost. There is evidence that exterior discretization matters,
but no evidence here that all 400,128 PML triangles are necessary.

This report reads existing completed meshes, inputs and native logs. It launches
no Gmsh/GetDP processes, changes no production code or manual-runner controls,
and leaves the concurrent performance benchmark alone. Count projections below
are arithmetic, not remeshing results, accuracy qualification or speedup claims.

## Reproduction and evidence

Diagnostic directory, relative to the repository:

`.linecablemodels/fem/pml-resolution-audit-20260929/`

Run `python3 .linecablemodels/fem/pml-resolution-audit-20260929/analyze.py` from
the repository root. It needs NumPy and reads ASCII MSH 4.1 directly. It checks
the ten mesh SHA-256 values against the completed native-attempt records. It
counts each triangle once, even when a surface belongs to several physical
groups. The existing `model.msh` is the tenth frequency; the other meshes are
`frequency_0001.msh` through `frequency_0009.msh`.

The primary completed scan is `.linecablemodels/fem/runs/run-vEE4tu/`:

- Two copper disks of radius 0.0425 m, centres (0, +1) and (1, -1) m.
- Flat interface y=0; upper air, lower earth, resistivity 100 ohm m.
- Both surrounding media have relative permittivity and permeability 1.
- Ten frequencies from 0.1 Hz to 1 MHz; quasi-fw; both source columns.
- `domain_skin_depths=24`, `pml_layers=192`, `mesh_size_factor=3`,
  `exterior_mesh_size_factor=8`, `pml_reflection=1e-10`.
- PML thickness equals physical half-width in all directions.
- Conductor controls: area tolerance 1e-3, six elements per skin depth,
  normal growth sqrt(1.25), five skin depths, four thickness elements.
- Gmsh metadata: 4.15.2-git. Native executable: GetDP 3.5.0 Linux64.
- First-order triangles; nodal axial potential, nodal scalar potential,
  edge transverse potential with tree gauge; 12-point triangle integration.
- PETSc/MUMPS unsymmetric direct LU, diagonal scaling, two frequency workers,
  four solver threads. One factorization per frequency, reused for source 2.
- Field maps were disabled; actual PML field amplitudes cannot be inferred
  from these saved meshes and raw terminal results alone.

The directory contains:

| File | Contents |
|---|---|
| `mesh-10-frequencies.csv` | Measured disjoint regional counts at every frequency |
| `surface-blocks.csv` | Every surface tag, bounding box, area, triangle count and corner cell aspect ratio |
| `solver-10-frequencies.csv` | DOFs, assembly/solve times, reported memory and sparse matrix storage |
| `parameters-70-points.csv` | Resolved geometry, strengths and grading for all seven recent scans |
| `normal-wave-cell-diagnostics.csv` | Per-cell phase/attenuation increments for three analytic normal-wave probes |
| `count-projections-not-solves.csv` | Counts for alternative N, holding physical mesh and tangential counts fixed |
| `mesh-sha256.csv` | Paths and hashes of the ten original meshes |
| `source-snapshot/` | Current geometry/model code, retained solver sources and primary resolved inputs |
| `pml-resolution.png`, `.svg` | Measured block counts, grading and arithmetic projections |

Other parameter records cover rho=0.1, 1, 1000 ohm m at the standard radius,
and radii 0.001, 0.01, 0.085 m at rho=0.1. Their mesh files were not recounted
in this audit; the ten measured meshes are all the 100 ohm m primary case.

## Exactly where the triangles go

| Disjoint region | 0.1 Hz triangles | 1 MHz triangles |
|---|---:|---:|
| Physical air | 16,197 | 14,515 |
| Physical earth | 17,101 | 14,897 |
| Both metal interiors | 626 | 6,004 |
| Left/right PML strips, air halves | 36,864 | 36,864 |
| Left/right PML strips, earth halves | 36,864 | 36,864 |
| Top PML, excluding corners | 15,744 | 15,744 |
| Bottom PML, excluding corners | 15,744 | 15,744 |
| Two air corners | 147,456 | 147,456 |
| Two earth corners | 147,456 | 147,456 |
| **Total** | **434,052** | **435,544** |
| **PML subtotal** | **400,128 (92.18%)** | **400,128 (91.87%)** |

Across all ten frequencies the PML count is exactly 400,128. Total triangle
counts range from 433,602 to 437,398, so the PML share is 91.48–92.28%.
Corners alone are 67.94% of the whole 0.1 Hz mesh and 73.70% of its PML.
Surface area is a different measure: for L=D, PML occupies 75% of the enclosing
box area, and corners occupy 25%. Mesh counts are not proportional to area.

Each transfinite rectangle is split into two triangles per cell:

- Four corners: each 192 by 192 cells, hence 73,728 triangles each.
- Four side half-strips: each 192 normal by 48 tangential cells, hence 18,432 each.
- Top and bottom: each 192 normal by (21+20) tangential cells, hence 15,744 each.
  The horizontal split at x=1 continues the buried terminal's native voltage path.

For these tangential counts the exact recipe is

`PML triangles(N) = 8 N^2 + 548 N`.

This explains why shrinking the domain does not remove the quadratic corner
cost: four N-by-N corners remain. Shrinking the domain may change strip tangential
counts and physical-region triangles; it does not make the corner count smaller.
Ordinary mesh-size factors cannot override the prescribed transfinite counts.
See the native [Gmsh structured-mesh controls](https://gmsh.info/doc/texinfo/#Structured-grids).

## Actual scale, stretching and cell placement

For the 100 ohm m endpoints:

| Resolved quantity | 0.1 Hz | 1 MHz |
|---|---:|---:|
| Conductive earth skin depth | 15,915.49 m | 5.032921 m |
| Physical half-width D | 381,971.86 m | 120.7901 m |
| PML thickness, each direction | 381,971.86 m | 120.7901 m |
| Entire enclosing box width, 4D | 1,527,887.45 m | 483.1604 m |
| First normal PML interval | 9.377328 m | 0.002965372 m |
| Last normal PML interval | 14,403.58 m | 4.554811 m |
| Side/top strength C | 57,524.80375 | 18.19094018 |
| Bottom strength C | 1.918820910 | 1.913490914 |
| Air phase across physical half-width, k_air D | 0.0008005538 | 2.531573 |
| Earth attenuation scale alpha D | 24.000000 | 23.933334 |

The very large low-frequency dimensions are resolved input values, not a unit
conversion mistake. Physical domain size changes by a factor sqrt(10^7), while
the PML triangle count remains exactly the same. At rho=0.1 and 1 MHz the layout
floor sets D=5 m; the skin-depth rule does not always dominate.

The implemented mesh progression and normalized node locations are

`q = (8N)^(1/(N-1))`, `u_i = (q^i-1)/(q^N-1)`, `i=0,...,N`.

At N=192:

- Adjacent normal widths grow by q=1.0391606108 (3.916%).
- First width is 0.00002454979 L; last width is 0.03770847 L.
- Last/first width is 1,536.
- 48 intervals end at u=0.00333558: 25% of the cells in 0.334% of the thickness.
- 96 end at u=0.02441892: 50% in 2.442% of the thickness.
- 144 end at u=0.15768141: 75% in 15.768% of the thickness.
- 173 end at u=0.48165540: 90.1% in 48.166% of the thickness.

The actual corner nodes agree with this prescription. Actual maximum Cartesian
cell aspect ratios are about 1,536:1. This does not by itself establish numerical
failure: the material tensor is anisotropic too. It does make ordinary
isotropic-mesh quality or physical wavelength-per-cell reasoning incomplete.

The stretched coordinate normal to a layer is

`F(u) = L [u + (1-j) C u^4 / 4]`, for the exp(+j omega t) convention.

The cubic derivative profile, equal real/imaginary strength, N, and the law
q(N) are prescribed. The node distribution does not depend on the frequency,
material wave number or C except through the overall physical length L. Even
changing N also changes the grading ratio, so it is not a pure uniform sampling
refinement. There is presently no public independent PML grading parameter or
separate layer count for side/top/bottom/corners.

This is relevant because the side strength spans 5.75248 to 1.81909 million
across the four-resistivity, standard-radius scan. At rho=100, the point C*u^3=1
moves from u=0.025905 at 0.1 Hz to u=0.380232 at 1 MHz. The current mesh puts
97 and 166 whole intervals, respectively, before that point. Resolving this
transition is not the only requirement, but a fixed distribution cannot be
assumed equally efficient across these scales.

Useful research diagnostics are the real and imaginary parts of
`gamma_m * (F(u_{i+1})-F(u_i))`, where gamma=alpha+j*beta. The CSV evaluates
these exactly for an air-side, earth-side and earth-bottom normal wave at every
frequency, and records log10 amplitude relative to the PML entrance. They are
analytic probes, NOT measured cable fields or error estimates. For example,
the last air-side cell has phase increment 1.641 rad at 0.1 Hz and 1.736 rad at
1 MHz despite N=192. At 1 MHz, cells entered with normalized normal-wave amplitude
at least 1e-3 can still have a phase increment of 1.120 rad. Thus a large N can
coexist with relatively coarse resolution farther into the transformed layer.

## Solver burden associated with this mesh

Across these ten original completed solves:

- 217,709–219,607 mesh nodes.
- 860,179–865,262 solver-reported DOFs. Mesh nodes and coupled-system DOFs
  are distinct counts; the fields use both node and edge bases.
- 11.14–11.22 million actual sparse matrix entries.
- 89.46–89.99 million allocated matrix entries in the native PETSc diagnostics.
  This is a separate storage observation, not a demonstrated assembly bottleneck.
- GetDP reported maximum memory 3,937–4,076 MB per worker; this is the native
  log metric, not a fresh process-tree RSS measurement or total Julia-job memory.
- At 12 quadrature points per triangle, PML alone supplies 4,801,536 integration
  locations per volume term spanning it. This is not the total expression-call
  count: basis-pair loops and multiple coupled terms add work.
- Accumulated worker wall time 1,031.60 s; assembly 858.27 s (83.2%);
  solves 116.48 s (11.3%). Worker time is not elapsed study time.

The original native runs took about 100–106 s per frequency with two source
columns. Later controlled pilot processes were slower on the host; do not
combine those times into a speedup. The same-mesh coefficient pilot paired its
own baseline and candidate: 243.26→159.19 s at 0.1 Hz and 239.55→155.48 s at
1 MHz, with identical saved Z/P entries and derived Y. It leaves this mesh cost
in place.

The three public-route benchmark scans finished during this audit. Baseline
2 workers/4 threads took 1,254.38 s, including 66.95 s reported compilation.
The warmed coefficient candidate took 650.11 s at 2 workers/4 threads and
383.84 s at 4 workers/1 thread. All mesh hashes matched. The first candidate's
Z/Y entries were identical; the second's maximum relative change in any raw G
entry was 2.4644e-6 (0.00024644%), with no R/X/G/B sign changes. The warmed
worker comparison reduces elapsed time by 41.0%; do not attribute the cold
baseline's compilation time to the optimization. These are prototype wrapper
results for one ten-frequency case, not the entire user's parametric grid or
an integrated production change. Evidence:
`.linecablemodels/fem/performance-public-20260929/comparison.csv` and
`component-comparison.csv`.

## Counts if N alone were smaller

These arithmetic projections preserve the measured 0.1 Hz physical mesh
(33,924 triangles), tangential counts, geometry and two triangles per cell.
They do not estimate LU fill, runtime or accuracy, and a real remesh could change
the physical interior even with unchanged boundary targets.

| N | Corner triangles | All PML triangles | Projected total | Reduction in total triangles |
|---:|---:|---:|---:|---:|
| 192 | 294,912 | 400,128 | 434,052 | baseline |
| 128 | 131,072 | 201,216 | 235,140 | 45.8% |
| 96 | 73,728 | 126,336 | 160,260 | 63.1% |
| 64 | 32,768 | 67,840 | 101,764 | 76.6% |
| 48 | 18,432 | 44,736 | 78,660 | 81.9% |

Increasing `mesh_size_factor` affects the small physical fraction and may undo
conductor/interface accuracy. `exterior_mesh_size_factor` can reduce strip
tangential counts but cannot remove the 8N^2 corner term. `domain_skin_depths`
changes PML onset, physical-domain size and the resolved stretch strength; it
is not a direct control over that term. `pml_thickness` and `pml_reflection`
change the transformed scales while preserving prescribed normal counts.

## What the previous evidence establishes, and what it does not

The historical record is `test/manual/calculations/fem_pml_execution.md`:

1. Misaligned nested rectangular rings gave approximately 18% corner-field
   error even after raising N from 64 to 256. Conforming Cartesian patches
   reduced it to 0.4585% at N=64 and 0.4531% at N=256. This establishes the
   importance of mesh alignment, not a requirement for N=192.
2. One older aerial case, rho=1 ohm m at 166,810 Hz, changed its G11 error
   from -2.5344% at N=128 to -0.1872% at N=256. It demonstrates real PML
   discretization sensitivity in that old configuration, not a universal
   lower bound for the current enlarged domain and normalized extraction.
3. The 24-skin-depth/N=192 preset was selected with domain and physical
   mesh controls as well. Its sign checks do not isolate the required N.

These are documented historical measurements, not newly reproduced controls.
Their referenced raw directory `/tmp/fem-pml-execution/` is no longer present
in this environment. The current user's accepted results supersede the old
report's stale statements about overall feature completeness.

It is therefore justified to call the present prescription expensive and
insufficiently optimized. It is not yet justified to declare a particular lower
N safe. The concrete research target is fewer cells through better resolution
of transformed phase, attenuation and coefficient variation, especially in the
tensor-product corners. Preserving corner/interface alignment is independently
motivated. Replacing it with indiscriminate unstructured coarsening is not
supported by the old control.

A bounded later comparison should hold conductor geometry/skin-depth controls,
physical interface targets, D and the PML profile fixed initially; record actual
mesh changes; compare raw individual R/X/G/B entries and signs at low/high and
previously troublesome frequencies, with aerial/buried/mixed coverage. If grading
is redesigned, vary it separately from N. Thickness comparisons under refinement
help separate truncation from discretization. No production self-refinement or
scientific acceptance machinery is implied by these offline comparisons.

## References and precise research questions

- [Gmsh structured grids](https://gmsh.info/doc/texinfo/#Structured-grids):
  native `setTransfiniteCurve(..., N+1, "Progression", q)` and transfinite
  surface interpolation. Implementation owner: `geometry.jl`, approximately
  lines 1600–1658; the actual retained construction is represented by the CSV.
- [Johnson, Notes on Perfectly Matched Layers](https://math.mit.edu/~stevenj/18.369/spring09/pml.pdf),
  sections 3.4–3.5 and 6–7: corners, tensor interpretation, evanescent fields,
  numerical reflection and angular dependence. The notes use the opposite
  phasor sign from this implementation.
- [Oskooi and Johnson, JCP 230 (2011), 2369–2377](https://math.mit.edu/~stevenj/papers/OskooiJo11.pdf):
  comparing different PML thicknesses under resolution refinement distinguishes
  discretization effects from a layer that merely has small reflection at one
  resolution. A test protocol, not a mandated layer count.
- [COMSOL 6.4 PML Implementation](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_definitions.21.137.html):
  its polynomial-stretch guidance starts at at least eight elements across the
  PML, and discusses curvature/scale allocation. This is useful evidence that
  192 is not an intrinsic PML requirement. It is NOT a transferable prescription:
  its profile, element order, formulation, modes and accuracy target differ.

Specific questions to take to the literature: how to distribute Cartesian P1
elements in complex stretched coordinates; whether side, top and bottom require
the same normal resolution; how to retain conforming corner transitions without
an excessive tensor product; how strong shared air/earth side stretching changes
local resolution needs; and how much round-trip attenuation is useful relative
to discretization error when extracting very small signed conductance entries.
