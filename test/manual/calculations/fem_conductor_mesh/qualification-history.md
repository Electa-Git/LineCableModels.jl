> Historical qualification procedure, archived on 2026-09-29. Executable
> prototypes and Python tooling referenced below have been retired. These
> descriptions and results are retained for interpretation of saved evidence;
> use [README.md](README.md) for current Julia entry points. The original
> source files are in the local [cleanup archive](../fem_development_cleanup.md).

# Conductor mesh qualification experiments

These are manual qualification scripts, not FEM engine code. No production
defaults or scientific acceptance rules are changed by these experiments.

The implemented conductor mesh and current delivery status are described in
[delivery-assessment.md](delivery-assessment.md). The original bare-wire and
three-cable selections completed. `run_interface_delivery.jl` is the bounded
public-route check after the qualified foil-edge/insulation correction: eight
1 MHz cases, both formulations, existing fixture definitions. It retains source
identity and skips completed cases. All eight cases and their saved-matrix
assessment are complete; see `interface-delivery/assessment.csv` under the
qualification output root. This closes the conductor mesh delivery scope, while
the separately recorded bare-wire G sign discrepancies remain open. No repeated
historical grids are needed.

[`fixtures.jl`](fixtures.jl) contains passive cable constructors; inclusion starts
no computation. [Dimensions and coverage](fixtures.md) accompany the code.
The older full-line/path-preparation experiments below used historical exports;
they are not entrypoints for the normalized native backend. Their local conductor
qualification results remain evidence for the current construction.

Current environment on the investigation host:

```text
python3 with PYTHONPATH=.linecablemodels/fem/conductor-mesh-qualification/python-env/lib/python3.11/site-packages
Gmsh 4.15.2
NumPy 2.4.6
```

The original `/tmp/lcm-onelab-export-env` used in the historical commands below
was lost during a reboot. Gmsh 4.15.2 is now retained in the workspace prefix
above; the ordinary Python environment supplies NumPy and mpmath. Completed
results are preserved and should not be regenerated to restore this environment.

`summarize_screen_air.py --root ROOT --output REPORT` assesses the completed
corrected screen batch, including its final recovery and the six reused controls,
without meshing or solving. It writes raw tables, refinement/reference comparisons
and standalone PNG/SVG plots; it does not claim warmed whole-engine performance.

`qualify_sector_proximity.py` places two reflected copies of captured Julia-owned
sector contours across a 0.5 mm gap. It uses the same +1/-1 A terminal-current
problem as the screen controls, with an A_z=0 circular boundary at 25 mm radius.
It preserves the experimental transfinite conductor partitions and extends mesh
sizes into air. `--levels` explicitly selects independent normal/tangential
refinements or the two locally isotropic Delaunay references. Artificial partition
edges are checked to be absent from the combined terminal contour. `--maps`
writes native current-density maps. This is a coupled local proximity control,
not the full earth/PML Z/Y formulation or a replacement for material-contact tests.

`qualify_gmsh.py` constructs small built-in GEO fixtures, measures conductor area,
actual first-layer heights on both material sides, shared-edge coverage and
signed triangle Jacobians, and saves meshes plus native GEO files. Its shape
dimensions and refinement parameters are controlled experimental inputs. It
does not solve Maxwell's equations or establish impedance convergence.

For example, reproduce the independent-layer and annulus checks in a fresh
output directory:

```sh
/tmp/lcm-onelab-export-env/bin/python -u \
  test/manual/calculations/fem_conductor_mesh/qualify_gmsh.py \
  --case pair_layers --case annulus_0.2 --case annulus_2 --case annulus_20 \
  --roundtrip --output /tmp/conductor-mesh-circles
```

`--roundtrip` compares API meshing, the native writer's unmodified output, and a
native file with explicit `BoundaryLayer Field = ...;` activation declarations.
The raw-writer route is a deliberate negative control. A failed geometric recipe
is recorded and not retried automatically; independent selected fixtures still run.

`--fit-extent` adjusts the prescribed first size downwards so an integer number
of geometric layers reaches the requested depth. The raw normal-column CSV files
record measured spacing; unstructured controls may have no aligned column.

`qualify_pml.py` reuses the retained mixed-pair geometry and exterior constraints,
then compares geometry-only refinement against two independent conductor layers.
It verifies unchanged physical group membership, PML element counts and PML node
coordinates. The skin depths in this mesh-only control are prescribed probe
values, not the material skin depths at the retained file's 0.1 Hz frequency.

```sh
/tmp/lcm-onelab-export-env/bin/python -u \
  test/manual/calculations/fem_conductor_mesh/qualify_pml.py \
  --source .linecablemodels/fem/conductor-geometry-evidence/case.geo \
  --output .linecablemodels/fem/conductor-mesh-qualification/q1-pml \
  > .linecablemodels/fem/conductor-mesh-qualification/q1-pml.log 2>&1
```

Use a new output directory for a new experiment. Existing evidence is preserved.
See [the execution notes](../fem_conductor_mesh.md) for measured results and
remaining qualification work. All scripts are excluded from automated test
discovery by the existing `test/manual/**` exclusion.

## Isolated electrical qualification

`qualify_impedance.py` uses `internal_impedance.pro` for a small conductor-only
GetDP diffusion problem. With exp(+j omega t), `q=sqrt(j omega mu sigma)` and
`-laplacian(Ez)+q^2 Ez=0`. Set Ez=1 at the outside; for a tube set its inner
normal derivative to zero. Compute `Zint=1/integral(sigma Ez dS)` and retain
integrated Joule loss as an independent normalization check. No exterior is
solved, so these results qualify internal impedance only.

The exact circular solution uses I0(qr); the tube uses
`I0(qr) K1(qb) + K0(qr) I1(qb)`, which has zero derivative at inner radius b.
Its area integral gives the current. These follow the
[NIST Bessel derivative identities](https://dlmf.nist.gov/10.29).
At DC the reference is rho divided by exact area. mpmath avoids exponential
overflow; the package's analytical implementation is not used.

The initial levels are candidate (96 segments, delta/3, growth 1.25), tangential
refinement (192 segments), normal refinement (delta/6, growth sqrt(1.25)), and
both refinements combined. Layer extent is fitted as above. Each case retains
its mesh, GEO declarations, GetDP input/log/solution, measured integrals and cost.
CSV target flags are harness evidence; failed accuracy targets do not trigger
automatic refinement or imply an engine policy. Genuine execution errors stop
the batch. GetDP costs include process startup; they are not Julia or warmed
production performance measurements.

Run five solid radii and two tubular sections across DC and the eight specified
frequencies (252 cases), serially:

```sh
env PYTHONPATH=/tmp/lcm-onelab-export-env/lib/python3.11/site-packages \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 -u test/manual/calculations/fem_conductor_mesh/qualify_impedance.py \
  --getdp /home/amartins/.julia/artifacts/51e049a32eeb2ffebcc8a5945bb00f70470a3d29/getdp-3.5.0-Linux64/bin/getdp \
  --tubes-mm 42.5:34,42.5:41.65 \
  --output .linecablemodels/fem/conductor-mesh-qualification/q2-round-spectrum \
  > .linecablemodels/fem/conductor-mesh-qualification/q2-round-spectrum.log 2>&1
```

The main Python environment supplies mpmath 1.3.0; the explicit PYTHONPATH supplies
the already installed Gmsh/NumPy environment. There are no package installations.
Monitor the log or `measurements.csv`; the agent ends its turn after confirming
startup and resumes only when the user asks. The full line formulation and cost
qualification are subsequent work, not covered by this isolated reference.

The completed spectrum and its errors can be summarized without another solve:

```sh
env MPLCONFIGDIR=/tmp/lcm-mesh-qual-mpl /usr/bin/python3 \
  test/manual/calculations/fem_conductor_mesh/summarize_impedance.py \
  .linecablemodels/fem/conductor-mesh-qualification/q2-round-spectrum/measurements.csv \
  --output .linecablemodels/fem/conductor-mesh-qualification/q2-round-summary
```

The system Python supplies Matplotlib for standalone report figures. The optional
`normal_fine` level adds delta/12 normal spacing with growth `1.25^0.25`, retaining
96 boundary segments; `combined_fine` also uses 192 segments. These levels are
explicit requested experiments and never automatic retries.

## Sector qualification

`export_shape_controls.jl OUTPUT` runs under the repository's Julia project and
exports five sharp/offset/rounded sector fixtures using the existing shape owner.
It writes ordinary native GEO plus a CSV of exact area and centroid, without
meshing or solving. `qualify_owned_shapes.py --source OUTPUT --output NEW_OUTPUT`
then tests native transfinite strips and a small core within each convex shape.
The physical conductor contour is retained and artificial cuts stay internal.

`--getdp PATH` additionally runs the conductor-only diffusion control. Here Ez=1
on the complete exterior conductor contour. The resulting response is a
prescribed-boundary diffusion benchmark; it is not the self impedance of the
complete line in an exterior magnetic field. Reference it to DC rho/area and
independently refined meshes. Full line/proximity tests remain required.

`--ratios` gives delta/back-radius (`inf` means DC); `--shapes` selects fixtures.
The levels `normal`, `normal_fine`, `combined`, `combined_fine` independently
vary normal and tangential resolution using the round-study parameters above.
The explicit `isotropic` control uses size <= delta/12 inside the conductor and
a finer boundary, with coarse exterior sizing. It is an independent reference
for measuring electrical error and the cost of the anisotropic construction.
Native boundary-layer failures remain recorded and do not select another recipe.

The completed serial invocation is retained in the ignored evidence file
`.linecablemodels/fem/conductor-mesh-qualification/q2-followup.sh`; its output is
`q2-followup.log`. It completes 63 round refinement controls, 72 new sector
controls and five isotropic references. The eight completed sector pilots are
reused as evidence instead of rerun. No supervisor or notification service is used.

## Current native three-wire integration and cost check

`export_native_wires.jl OUTPUT` creates a fresh native export of the existing
mixed three-wire fixture and a Gamma=0 analytical comparison. It performs no FEM
mesh or solve. `qualify_native_wires.py --input INPUT --output OUTPUT --getdp PATH`
runs the prescribed baseline, normal and combined-fine conductor meshes at
0.1 Hz, 10 kHz and 1 MHz, with both formulations. It uses Gmsh/GetDP directly;
there is no Python voltage preparation or detached solver implementation.

The batch has 18 cases / 54 source columns. Three excitations share one field
factorization at each frequency/formulation/mesh. Fixed native AMD ordering and
one solver thread apply to every compared level. Mesh and solver costs are
recorded separately; these are not warmed whole-Julia-engine timings. The
current public export supplies all equations, field-edge measurements and matrix
normalization. Only experimental conductor constraints differ between levels.

Each mesh retains its constraints, region measurements and PML coordinates.
Each solve retains native matrix tables, actual GetDP output and `/usr/bin/time`
wall/RSS data. Running the same command again skips completed meshes/solves;
the captured harness must remain unchanged for that resumption. The active
invocation is retained in the ignored evidence directory:

```sh
bash .linecablemodels/fem/conductor-mesh-qualification/run-native-three-wire.sh
tail -f .linecablemodels/fem/conductor-mesh-qualification/native-three-wire.log
```

The launcher prevents a duplicate batch with a nonblocking file lock. It is an
ordinary serial foreground process, with no supervisor or timer. The log streams
actual mesh/solver output. The agent returns after confirming startup.

## Historical complete two-wire formulations

This section records the completed campaign's commands. `qualify_line_pair.py`
and `run_line_followup.sh` require the old bundle's `driver.py` and prepared-path
contract. Fresh exports have no such driver. Do not use these launchers with a
current export, restore the deleted preparation code, or interpret their stored
Z/P/Y as current-route validation. The result collectors remain useful for
reading retained evidence. Current full-line checks must use public `compute`
or fresh native Gmsh/GetDP exports and the conforming field-edge voltage paths.

`export_line_controls.jl OUTPUT` exports the twelve placement/resistivity fixtures
through the existing public API and computes matched primitive analytical
references. It does not mesh or solve FEM. Run it with the repository project:

```sh
julia --compiled-modules=existing --project=. \
  test/manual/calculations/fem_conductor_mesh/export_line_controls.jl OUTPUT
```

`qualify_line_pair.py` applied the experimental round-wire constraints to these
captured exports, then called the former exported detached driver for solves.
It does not reimplement the equations or measurement paths. For the mixed-pair
pilot, use a new result directory:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /tmp/lcm-onelab-export-env/bin/python -u \
  test/manual/calculations/fem_conductor_mesh/qualify_line_pair.py \
  --entry OUTPUT/air_earth-rho100.0/study.pro \
  --indices 6 --levels normal,baseline,combined_fine \
  --getdp /path/to/getdp --output RESULTS > RESULTS.log 2>&1
```

Indices 1–8 correspond to 0.1 Hz through 1 MHz in decades. Both formulations
run serially by default (`--physics 1,0`). `--maps` requests the existing native
field outputs. Raw Z, P and Y remain in each detached run; `runs.csv` records
their paths and costs. Each mesh has `mesh.csv`, native constraints and a direct
`case.geo` entrypoint. The native copper sizing expression responds to changes
in the exported sigma/mu/frequency arrays; this bounded circular-fixture harness
is not a generic production implementation. PML coordinates must agree across
selected refinement levels.

The wrapper records actual GetDP process wall time and peak RSS; mesh time is
measured separately. These are not first-use/warmed Julia cost measurements.

`summarize_spectrum.py --root CAMPAIGN --output SUMMARY` collects the completed
initial 192 cases and retained pilot results, including raw signs, crossing
brackets, resistance errors, costs and algebraic residuals. It performs no solve.
`plot_spectrum.jl SUMMARY/matrices.csv NEW_PLOT_DIRECTORY` renders R/G/B through
the public observation/plot API with `clip=false`. GLMakie is the default;
`FEM_PLOT_BACKEND=cairo` explicitly selects static report generation.

`export_line_followup.jl NEW_INPUT` exports the fixed intermediate-frequency
cases selected from the initial spectrum. `run_line_followup.sh ORIGINAL_INPUT
NEW_INPUT NEW_RESULTS GETDP` performs the 108 fixed crossing/refinement/baseline
controls recorded in the execution notes. It runs one solver at a time and
appends actual solver output to `NEW_RESULTS/solver.log`. The user retained the
then-existing path preparation; that execution decision is historical and has
been superseded by the current native integration contract.

`summarize_line_followup.py --root CAMPAIGN --output SUMMARY` assesses the
completed 108 follow-up cases, including recovered subdirectories, against
retained original results. `plot_matrices.csv` can be passed to the plotting
script; each case keeps its actual frequency grid.

## Coupled-current screen controls

`qualify_screens.py` uses the test-only `coupled_conductors.pro` to resolve a
core and explicit 3 mm screen strands in air, with +1/-1 A terminal constraints.
It records integrated losses, terminal voltage/current, residuals, mesh geometry
and cost. This local magnetoquasistatic benchmark does not solve the full layered
Y problem. `--strands 0` selects the independent single-wire reference control.

`run_screens.sh NEW_OUTPUT GETDP` runs the 56 fixed serial follow-up controls
described in the evidence notes, reusing the two completed setup pilots. It
requires the documented Python/Gmsh/mpmath environments and writes native solver
output to `NEW_OUTPUT/solver.log`. The isotropic references use conductor-local
native distance sizing. All decisions about accuracy and feasibility remain in
the qualification harness.

The screen schedule explicitly uses `--mumps-ordering 0` (AMD). The default
PORD selection in this GetDP build spent over 35 minutes in symbolic ordering
on one isotropic reference; the byte-identical mesh/equations completed in
47.39 seconds with AMD. The harness now emits native MUMPS phase messages.
This is a fixed test setting, not a production default or automatic fallback.
The recovery schedule and same-mesh preservation checks are recorded in the
execution notes and `q2-screen-ordering` artifacts.
The high-resolution isotropic reference groups also prescribe
`--petsc-prealloc 1024`: their air triangle fans exceeded GetDP's default 100
entries per ordinary row and caused expensive sparse-storage reallocations.
This consumes more memory and changes neither mesh nor equations. The
same-mesh preservation check is in `q2-screen-preallocation`; this remains a
qualification-only setting. GetDP reports progress at 1-percent increments
with operation timings (`-p 1 -cpu`).
`summarize_screens.py --root QUALIFICATION_ROOT --output SUMMARY` reads the 56
scheduled controls and two pilots, writes numerical refinement/reference tables,
and renders resistance and cost plots without solving again.

The completed screen batch did not pass the refinement target. Its air mesh was
too abrupt to isolate conductor accuracy. `--air-growth VALUE` applies a native
Gmsh `Extend` field, restricted to air with `IncludeBoundary=0`; it retains the
same remote size cap and uses the actual contour-edge lengths. The first pilot
preserved conductor coordinates/connectivity exactly after accounting for node
renumbering, while changing loop R by 20.1345%. See its separate
`conductor-preservation.csv`, which compares saved meshes; its original per-row
experimental hash included node numbers and is not the preservation criterion.
`run_screen_air.sh NEW_OUTPUT GETDP` runs 17 serial follow-ups on the 24-strand
fixture, reusing this pilot. It separates normal/tangential conductor refinement,
three air growth levels, isotropic references, and low/transition frequencies.
This investigates the qualification fixture; it does not change production air
meshing or establish a cause for unrelated full-line Y discrepancies.

The 17 follow-ups completed. See `q2-screen-air-summary`: the 24-strand layered
normal/combined-fine R changes are below .25% at 50 Hz, 10 kHz and 1 MHz with air
growth .125. Exact conductor preservation across air growth holds for the normal
and isotropic levels. Fine layered reconstructions have some interior remeshing
differences; their comparisons are not strictly air-only.
The finer isotropic reference is unqualified: its Frontal-Delaunay core mesh
has fewer triangles than the coarser reference and loses skin resolution.
The mesh-only check in `q2-screen-reference-mesh` prescribes native Delaunay
(algorithm 5), which restores the requested core grading. The explicit
`--isotropic-algorithm` selects that algorithm on reference metal surfaces only.
`run_screen_air_spectrum.sh NEW_OUTPUT GETDP` runs 42 serial controls: 30 missing
graded layout/count/frequency cases and both Delaunay reference levels for all
six layout/count combinations. It reuses the six completed 24-strand fixed-gap
graded controls. No scientific acceptance or recovery policy enters production.

## Historical full-line logging and preparation controls

The following launchers and timing probes refer to the retired detached-driver
campaign. They remain useful for interpreting its retained files, not for running
the normalized MVP.

Qualification stopped on execution or mesh-integrity failures, with no automatic
retry or refinement. The outer batch log reports meshing and stage transitions;
it does not contain GetDP's intermediate output. New invocations also append
line-buffered GetDP output to `RESULTS/solver.log`, while preserving each native
run's `getdp.log`. Follow it in a terminal with `tail -f RESULTS/solver.log`;
opening a file link in the chat does not provide terminal streaming. During
voltage-path preparation there is no GetDP output yet, and the last line says
`PREPARE`. Resume the agent when the batch ends.

The first pilot predated this logging adjustment; its native
`work/run-*/getdp.log` files retain solver progress from that run.

`qualify_path_bounds.py --entry ENTRY --pilot RESULTS --output NEW_OUTPUT`
measures first/repeated path preparation on the retained pilot meshes and requires
byte-for-byte agreement with their original `paths.pro` files. It performs no
GetDP solve. The explicit `qualify_line_pair.py --reuse-path-bounds` option uses
this test-only experiment; it does not alter the exported driver or its clipping
arithmetic. `--progress-log FILE` appends solver output from successive serial
invocations to one file, so `tail -f FILE` continues across cases.

`summarize_line_pair.py --input RESULTS --reference ANALYTICAL_CSV --output SUMMARY`
tabulates raw Z/P/Y, analytical references and candidate/fine changes without
another solve. Relative differences do not classify a near-zero sign as resolved.

`run_line_spectrum.sh INPUT OUTPUT GETDP` starts the complete serial spectrum
from the repository root, with a fresh output directory. It reuses the two
normal-level mixed/100 ohm m/10 kHz pilot solves and computes the other 190.
Both formulations, three placements, four resistivities and eight frequencies
are covered. The mixed 100 ohm m case additionally writes maps at 1 MHz.
Follow native solver output with `tail -f OUTPUT/solver.log`. This is an ordinary
foreground batch, with no supervisor, polling agent or parallel solver workers.

`native_line_integral.geo` / `.pro` are an independent 232-DOF manufactured
complex edge-field example. Mesh the GEO file, then run GetDP's `Check` resolution
with that mesh. Native line-group integrals should be `-1.65-3.3j` on the internal
vertical path and `2.8+5.6j` around the exterior. The example proves integration
over conforming line groups; it is not a replacement for the cable voltage-path
convention. It uses no custom path preparation.
