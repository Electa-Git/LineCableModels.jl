# Detached Gmsh/GetDP study

Keep this directory together when moving it. Open the main `.pro` in Gmsh:

```sh
gmsh study.pro
```

Use a native GetDP 3.5 or compatible newer build with Gmsh support. The ONELAB
distribution includes suitable Gmsh/GetDP executables and libraries. Select
GetDP in Gmsh's standard Solver settings if it is not already configured.
There is no additional driver, environment installation or executable setting.

Select a frequency case, formulation, basis and mesh settings, then Run. Basis
zero computes the complete matrices; a positive terminal index computes only
its diagnostic columns. Run action also offers mesh-only execution. Check
parses inputs without solving. Gmsh's Stop controls its native GetDP process.
The GetDP thread control sets GetDP's native `-nt` option; linked numerical
libraries can have their own threading settings.

Enable **Run frequency scan** to process every exported frequency in order with
one Run. ONELAB meshes each case before invoking GetDP, using the selected
formulation, basis and mesh settings throughout. The frequency dropdown is hidden
during scanning; switching the checkbox off restores its previous selection.
Every scan starts from the first case. Stop interrupts the scan; another Run
starts it again. The Results panel and displayed fields show the last processed
case; each case's tables and maps remain in its own result directory. The working
`study.msh` is replaced at each frequency. Mesh-only scans likewise visit every
case, leaving the last mesh.

Direct native commands also work:

```sh
gmsh study.geo -setnumber BuildMesh 1 -0
getdp study.pro -msh study.msh -solve LineCableModelsFEM
```

Replace `study` with the exported entry name. Select another exported case with
`-setnumber FrequencyIndex 2` on **both** commands. Solver-only selections include
`-setnumber Physics 1` (quasi-full-wave, the default),
`-setnumber BasisTerminal 1`, and `-setnumber PlotFieldMaps 0`. An existing
compatible mesh can be selected with GetDP's `-msh` option. Mesh refinement
requires regenerating the mesh. `MeshSizeFactor` is the Julia `mesh_size_factor`;
`ExteriorMeshSizeFactor` is `exterior_mesh_size_factor`. Both start at the
exported values. Changing them updates the local size fields, graded voltage
paths and PML tangential divisions on the next mesh. The exterior factor only
coarsens the remote buffer; it retains the conductor and central bulk targets.
PML normal strips and domain dimensions remain prescribed at export.
`InterfaceRefinementFactor` matches Julia's `interface_refinement_factor`
(default 1; at least 1). It widens the projected cable footprints retained near
the interface. A medium uses these local footprints only when its remote target
already meets the exported frequency/material/Gamma wave-size bound; otherwise
the full interface remains a refinement source. Physical/exterior size edits
reevaluate this condition. Constant fields remain unchanged. This control does
not alter the PML, conductor targets or prescribed voltage-path divisions.
`MeshRefinements` applies uniform refinement. Voltage paths share the field
mesh edges and refine with that mesh.

The Mesh panel also exposes conductor geometry tolerance, elements per skin
depth, normal growth, graded skin depths and wall divisions. These prescribe
local sizes for disks, annuli and convex cable sectors, including individual round screen wires.
The default minimum circle resolution is 96 segments; normal layers start at no more
than one sixth of the material's skin depth. Local characteristic lengths also
cap round-conductor edge lengths, which can add edges on thin foil. `MeshSizeFactor`
scales this local cap and the extension of boundary sizes into passive cable
media; it does not relax the angular bound or the normal skin-depth spacing.
The extension is restricted to the cable's passive surfaces. There is no
automatic convergence test or refinement loop.
Sectors use native triangular transfinite strips around a 0.1-scale core, with
192 angular subdivisions at the default and a matching straight-side cap.
Material/frequency edits recompute the radial counts. The sector strips extend
to the core; the graded-skin-depth limit applies to disk/annular boundary layers.

These two direct commands execute one case. To run ONELAB's complete native
frequency loop without opening a window, use its configured GetDP solver:

```sh
gmsh study.pro -setnumber RunFrequencyScan 1 -run
```

## Read and edit the model

- `study_data.pro`: named material coefficients with SI units, evaluated
  frequency cases, terminal names/tags, phase connections, source amplitudes,
  normalization and ONELAB controls.
- `study.geo`, `geometry/case-*.geo`: native geometry primitives, physical
  memberships, measurement paths/reference points, mesh fields and transfinite
  PML constraints.
- `formulations/quasi-full.pro`: regions/domains, terminal and
  boundary constraints, spaces, field equations, native measurement operations
  and field postoperations.
- `formulations/materials.pro`, `pml.pro`, `integration.pro`: constitutive
  assignments, coordinate stretching and quadrature.
- `formulations/line-parameters.pro`: connection ordering, bundle transforms,
  Schur-complement solves, transposition and the native `P Y = I` solve.
- `formulations/onelab.pro`: selected-case setup, the complete native resolution,
  table output and ONELAB publication.
- `views.geo`: clears this project's derived views before each run or input
  check, including repeated runs with unchanged geometry.

Geometry and evaluated material frequency cases are fixed at export. Select
these cases; changing the frequency alone does not evaluate a new material law.
The case also selects its prescribed complex Γ from `GammaReValues` and
`GammaImValues`, together with its matching PML strength and stretch ray.
Nonzero Γ requires quasi-full-wave. Re-export after changing Γ so the field
operator, mesh and PML coefficients remain consistent.
Conductor sizes are recalculated from the selected frequency and editable
material coefficients when meshing. The exported exterior and PML dimensions
remain fixed; re-export when changed inputs require different domain geometry.
Named coefficient and constraint edits are read on the next run. Execution
never rewrites source files. ONELAB remembers `DefineConstant` controls: reset
its database to restore file defaults. When relocating a previously opened
project, discard its derived `.db` if it contains stale absolute session paths.

A voltage path runs from the air/earth interface for an overhead terminal, or
from the bottom outer boundary for a buried terminal, to its lowest terminal CAD
vertex. Each path interval shares the field mesh edges, including through the
PML. GetDP evaluates the scalar reference and integrates the native `BF_Edge`
trace of its pulled-back transverse vector potential, using line rule `I2` in
`integration.pro`. Metal intervals are excluded. The quasi-full-wave voltage
includes both scalar and inductive contributions. Field-mesh and PML refinement
control discretization error; no separate path sampling is used. PML fields are analytic
continuations, not physical fields outside the interior domain.

## Solver controls

The **Numerics** controls include MUMPS ordering and sparse row allocation.
Ordering −1 retains the solver default; 0 requests AMD. Codes 2–7 are native
MUMPS alternatives whose availability depends on the build. Row allocation 0
retains GetDP's default; positive values set `-petsc_prealloc`. Larger estimates
can reduce reallocations while increasing memory use. These controls prescribe
one execution; no alternate ordering or refinement is tried automatically.

Meshes use binary MSH 4.1. Physical groups, all boundary elements and numerical
coordinates are retained; the bundle still accepts existing ASCII meshes.

## Results

Derived files are under `results/fNNNN-physics-bNNNN/`. `matrices/` contains named
primitive Z and P tables and, for a complete basis, reduced Z, P and Y. Entries
use response rows, excitation columns and the `exp(+j omega t)` convention.
Primitive Z is in ohm/m, P in ohm m, and Y in S/m. Reduction follows the visible
connection map and bundle/Kron/transposition settings. Total output multiplies
Z and Y by line length; P remains explicitly per unit length. No extra `j omega`
factor is used in the FEM admittance inversion.

Source amplitudes also normalize the reported coefficients. Scaling an amplitude
scales the corresponding fields while leaving normalized matrices unchanged.
Changing constraint patterns may change their interpretation; use a diagnostic
basis and inspect the equations for such studies.

`raw/jobs/` holds indexed measurement columns. `maps/` contains native `.pos`
fields, which Gmsh can merge/view and independent artifact readers can inspect.
`completed.txt` is written last and is invalidated when the model is checked or
run again. An incomplete/stopped run may leave partial files; it has no valid
completion marker. The ONELAB Results status distinguishes full and diagnostic
completion. A rerun replaces the derived outputs for that selection.
"Not completed" means the current run has not published a successful result;
it also remains the state after Stop or a solver failure.
