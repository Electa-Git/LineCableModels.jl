# Detached Gmsh/GetDP model

The entry filename is chosen by the caller's `file_name` argument when exporting;
there is no required basename. This export uses `{{MODEL_NAME}}.pro`, with
matching `{{MODEL_NAME}}.geo` and `{{MODEL_NAME}}_data.pro` files. The commands
below use those actual filenames. Run them from this directory, and keep the
directory together when moving it.

Open the exported entry file in Gmsh:

```sh
gmsh '{{COMMAND_STEM}}.pro'
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

The Numerics panel selects direct MUMPS (the default) or GMRES with LU/MUMPS
preconditioning. In direct mode it includes native MUMPS error analysis
(0 off, 1 full, 2 backward errors), maximum refinement iterations,
backward-error target, and a forward-error comparison budget. Julia defaults
are 2, 2, 1e-12, and 0.01; exported values follow the caller's solver options.
Refinement reuses the LU factors and may stop before the prescribed maximum.
Set analysis and refinement steps to zero to disable both. No regularization,
automatic remeshing or result rejection is introduced.

GetDP output shows native MUMPS estimates after each solve. Compare the sum of
`RINFOG(7)` and `RINFOG(8)` against the backward target. Full analysis also
reports `RINFOG(9)`, an estimated forward error for the scaled field solution;
compare it with the forward budget. These estimates concern the field system. Terminal conductance and admittance
require independent convergence checks. Analysis 2 omits forward-error and
condition estimates, whose output fields contain unused zeros. GetDP 3.5 cannot access these MUMPS values from `.pro`
expressions, so detached execution reports them for inspection; automatic
target-based warnings are specific to managed Julia execution.

GMRES starts from zero, with right LU/MUMPS preconditioning; no preliminary
direct solution is computed. MUMPS refinement and error analysis are disabled
in this mode. The separate GMRES controls prescribe maximum total iterations,
relative residual tolerance and absolute residual tolerance (Julia defaults:
20, 1e-12 and 0). The solver reuses factors for every iteration and subsequent source
right-hand side, and returns its computed result at the iteration limit.
This solver choice applies to the field system; the small terminal-matrix
inversions and reductions remain direct solves.

Native GMRES output reports its convergence reason, iteration count and both
estimated and recomputed true residual norms in PETSc's diagonally scaled
coordinates. The true-residual monitor adds a matrix-vector product and norm
calculations per iteration. Compare its final relative residual with the
relative tolerance, or its true residual norm with the absolute tolerance.
The existing `FEM algebraic residual` line, visible at higher GetDP verbosity,
is in the original coordinates after PETSc restores the system, and is a
separate diagnostic. Detached mode
leaves interpretation of these native diagnostics to the user.

Enable **Run frequency scan** to process every exported frequency in order with
one Run. ONELAB meshes each case before invoking GetDP, using the selected
formulation, basis and mesh settings throughout. The frequency dropdown is hidden
during scanning; switching the checkbox off restores its previous selection.
Every scan starts from the first case. Stop interrupts the scan; another Run
starts it again. The Results panel and displayed fields show the last processed
case; each case's tables and maps remain in its own result directory. The working
`{{MODEL_NAME}}.msh` is replaced at each frequency. Mesh-only scans likewise visit every
case, leaving the last mesh.

Direct native commands also work:

```sh
gmsh '{{COMMAND_STEM}}.geo' -setnumber BuildMesh 1 -0
getdp '{{COMMAND_STEM}}.pro' -msh '{{COMMAND_STEM}}.msh' -solve LineCableModelsFEM
```

Select another exported case with
`-setnumber FrequencyIndex 2` on **both** commands. Solver-only selections include
`-setnumber Physics 1` (Helmholtz, the default),
`-setnumber BasisTerminal 1`, and `-setnumber PlotFieldMaps 0`. An existing
compatible mesh can be selected with GetDP's `-msh` option. Mesh refinement
requires regenerating the mesh. `MeshSizeFactor` is the Julia `mesh_size_factor`;
`ExteriorMeshSizeFactor` is `exterior_mesh_size_factor`. Both start at the
exported values. Changing them updates the local size fields, floating measurement
lines and PML tangential divisions on the next mesh. The exterior factor only
coarsens the remote buffer; it retains the conductor and central bulk targets.
The `Boundary` panel supplies `DomainSizeFactor`, `PmlReflection` and directional
`PmlSide/Top/BottomThicknessFactor`, `Layers` and `Grading`. These are the native
counterparts of Julia's relative thickness, interval and exponent controls.
Dimensions, strengths, roots and physical mesh targets are reevaluated together
for the selected frequency/Γ/material case. The physical conductor CAD stays
fixed. `InterfaceRefinementFactor` matches Julia's interface footprint control;
its native distance fields use the current transverse wave and decay scales.
There is no uniform-refinement or factorization-reuse UI switch. Field-map
output defaults to false.

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
gmsh '{{COMMAND_STEM}}.pro' -setnumber RunFrequencyScan 1 -run
```

## Read and edit the model

- `{{MODEL_NAME}}_data.pro`: named material coefficients with SI units, evaluated
  frequency cases, terminal names/tags, phase connections, source amplitudes,
  normalization and ONELAB controls.
- `{{MODEL_NAME}}.geo`, `geometry/physical.geo`: native entry and serialized
  frequency-independent physical CAD and material and terminal ownership.
- `formulations/parameters.pro`, `geometry.geo`, `mesh.geo`: shared native
  coefficients and extents, exterior and measurement-line topology and mesh prescriptions.
- `formulations/helmholtz.pro`: regions/domains, terminal and
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

Physical cable CAD and evaluated constitutive laws are fixed at export.
Indexed frequency/Γ and material coefficients are editable under Inputs and
Materials. Frequency edits do not evaluate a new Julia material law: edit the
selected case's conductivity, permittivity and permeability when appropriate.
`formulations/parameters.pro` derives the transverse roots, need-based eta per direction,
strengths, domain dimensions and physical mesh targets. The stretch has a cubic
real part and degree-eight imaginary part. Air conductivity is retained in both
roots and equations, and air and soil permeabilities are independent. The sizing-rate
floor at 0.1*b0 is a numerical crash guard. Near-cutoff G remains unqualified. Each medium uses Im(q)≥0 when Im(Γ)<b0_m, otherwise Re(q)≥0.
For b>0, its imaginary-stretch need is clamp((b−a)/b,0,1); for b≤0 it is zero.
Each direction uses the maximum need of its participating media, capped by
the existing negative-b eta bound. Side eta is shared by air and earth. The read-only net exponents Re(q·x̃_end), with
x̃_end=D+L(1+A/4−jηA/4), include the physical domain and PML thickness. A
non-positive net exponent outside exact cutoff rejects an unsupported Γ
prescription. Exact cutoff solves like its neighbours and warns “attenuation
below target; G not qualified”. Sizing-floor activation is informational. G is unqualified when
any net exponent is below (1−1e−9)T, at exact cutoff, or when the resistive sizing ceiling binds. Each frequency writes its native flags, target and exponents to
`raw/jobs/pml-fNNNN.tsv`; managed Julia retains these observations in
`details(result).data.fem.pml_observations` and warns once per run.
The observations also retain `pml_eta` and `effective_pml_layers` (side, top, bottom),
published as read-only derived ONELAB values and in mesh metadata.
Near-cutoff G (approximately |c−1|≤1e-3 for Γ=c·j·k0) is unqualified; exact-
cutoff accuracy and the analytical reference's branch convention require
separate user decisions. All native geometry and mesh files consume the
selected case's current values.
A prescribed Gamma can be non-passive. Native prescriptions whose net exponent
is nonpositive outside cutoff are unsupported. The analytical reference's root
convention is unchanged. Its mean-field single-line-source receiver requires
`abs(kappa_m*r_p) <= about 0.1`; thick receiving conductors near the interface
are outside its declared scope. Exact air cutoff is singular.
Remote/far sizes, including PML tangential sizes, use a medium's wave cap only
when the physical box lies inside its decay footprint at the outer boundary; its fine wave target
inside that footprint is unchanged. Native soil frame and voltage-path counts respect the existing wave target
when its refinement footprint covers the physical box. The footprint width is
computed once in `parameters.pro` and reused by `mesh.geo` as a read-only derived value. The physical half-width is
`D=max(layout_radius,domain_size_factor*L*abs(q0_e)/max(abs(q0_e),abs(q_e)))`.
Here `L=min(1/Re(q0_e),2*pi/abs(q0_e),sqrt(2e5/(omega*mu_earth)))` in metres.
Zero real part means infinite decay length without evaluating division by zero.
`q0_e`, `q_e` are the earth roots at zero and prescribed Gamma.
The rho=1e5 ohm m reference ceiling publishes an active flag and warns
`earth too resistive for FEM domain sizing; results not qualified` when binding.
a larger transverse-root magnitude shortens its numerical scale.

`pml_layers` sets minimum normal interval counts. Native coefficients raise
each count with the fixed constant `PmlPointsPerWavelength=10`:

```text
X_d = L_d*(1+A_d/4-j*eta_d*A_d/4)         # layer only, excluding D
E_m,d = Re(q_m*X_d)
Phi_m,d = |q_m|*|X_d|*min(1,T/E_m,d)   # if E_m,d > 0
Phi_m,d = |q_m|*|X_d|                 # otherwise
PPW_m = 10*clamp(sqrt(0.1*|q_m|/Re(q_m)),1,3) # Re(q_m)<=0 -> 30
N_d = max(N_min,d, ceil(max_m(PPW_m*Phi_m,d)/(2*pi)))
T = -log(pml_reflection)/2
```

Side uses both media, top uses air and bottom uses earth. The decay factor
restricts the interval bound to the stretched-layer portion before attenuation T, so
an already attenuated medium cannot demand excessive counts. Exact cutoff
gives Phi=0. The interval floor remains 48. PPW clamp bounds are read-only
ONELAB values. Effective counts control both
transfinite divisions and `exp(g/N)`. The grading exponent is unchanged.
Named coefficient and constraint edits are read on the next run. Execution
never rewrites source files. ONELAB remembers `DefineConstant` controls: reset
its database to restore file defaults. When relocating a previously opened
project, discard its derived `.db` if it contains stale absolute session paths.

A floating measurement line runs upwards from the interface (overhead) or
outer earth boundary (buried) to an interior anchor in its own terminal metal.
Its native horizontal offset is `min(1e-5,0.01*metal_dimension)` metres. Lines
have independent nodes and do not constrain the field mesh or gauge tree.
GetDP stores `grad v + j*omega*bt` on all media and PML once per excitation,
in memory, then integrates its vertical component with line rule `I2`.
There are no endpoint scalar terms or extra PML stretch factors. Metal stores
no field. Physical line targets are `0.25` times the native exterior local targets. Bottom PML
line intervals are four times the effective count, with the same progression.

An earth-only interface BoundaryLayer field uses `MeshWaveEarth` and growth 1.4.
Its candidate thickness is `min(MeshDecayEarth,0.5*min_buried(abs(CableY)-CableRadius-CableSize))`.
without buried cables it is `MeshDecayEarth`. It is active only when
`MeshWaveEarth < MeshRemoteEarth` and thickness is at least `2*MeshWaveEarth`.
Otherwise effective thickness is zero. Read-only values publish thickness and
the status of a clipped or omitted layer. clipping alone does not warn.
Conductor skin layers remain active. The horizontal background law and
interface footprint factor are unchanged. Vertical box/PML edges grade from
`min(last,bulk,wave)`. Interface seeds are merged natively using domain-scaled
absolute tolerance and 0.1 times local size, with a hard short-segment check.

## Solver controls

The **Numerics** controls include MUMPS ordering and sparse row allocation.
Ordering −1 retains the solver default. 0 requests AMD. Codes 2–7 are native
MUMPS alternatives whose availability depends on the build. Row allocation 0
retains GetDP's default. Positive values set `-petsc_prealloc`. Larger estimates
can reduce reallocations while increasing memory use. These controls prescribe
one execution. No alternate ordering or refinement is tried automatically.

Meshes use binary MSH 4.1. Physical groups, all boundary elements and numerical
coordinates are retained. The bundle still accepts existing ASCII meshes.

## Results

Derived files are under `results/fNNNN-physics-bNNNN/`. `matrices/` contains named
primitive Z and P tables and, for a complete basis, reduced Z, P and Y. Entries
use response rows, excitation columns and the `exp(+j omega t)` convention.
Primitive Z is in ohm/m, P in ohm m, and Y in S/m. Reduction follows the visible
connection map and bundle/Kron/transposition settings. Total output multiplies
Z and Y by line length. P remains explicitly per unit length. No extra `j omega`
factor is used in the FEM admittance inversion.

Source amplitudes also normalize the reported coefficients. Scaling an amplitude
scales the corresponding fields while leaving normalized matrices unchanged.
Changing constraint patterns may change their interpretation. use a diagnostic
basis and inspect the equations for such studies.

`raw/jobs/` holds indexed measurement columns. `maps/` contains native `.pos`
fields, which Gmsh can merge/view and independent artifact readers can inspect.
`completed.txt` is written last and is invalidated when the model is checked or
run again. An incomplete/stopped run may leave partial files. It has no valid
completion marker. The ONELAB Results status distinguishes full and diagnostic
completion. A rerun replaces the derived outputs for that selection.
"Not completed" means the current run has not published a successful result.
It also remains the state after Stop or a solver failure.
