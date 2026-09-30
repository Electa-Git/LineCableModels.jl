# Native ONELAB export: execution evidence

Date: 2026-09-28. Implements the [execution plan](fem-onelab-export.md).

`export_data(:onelab, problem_or_system, formulation; ...)` now emits a detached
native Gmsh/GetDP project. Julia performs export only. The exported runtime has
no interpreter, package installation, driver, shell handoff, external path
preparation, or external matrix processing.

## Delivered example

The fresh, solved mixed overhead/buried example is
`.linecablemodels/fem/onelab-native-two-wires/study.pro`. It contains two
5 mm-radius wires at `(0,0.1)` and `(0.2,-0.1)` m, 50 Hz and 10 kHz cases,
and eight PML layers. This is an execution example, not a converged engineering
reference. Open it using the native ONELAB distribution:

```sh
/home/amartins/Applications/onelab-Linux64/gmsh \
  .linecablemodels/fem/onelab-native-two-wires/study.pro
```

Inside the exported directory the verified direct commands are:

```sh
gmsh study.geo -setnumber BuildMesh 1 -0
getdp study.pro -msh study.msh -solve LineCableModelsFEM
```

The supplied README documents frequency, formulation, basis, solver, mesh,
field-map and total-length controls. Add `-setnumber FrequencyIndex 2` to both
commands for the second case. `BuildMesh` applies the same uniform refinement
sequence as native ONELAB Run. Check does not solve or generate a mesh.

The earlier `.linecablemodels/fem/onelab-two-bare-wires` is preserved, including
manual edits and derived files. Its native sources describe 42.5 mm-radius wires
at `(0,-1)` and `(1,1)` m, soil conductivity 10 S/m, and ten frequencies from
0.1 Hz to 1 MHz. The current manual sweep source instead places both wires at
`installation_z=-1`. Replacing that directory from the current source would
silently change the study. The new example is deliberately a separate project;
existing exported projects are not automatically rewritten.

## Native frequency scan

The exporter now provisions **Run frequency scan**, disabled by default. When
enabled, ONELAB's native `Loop` metadata advances a hidden index through all
exported frequency cases. Each iteration executes the existing Gmsh mesh action
before GetDP. The independent manual dropdown keeps its selection when scanning
is switched off. The hidden index also lets the checkbox authoritatively enable
and disable looping: ONELAB preserves user loop metadata on visible parameters.
No runtime driver or solver formulation change was introduced.

The delivered example above was refreshed after verifying every owned source
against its recorded hash, and both cases were solved using one native command:

```sh
gmsh study.pro -setnumber RunFrequencyScan 1 -run
```

The maintained frequency-scan test covers native Gmsh/GetDP execution for both
formulations at 50 Hz and 10 kHz, ordered mesh/save/solve events, exact table
agreement with independent manual mesh/solve runs, checkbox metadata transitions
within one ONELAB server, single-case execution, mesh-only scanning and a
one-frequency export. Test processes use distinct native sockets so they can
coexist with an interactive Gmsh session.

Actual GUI checks verified enabling the checkbox after selecting 10 kHz,
both-frequency completion with field maps, restoring the 10 kHz dropdown when
disabled, and a manual Run that left the 50 Hz completion marker unchanged.
Stop during the second frequency left the first complete and the second without
a completion marker. Restart began at the first frequency and completed both.
Results showed the last processed frequency; per-case files remained separate.
The working mesh and displayed maps represent the last case only.

Evidence is under `.linecablemodels/fem/native-export-evidence/frequency-scan/`.
The full exporter run passed 157/157 checks (116 existing and 41 initial scan
checks). The final focused scan run passed 44/44 checks, including the additional
mesh-only scan coverage, and is recorded alongside the native batch and GUI logs. These checks concern scan
orchestration and numerical preservation, not a new physical convergence claim.

## Model ownership

| Source | Visible content |
|---|---|
| `study.geo`, `geometry/case-*.geo` | Native primitives, physical tags, mesh fields, measurement references/curves, grading and refinements |
| `study_data.pro` | Named materials and SI units, frequency tables, terminals, connections, amplitudes and ONELAB controls |
| `formulations/quasi-*.pro` | Domains, constraints, spaces, equations, native voltage integration and fields; shared with the Julia worker |
| `formulations/line-parameters.pro` | Permutation, bundle transforms, Schur-complement solves, transposition and arbitrary-size complex `P Y = I` |
| `formulations/onelab.pro` | Native resolution, selected case, outputs, publication and completion lifecycle |
| `views.geo` | Gmsh view cleanup, invoked on checks and through native GetDP `SendMergeFileRequest` before every solve |

The original worker resolution remains headless. Export uses the same field
system and scan macros; no copy of the physical equations was introduced.
Read-only ONELAB entries and named TSV tables expose primitive Z/P and reduced
Z/P/Y. Single-basis output is explicitly diagnostic and has no full Y table.
A native completion marker is written last, invalidated before input validation,
and absent after errors or Stop. Old raw rows and maps are removed before rerun.
Gmsh's normal GetDP registration supplies the executable; there is no second path
control. Existing global Gmsh solver registrations belong to the user's settings
and are not overwritten by the project.

## Acceptance evidence

| Gate | Verified evidence |
|---|---|
| A01 | Problem/system dispatch, export without a solver/display, native geometry membership, caller models/views/options preserved |
| A02 | Relocated path with spaces, different cwd, empty environment, Linux mount namespace exposing only two native executables and system libraries; no `/home`, repository, Julia or Python executable |
| A03 | Actual native GUI opening/Check, mesh-only Run, full Run, selected basis, Stop, rerun, second-frequency selection and reopened second-frequency full run; configured absolute GetDP path |
| A04 | Shared-mesh primitive and reduced parity at 50 Hz and 10 kHz for both formulations; existing real/imaginary component tolerances unchanged |
| A05 | All reduction flag combinations, reordered/repeated phases, grounds retained without Kron, plus physical grounded and insulated-bundle cases |
| A06 | Production algebra tested with complex nonsymmetric 1-, 2- and 4-terminal matrices and changed coefficients in one process; inverse/reference and identity checks at `1e-13` |
| A07 | Singular algebra fails; zero normalization invalidates prior completion; stopped GUI run has no marker; selected columns have no full Y |
| A08 | Analytic open/closed complex integrals and off-mesh interpolation, explicit overhead/buried references, manufactured excitation-dependent extraction with maps disabled |
| A09 | Fourfold path refinement leaves volume-element coordinates/connectivity identical; one uniform refinement gives fourfold triangle count and doubles line segments; extraction convergence recorded separately below |
| A10 | Native conductivity, source amplitude, constraint sign and connection edits; same-mesh numerical references; mesh-scale remeshing and total-length output; input files remain authoritative |
| A11 | Native GUI fields, compatible list-based POS output, independent mesh/field import, Makie field rendering and geometry/mesh overlay; units and ordering checked |
| A12 | Unchanged full reruns retain 28 views and six lines per two-terminal matrix table; full/diagnostic selections, Stop/restart, source validation failure, maps-disabled cleanup and relocated/reopened sessions |
| A13 | Native-only asset inventory; obsolete owned asset deletion, unrelated-file preservation, unsafe ownership-path rejection; exporter dependency/source/CI audit |
| A14 | Focused export, matrix algebra, measurement, both-formulation, mesh grading, worker, resume and import checks; separate Makie checks |

The primitive comparison remains
`abs(a-b) <= 2e-9*abs(b) + 100eps(Float64)*maximum(abs,b)`, separately for real
and imaginary components. No tolerances were loosened. The delivered 50 Hz
native solve reports `max(abs(P*Y-I)) = 1.0428649832939931e-17`.

The integration refinement experiment on a fixed coarse mixed-wire solution at
10 kHz gave relative extraction differences of 1.122896%, 0.072016%, 0.291954%,
0.079424%, 0.070371%, and 0.018727% for 1, 2, 4, 8, 16 and 32 line subdivisions
relative to the original grading. This is a comparison to an independent
triangle-partition integration reference, not a physical-model error bound.
Convergence is not monotone across element boundaries. See
[fem-native-measurements.md](fem-native-measurements.md) for the conventions and
numerical limits. The production runtime contains no triangle-clipping helper.

## Reproduction and checks

```sh
JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
julia --project=test --startup-file=no test/runtests.jl \
  fem_export fem_native_algebra fem_native_measurements fem_quasi_full \
  fem_mesh_grading fem_resume fem_workers fem_import

JULIA_DEPOT_PATH=/tmp/lcm-onelab-depot:/home/amartins/.julia \
julia --project=test/visual --startup-file=no test/runtests.jl native_makie_fem

julia --project=test test/manual/fem/validate_onelab_export.jl /new/evidence/directory
julia --project=test test/manual/fem/export_onelab_toy.jl /new/toy/directory
```

The first broad focused run passed all numerical/inspection cases but caught two
stale assertions listing the previous source inventory. Both were updated to
include the two new owned `.pro` assets. Final export/resume/algebra verification passed 517/517 checks (116 export,
80 resume, 321 algebra). The other focused scopes passed 198/198 checks:
48 import, 45 grading, 15 native integration/reference, 57 both-formulation and
manufactured extraction, and 33 worker/recovery checks. The Makie test passed
27/27 checks. These are scoped runs, not a claim that the full repository suite
was executed. The full physical edit/reduction runner ends
with `ALL EXTRA VALIDATIONS PASSED`.

Isolation used `bwrap` with an empty root filesystem, read-only `/lib64` and
`/usr/lib64`, only the native Gmsh/GetDP executables bound under `/bin`, a writable
study mounted at `/study with spaces`, fresh `/proc`, `/dev` and `/tmp`, and
`--clearenv --setenv PATH /bin --setenv OMP_NUM_THREADS 1 --chdir /elsewhere`.
Both meshing and solving completed there. No launcher/interpreter executable or
repository path was mounted. A PATH-only test was not used as independence proof.
GUI automation tools used during development are outside the exported runtime.

Tested: Linux, Julia 1.12.7 for export/tests, Gmsh library 4.15.2, native GUI
Gmsh 4.14.0-git-67db5bd93, GetDP 3.5.0 and 3.6.0-git-1cf7fa06.
Windows/macOS were not tested. Physical convergence of arbitrary user models is
outside these execution/preservation checks.

Local evidence is under `.linecablemodels/fem/native-export-evidence/`; it
includes solver/test logs, native GUI capture, separately rendered mesh/field
images and the task-specific changed-file manifest. No commit was created.
