> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# FEM shunt-voltage extraction correction — 2026-09-25

This is a historical contour-averaged extraction study. Its preparation-based solve modes are retired; saved-result plotting remains available. See the [native voltage contract](../../../docs/plans/fem-native-voltage-extraction-plan.md).

## Conductance plot investigation

The two-wire CSVs written at 14:48 on September 25 came from the old backend.
For example, the retained `run-MlAcaK` under the main `LineCableModels` checkout
contains no surface-reference code in `input/getdp/quasi-tem.pro`. Inverting its
1 MHz raw P reproduces the CSV: G11 = 3.2480918586979757e-6 S/m and
G12 = 2.7939915062752223e-6 S/m. Including a script from this investigation
worktree in an existing Julia session does not change the loaded package.
Both manual runners now check the loaded package directory before solving.

The CSV values themselves are nonzero. `ObservedResult` clips G to zero below
the default engineering floor of 1e-12 S/m. At rho = 0.1 ohm m, 278 of 324
entries in the saved 81-frequency sweep fall below that floor. The two-wire
script now exports and displays unclipped, signed G plots before XLSX export;
its retained G/B plots and spreadsheets also use `clip=false`. The package's
general observation policy is unchanged.

`using JSON3` reproduces `Package JSON3 not found in current path` under the
user's default `@v1.12` environment. The manual study no longer imports JSON3
directly; checkpoint decoding uses the FEM backend's existing dependency.
No project or manifest was changed for this investigation.

To plot retained study matrices without solving or reading JSON, run from this
worktree in a fresh Julia process:

```sh
julia --project=test test/manual/calculations/run_fem_voltage_reference.jl --plot /tmp/linecablemodels-manual/fem-voltage-reference-1ftwgM
```

Each case produces `conductance.png` and `conductance.pdf`, showing historical
point, matched deep-earth average, corrected surface reference, and analytical
conductance for all four matrix entries. Panels spanning more than three decades
use a signed logarithmic vertical scale, linear below their smallest decade;
other panels use a linear scale. No values are clipped or replaced by magnitudes.
Four-frequency plots explicitly identify their sample count.

The follow-up also executed one 81-frequency quasi-TEM case: radius 0.0425 m,
positions (0,1) and (1,1) m, rho = 1000 ohm m, unity soil relative permittivity
and permeability, copper at 20 C, line length 1 m, and no reductions. The
retained run is
`/tmp/linecablemodels-manual/conductance-investigation/overhead-rho1000-quasi_tem/runs/run-pGKEpL`.
Its `matrices.csv`, `conductance.png`, and `conductance.pdf` sit in the case
directory two levels above the run. G passes 139/324 entries; B passes 324/324
under the unchanged `1e-12 + 0.01*abs(reference)` criterion. The worst G tolerance
ratio is 48.4103 at entry (1,2), 1 MHz: error +3.963498416726682e-7 S/m against
reference -8.186303916668863e-7 S/m. The maximum absolute G error is
4.0798829137660984e-7 S/m. This is a single dense case, not the full manual study.

At 1 MHz, G11 is +3.2480918586979757e-6 S/m with the historical extraction,
-1.3845190549383485e-6 S/m with the corrected extraction, and
-1.792433120317438e-6 S/m analytically. The reference correction changes the
computed conductance; it does not establish the requested scientific agreement.
Across all 324 complex entries, the historical admittance from this new solve
reproduces the saved old-backend rho = 1000 sweep to a maximum difference of
3.39e-21 S/m. The new run completed 162 source columns with 81 factorizations;
its 81 native GetDP logs contain no warning or error lines.

The saved old-backend resistivity and radius sweeps were also plotted from their
raw CSVs against freshly computed analytical results in
`/tmp/linecablemodels-manual/conductance-investigation/saved-old-backend/`.
All 21 checks passed: nonzero raw G, exact preservation with `clip=false`, and
the expected number of zeros with the default floor, for all seven parameter
points. Headless command-line plotting saves files without opening a desktop
viewer; interactive Julia sessions can display the figures.

## Original implementation and smoke run

Baseline worktree revision: `94e4b058706031e937954802cec3e2921c4dd751`.
The checkout was clean before this work. No analytical formula, material law,
longitudinal wavenumber, PDE, gauge, electrode excitation or matrix reduction was
changed. Julia is 1.12.7. GetDP is the package artifact, version 3.5.0 (PETSc 3.14.4). Its
[manual](https://getdp.info/doc/texinfo/getdp.html#Types-for-PostOperation)
documents local `OnPoint` traces, global `OnRegion` quantities and
`StoreInVariable`; native residual operations are described under
[Resolution](https://getdp.info/doc/texinfo/getdp.html#Types-for-Resolution).

## Implemented measurement

GetDP coordinates are x horizontal, y vertical, z axial; the interface is y=0.
Rows identify receivers and columns identify source terminals. Air receivers use
the local surface projection; buried receivers retain deep earth.

- Quasi-TEM: `Pᵢⱼ = (Vᵢⱼ − ⟨v_surface⟩ᵢⱼ)/qⱼ`, with outward transverse
  source `q` in A/m. Its field is `−grad(v)`.
- Quasi-fw: `Pᵢⱼ = ⟨vᵢⱼ − v_surface + jω∫surface→C bt·dl⟩/Iⱼ`, with
  axial source `I` in A, normalized `v` in V m and `bt` in T m².
  No numerical Gamma is introduced. Buried scalar references and vector paths
  retain their previous definition.
- Raw complex `P` is in ohm m. The existing reductions precede the full matrix
  solve `lu(P) \ I` in `results.jl`, giving `Y` in S/m. Analytical charge-based
  `Pe` corresponds to `jωP`. There is no further `jω` in FEM inversion.

Overhead bare disks use normalized two-point Gauss quadrature on each actual
first-order mesh-contour segment. The scalar reference and vector paths use the
same endpoints and weights. Other geometries and buried receivers retain the
lowest contour node (lowest x breaks ties). No equivalent circle is inferred.
Existing triangle clipping, zero transverse field inside metal and the buried
infinite-shell pullback are preserved. Native algebraic residuals and RHS norms
are logged after each solve; these do not estimate discretization error.

## Changed files and responsibilities

| File | Responsibility |
| --- | --- |
| `ext/LineCableModelsGmshExt/voltage_paths.jl` | Resolve host by terminal/cable ownership; contour quadrature; reference arrays; weighted surface paths; metadata for both physics modes. Scalar mode skips vector geometry preparation. |
| `ext/LineCableModelsGmshExt/workers.jl` | Supply `PathDataPath` for both physics choices through the existing command builder. |
| `ext/LineCableModelsGmshExt/getdp/model.pro` | Include measurement data once and validate flat reference ranges/weights. |
| `ext/LineCableModelsGmshExt/getdp/quasi-tem.pro` | Sample local scalar traces and subtract the receiver reference after source normalization; native residual logging. |
| `ext/LineCableModelsGmshExt/getdp/quasi-full.pro` | Consistent scalar subtraction and weighted vector extraction; preserve `Pscalar` meaning; native residual logging. |
| `src/engine/formulations.jl`, `ext/LineCableModelsGmshExt/results.jl`, `docs/src/fem.md` | Equations, coordinates, units, receiver averaging, inversion and model limitations. |
| `test/extensions/fem_quasi_full.jl` | Geometry/row controls, complex arithmetic, gauge cancellation, direct parsing and native two-column extraction with maps disabled. |
| `test/extensions/fem_resume.jl` | Existing source/input identity rejects old path or solver implementations. |
| `test/manual/calculations/run_fem_voltage_reference.jl` | Same-field historical/matched-deep/surface study; fixed component criteria, signed outputs, mixed/permuted controls and refinement. |
| `test/manual/calculations/run_two_bare_wires_fem.jl` | Correct the manual environment command and describe source/reference conventions. |
| This note | Scope, measured verification and handoff. |

The existing `FEM_ADAPTER_SOURCES` already includes `voltage_paths.jl`, and
`FEM_GETDP_SOURCES` captures all changed `.pro` files. Their hashes participate in
`_fem_input_record` and resume matching. No new cache, revision counter, public
reference selector or production diagnostic matrix was added. By inspection,
production has no analytical comparison, study criterion, model switching,
symmetry enforcement, sign clipping or reference-driven refinement.

## Software execution

Commands ran from the worktree root. This managed host used
`JULIA_DEPOT_PATH=/tmp/lcm-fem-depot:/home/amartins/.julia`, `DISPLAY=`,
`OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.

```sh
JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=test -e 'using Pkg; Pkg.offline(true); Pkg.resolve(); Pkg.instantiate(;allow_autoprecomp=false)'
julia --project=test --compiled-modules=no test/runtests.jl extensions/fem_quasi_full.jl extensions/fem_resume.jl
julia --project=test --compiled-modules=no test/runtests.jl 'native complex reference' 'receiver rows' 'complex receiver measurement'
julia --project=test --compiled-modules=no test/runtests.jl 'receiver rows'
julia --project=test --compiled-modules=no test/runtests.jl --list
git diff --check
```

- Focused suite: **195/195** assertions, including real native solves with maps,
  receiver mapping, gauge/path checks and resume protection.
- Follow-up after the scalar preparation change: **98/98**, including **22/22**
  native extraction checks. These totals overlap; they are not additive.
  Both assembled physics files passed GetDP `-check` and native two-RHS execution.
  Manufactured complex quantities exercised the actual GetDP expressions with
  non-unit source amplitudes, maps off, one invocation and one factorization.
- The final buried-path check passed **88/88**, adding direct checks that buried
  paths reach the deep boundary and their oriented weights span the full ray.
- Ordinary discovery listed **342 items in 135 files**, with zero manual paths;
  no bodies ran in the listing. `test/manual/**` remains excluded by
  `JuliaTestItems.toml`. Julia syntax parsing checked the changed Julia files.
- Exact source-block comparison against HEAD confirmed unchanged constraints,
  function spaces, formulations and excitation macros in both `.pro` files.
- Earlier development attempts exposed a GetDP constant-expression syntax error
  and a Julia range-line-break error; both were corrected before the passing runs.
- Logs: `/tmp/lcm-fem-focused3.log`, `/tmp/lcm-fem-manufactured.log`,
  `/tmp/lcm-fem-buried-paths.log`, `/tmp/lcm-fem-discovery.log`,
  `/tmp/lcm-fem-environment.log`.

## Historical manual scientific study

The former `--smoke` and `--full` solve modes are retired. They depended on
mesh-derived reference/path quadrature, which the native replacement removed.
The retained results and criteria below describe that historical experiment.
The script now accepts only `--plot STUDY_DIRECTORY` to inspect saved CSV data.
Current native extraction checks are in `test/extensions/fem_native_measurements.jl`
and `test/extensions/fem_quasi_full.jl`.

The fixed, provisional criterion is independently
`abs(component error) ≤ 1e-12 + 0.01*abs(reference)` for signed G and B in S/m.
It is not a measured error bound or a production setting. Full-case agreement
requires every entry and frequency of both components. Scientific acceptance
remains the user's decision.

The smoke subset uses copper from the existing material library, 20 C, radius
0.0425 m, x=0 and 1 m, length 1 m, eps_r=mu_r=1, disabled reductions and
frequencies 100 Hz, 10 kHz, 100 kHz and 1 MHz. The overhead case has y=+1 m,
rho=1000 ohm m; the buried control has y=−1 m, rho=0.1 ohm m. Both physics choices
use identical inputs. Effective copper/earth properties and source identities are
captured in each run's existing inputs and the manual study record. The library
copper resolves to sigma=58,001,276.028072625 S/m (rho=1.7241e-8 ohm m),
mu_r=0.999994 and eps_r=1; the specified unity relative properties apply to air
and soil, as in the user's existing script. Copper properties were not altered.

Manual instrumentation changes run-local postprocessing only. Every RHS yields
historical point, matched-averaging deep and surface measurements from the same
field solution. `matrices.csv` retains their P and Y matrices, S, analytical Y,
signed errors, and separate scalar/vector terms; native GetDP/PETSc logs are kept.
`Pscalar` remains a gauge-dependent scalar diagnostic, never the physical Y.

The initial completed subset exposed a manual-output retention bug at process
exit: Julia cleaned up the root temporary directory. The script now uses
`cleanup=false`; the retained rerun is the evidence used below. Its log is
`/tmp/lcm-fem-study-retained.log`.

A separate contour check remeshed the original CAD at mesh-size factors 1, 1/2
and 1/4 (no additional PDE solve). For receiver 1, the normalized second moment
`⟨(x−x_center)^2⟩` has the circular target `r²/2 = 9.03125e-4 m²`:

| Receiver samples | Absolute moment error, m² |
| ---: | ---: |
| 64 | 5.7844312191e-6 |
| 128 | 1.4495979088e-6 |
| 256 | 3.6261787157e-7 |

Both receivers behaved alike, with weight sums within 7e-16 of one. Error
approximately quarters per doubling. Subdivision of an imported mesh leaves its
polygon unchanged and did not reduce this geometric error; the manual full-study
refinement therefore remeshes CAD with a halved mesh-size factor. This check is
about contour quadrature/geometry, not solved-field convergence. The command was
`julia --project=test --compiled-modules=no /tmp/lcm-fem-validation/contour_refinement.jl`,
with the same environment variables as above. The script, meshes and
`/tmp/lcm-fem-contour-remesh.log` are retained locally.



## Retained smoke results

Output root: `/tmp/linecablemodels-manual/fem-voltage-reference-1ftwgM`. All four directories were checked **after Julia exited successfully**. The run contains 16 frequency jobs, 32 source columns and 16 factorizations. Every case has eight completed columns and four factorizations. Maps were disabled.

| Case | G satisfying criterion | B satisfying criterion | Maximum absolute G error, S/m | Maximum absolute B error, S/m |
| --- | ---: | ---: | ---: | ---: |
| overhead-quasi_tem | 4/16 | 16/16 | 4.079882914e-07 | 6.989307532e-07 |
| overhead-quasi_fw | 4/16 | 16/16 | 7.193673695e-07 | 8.247794371e-07 |
| buried-quasi_tem | 14/16 | 16/16 | 0.09777615176 | 0.1040235262 |
| buried-quasi_fw | 14/16 | 16/16 | 0.09604328449 | 0.1040264783 |

These counts cover the four requested frequencies only. **None is a full-case match under the fixed criterion**, because each retains conductance failures. The supplied earlier buried-comparison evidence is not relabeled as this run.

Worst locations below maximize absolute component error divided by its fixed componentwise tolerance. Signed errors and reference values are in S/m; ratio greater than one exceeds the provisional criterion.

| Case | Component | (receiver, source) | Frequency, Hz | Signed error | Signed analytical reference | Tolerance ratio |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| overhead-quasi_tem | G | (1, 2) | 1e+06 | 3.963498417e-07 | -8.186303917e-07 | 48.4103 |
| overhead-quasi_tem | B | (2, 1) | 10000 | -1.6355485e-09 | -1.982942249e-07 | 0.8243932 |
| overhead-quasi_fw | G | (2, 1) | 100000 | 1.261068411e-08 | -6.749506122e-09 | 184.11085 |
| overhead-quasi_fw | B | (2, 2) | 1e+06 | 8.247794371e-07 | 9.496838534e-05 | 0.86847698 |
| buried-quasi_tem | G | (2, 1) | 10000 | -0.03646246525 | -3.32022793 | 1.0981916 |
| buried-quasi_tem | B | (1, 2) | 100000 | 0.0095461408 | 1.050229051 | 0.90895798 |
| buried-quasi_fw | G | (1, 2) | 10000 | -0.03410892895 | -3.32022793 | 1.0273068 |
| buried-quasi_fw | B | (2, 1) | 100000 | 0.01009865012 | 1.050229051 | 0.96156644 |

The separate effects below distinguish changing reference from changing receiver averaging. The conductance differences are in S/m. The final column checks reconstruction of the exported P from its saved scalar/reference/vector components (ohm m); it is an arithmetic diagnostic.

| Case | max abs(Gsurface−Gdeep) | max abs(Gdeep−Ghistorical) | max P reconstruction difference |
| --- | ---: | ---: | ---: |
| overhead-quasi_tem | 4.632610914e-06 | 0 | 0 |
| overhead-quasi_fw | 4.783607897e-06 | 3.794285076e-12 | 0 |
| buried-quasi_tem | 0 | 0 | 0 |
| buried-quasi_fw | 0 | 0 | 0 |

For both buried cases, `Psurface == Pdeep == Phistorical` exactly in the saved tables. No surface samples were used.

The 32 native algebraic residual norms ranged from 8.971043295e-14 to 1.487834287e-09; RHS norms were 1.414213562 to 1.414213562. There were 0 native Warning/Error lines. These residuals are not discretization-error estimates or bounds on conductance accuracy.

Each case retains `matrices.csv`, `summary.json`, generated measurement inputs, source snapshots, `quadrature.csv`, raw columns and native GetDP/PETSc logs. `/tmp/lcm-fem-validation/summarize.py` performed the final retained-data checks.

## Explicitly unexecuted work

The 81-frequency full resistivity/radius study, its physical mixed/permuted comparisons and its solved-field CAD-refinement comparisons were created but not run. Mixed/permuted receiver mapping was exercised in software tests, and contour-only CAD refinement was measured separately. The complete repository test suite, display/UI tests and Documenter build were not run. Scientific acceptance and any tolerance revision remain with the user. No commit or remote publication was requested or performed.
