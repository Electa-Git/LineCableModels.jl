# PR 54 integration with the accepted earth checkpoints

## Final integration into `release/v0.2.0`

The target branch fast-forwards to the reconciled checkpoint `90e37b68` without
conflicts. Its earth kernels, spectral integration, uncertainty handling, CIM
fitting and retained numerical regressions match `779a7650` exactly.

At `90e37b68`, both core Julia jobs, deterministic FEM references, quality, Aqua,
clean installation and documentation passed in CI. Gauntlet passed 885 checks
but errored in the new spectral-reporting regression: `(samples,)` supplied a
tuple where integration requires a named tuple. Commit `ece24791` corrects that
argument to `(;samples)`.

Executing the corrected regression exposed a serialization issue: recursively
mapping selection records narrowed their declared named-tuple fields to each
point's options. Reloading a grid containing both `samples=nothing` and integer
budgets then produced incompatible result types. The shared record writer now
preserves declared field types whenever their portable values still fit. It
continues to lower callables and other runtime objects to plain records; hashes
still encode the same values. Numerical result types and the result-space
concreteness check remain unchanged.

The regression checks the concrete reloaded result-space type, live numerical
controls, portable records of the resolved method tags, matrix agreement with
quadrature, report labels and reuse. The focused record and result-space suites
pass 60 checks; the package CLI passes 10. The spectral run passes its numerical,
reload and reuse assertions; its method-tag assertion was corrected to expect
the portable `Base.Val` qualification and is rerun separately. Final spectral
and CI results are reported on PR #54.

## Reconciliation with the CIM checkpoint

Upstream base: `779a7650a79335a8ec6ceccc967c4f832d1c6ca1`, including
`2ae18afc5ed7f81d230db217265193b99cff0cdc` and
`fe09229281b65660418e9937952220ac26c56407`.
The package CLI work was checkpointed first at `2a188bc1`.

The merge into `release/v0.2.0-pscad-ext` requires no textual conflict resolution.
The upstream earth-return implementation, spectral integration, sampling,
analytic bounds, CIM workspace and their production regressions remain unchanged
from `779a7650`. This preserves:

- Adaptive construction by default, with explicit integer `samples` limiting
  construction kernel evaluations per scalar integral. Prototype, pilot and
  verification evaluations retain their separate accounting.
- Rectangular pencils using all samples, reusable factorization across rank
  trials, analytic tail bounds and fit reuse subject to the required certificate.
- Uncertainty-aware verification when a nominal envelope cannot bound derivative
  contributions, promoted envelope arithmetic and explicit numerical failures.

The PR's complete requested/resolved formulation records, shared comparison
reports, Gauntlet lifecycle and package CLI remain in place. Numerical controls
pass through the existing formula options and serialization; no adapter rewrites
or extra control objects are needed. The Gauntlet documentation now explains the
construction budget independently of frequency points and timing repetitions.

`test/gauntlet/spectral_reporting_tests.jl` exercises a zipped Gridspace of trapz
and CIM with adaptive and explicit construction budgets. It computes the same
physical model, saves and reloads the results, checks requested and resolved
controls in comparison reports, compares the matrices with quadrature of the same
equations, verifies saved-result reuse and requires insufficient budgets to fail.
The earlier exploratory earth audit and its `core_material` selector stay retired;
the accepted production earth fixtures remain intact.

### Local validation at checkpoint `90e37b68`

The local run was narrowed on 2026-09-10 at the user's request to conserve
battery. Broad engine, Gauntlet, FEM and documentation runs were stopped before
their final summaries. Their partial logs do not establish suite success.

- Package CLI: 10 checks passed.
- Cairo visual suite: 885 checks passed, including formulation overlays and SVG
  viewport retention.
- Quality suite: 1,193 checks passed; Aqua's persistent-task check failed at its
  precompilation shutdown timeout in both attempts. This remains unresolved; no
  quality check or timeout was relaxed.
- Cold GetDP installation reached the successful continuation after both compiled
  and source-only child processes. The enclosing FEM suite was interrupted later.
- The first engine run also hit the sandbox's read-only FEM run directory. Its
  rerun with write access was interrupted with the broader FEM suite; no numerical
  source change was made for that filesystem failure.
- The focused spectral-reporting regression reached its 60-second limit during
  execution, before a test summary. It remains pending.
- Julia syntax checks for the focused regression and test runner, Git whitespace
  checks and upstream spectral-source/fixture identity checks passed.

The Gauntlet runner now accepts file or test-name fragments. A focused follow-up
uses one process and one BLAS thread:

```sh
OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1 julia --project=gauntlet \
  --startup-file=no --compiled-modules=existing \
  test/gauntlet/runtests.jl spectral_reporting_tests
```

The reconciliation was committed as a checkpoint at the user's request, with
validation pending at that point. Existing CI selects this new regression as part of the
Gauntlet suite; wider numerical and documentation checks can run there. No native
PSCAD solve, publication or push was performed.

## Earlier integration with the earth checkpoint

Integration base: `fe09229281b65660418e9937952220ac26c56407`.
PR checkpoint: `6bbd22b4c33c641fc6254a7edf89529b0d8d85f3`.

The accepted earth equations and production references remain the baseline.
The integration resolves the case relocation, scalar selection records and CI
failures without changing the earth kernels or their numerical algorithms.

## Repairs

- Accept the relocated `build_case(::Val{:two_bare_wires}, p)` implementation.
  Retire the exploratory audit and its obsolete launch instructions. Remove the
  audit-only `core_material` case selector; retain the library PEC material,
  production earth tests and frozen references.
- Retain the complete existing `NamedTuple(formulation)` representation in
  scalar analytical, FEM and PSCAD results. `requested` records parameters and
  hooks; `methods` records resolved selections. Analytical execution records
  remain available through `effective`, `modified`, `numerical` and
  `equivalent_earth`.
- Bound the outer record field types with `NamedTuple`, preserving their actual
  values and callables. Varying a parameter tuple or hook type therefore does
  not create incompatible numerical result containers in a Gridspace. No new
  wrapper, registry or reporting conversion is introduced.
- Require mixed default admittance to compute in both conductor orderings.
  Retain rejection checks with the explicitly restricted Pollaczek selection.
- Install GetDP through the documented `Pkg.Artifacts.ensure_artifact_installed`
  API. Separate-process regressions install and execute GetDP with initially empty
  artifact depots under both compiled and source-only loading.
- Run the existing Aqua test from the independent CI job. Its documented
  dependency exceptions have one owner. Declare Base64 in the test environment,
  as required by the parser fixture's direct import.
- Derive the CIM workspace's two reference field types from the public `Ref`
  constructor. These are the same concrete types as before; this removes private
  `Base.RefValue` access exposed by the merged ownership checks.

The scalar regression changes the earth voltage reference from `:deep` to
`:interface`, verifies the resulting admittance difference, checks distinct
requested/resolved records and labels, round-trips the results, retains an invoked
hook, and compares scalar results with the corresponding zipped Gridspace.

## Validation

Validation uses Julia 1.12.7 with one OpenBLAS thread in the isolated integration
worktree. The following checks have completed:

| Check | Result |
|---|---|
| Scalar/Gridspace records, homogeneous selections and PSCAD boundary/resume tests | 326 passed |
| Frozen earth matrices, spectral integration and uncertainty | 4,234 passed |
| Gauntlet execution, reporting, recovery, archives and vault packages | 829 passed |
| FEM numerical references, constitutive laws and cold GetDP installation | 277 passed under source loading and coverage; one headless GUI check skipped |
| Independent Aqua job in a fresh environment | 11 passed |
| Isolated `Pkg.test` for PSCAD parser and FEM resume | 166 passed |
| CIM fit identities and trapezoid resolution after the type-annotation repair | 22 passed |
| Package quality, equation ownership, semantic economy and Gridspace architecture | 1,194 passed |
| Formula contracts, solver/reduction, FEM records and PSCAD constitutive export | 764 passed |
| Documentation doctests and full site build | Passed |
| Formulation matrix overlays, labels and SVG viewport retention | 89 passed |

The production earth kernels, integral algorithms, uncertainty regression,
unified-earth regression and frozen fixture are unchanged from `fe092292`.
PSCAD checks use retained fixtures and transport processes; no native PSCAD
calculation or benchmark publication was performed during this integration.

The first integrated CI run exposed an additional Julia 1.12.7 artifact-macro
failure with `--compiled-modules=no`: the LazyArtifacts module was present, but
the macro's earlier-world binding check could not see its installer. The initial
cold-cache regression forced compiled loading and missed this case. The direct
public installer replaces that macro path; the regression now covers both modes.
The full FEM suite then passed with CI's `--compiled-modules=no
--code-coverage=@.` settings, and all 1,194 quality checks passed again.

The full prerelease core run subsequently exposed three additional boundaries:

- Both full core runs inferred a missing-valued alternative for the workspace's
  `uses_earth_systems` flag. Its Boolean contract is now explicit; the existing
  workspace inference assertion remains unchanged.
- Packed sectors retain their member-local equivalent-area flattening checks.
  Their overlapping equivalent circles are explicitly rejected by the complete
  earth solver; separated sectors exercise successful default computation.
- The fixed FEM mesh digest depended on Julia's dictionary iteration order.
  Comparing the serialized records under Julia 1.12.7 and 1.13.0-rc4 confirmed
  identical parsed values with different JSON key order. The test now checks
  repeatability and version-sensitive mesh reuse alongside the existing geometry
  and ownership invalidation checks. Production hashing is unchanged; FEM resume
  already records the Julia version.

Both workspace/sector boundary tests pass on Julia 1.12.7 and 1.13.0-rc4
(42 checks on each), as do the revised FEM resume checks (76 on each).
All 1,194 quality checks pass after these repairs.

The stable core job also reached its previous 60-minute limit late in the expanded
suite. CI retains the complete test selection, uses the same single-thread BLAS
settings as local validation, and allows 90 minutes for core and 120 minutes for
the combined coverage job.

The complete CI run at `f5028326` passed every test step: 15,877 core checks
on each Julia version, 15,877 under coverage, 115 PSCAD worker checks, 42 unloaded
extension checks, 890 Cairo/visual checks, and five activation checks for each
of GLMakie and WGLMakie. The combined coverage check reported 16,458 of 17,361
production lines (94.80%), below the unchanged 95% requirement.

A complete local reproduction reported 16,459 of the same 17,361 production
lines. Adding separately measured quality-contract traces covered no additional
lines, so the quality job's arrangement remains unchanged.

Inspection identified unused modal coalescence helpers and history-weight
controls. No registered formula, public option, or production call uses them;
the sole active matching call uses eigenvector overlap. The cleanup removes those
paths and their unused SVD import. The registered Levenberg–Marquardt formula
and the earth equations remain unchanged. Repeated-mode, matched-fallback and
phase/modal round-trip tests exercise the retained modal route.

All 163 selected transformation checks pass after cleanup. Four public modal
calculations, including repeated eigenvalues and forced fallback, produce
exactly unchanged matrices, operators and fallback selections. A further
regression verifies scalar, range, reordered and full-frequency selections:
stored operators follow the same indices, inverse transformation reconstructs
the selected phase matrices, and edits to selected operators preserve the
original data. The expanded modal-tracking file passes 76 checks; all 1,194
quality checks pass. The unchanged coverage check passes at 16,459 of 17,318
production lines (95.04%).

CI now retains `lcov.info` as an artifact after a failed threshold check, so
uncovered code can be inspected without repeating a run solely to recover its
report. Coverage exclusions and the minimum are unchanged.

All nine CI jobs passed at `32125c9b`. The final production measurement is
16,460 of 17,318 lines (95.05%); both core versions pass 15,907 checks. The
additional Codecov changed-line check then reported 89.05% because it combined
the production and Gauntlet inventories. This differs from the established
requirement in `test/README.md` and `test/coverage.jl`, which enforces 95% on
`src/` and `ext/` and publishes Gauntlet coverage separately.

Recalculating changed-line coverage from the retained CI LCOV report and the
GitHub base reproduces Codecov's 3,107 of 3,489 lines (89.05%) exactly. Production
changes cover 2,151 of 2,220 lines (96.89%); Gauntlet changes cover 956 of 1,269
lines (75.33%). Codecov's existing project and changed-line checks now select
the same production paths as the local checker. Both 95% targets remain;
the complete report, including uncovered Gauntlet code, remains published.
