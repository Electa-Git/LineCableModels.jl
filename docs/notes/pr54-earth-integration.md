# PR 54 integration with the accepted earth checkpoint

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
