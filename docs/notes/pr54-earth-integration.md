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
