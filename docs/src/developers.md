# Commons invariants

LineCableModels has three global calculation roles: a complete problem and a formulation that selects how to calculate it, followed by a completed result. Package
modules may add concrete types below these roots, but they do not add parallel
calculation entry points.

The type trees below are generated from the loaded package during the
documentation build.

```@example grammar_type_trees
using LineCableModels

println("Problem definitions")
print(Main.DocumentationTrees.type_tree(AbstractProblemDefinition))
println("\nFormulations")
print(Main.DocumentationTrees.type_tree(AbstractFormulation))
println("\nResults")
print(Main.DocumentationTrees.type_tree(AbstractProblemResult))
```

`Combinatorial`, `LinearError`, and `MonteCarlo` occupy the same formulation
tree. UQ uses the same calculation supertypes. Its completed
collections are subtypes of `AbstractUncertaintyResult`, while deterministic finite
collections are subtypes of `AbstractParametricResult`.

`LineParamsDomain` is independent of this calculation grammar. It tags the
physical coordinate system of a completed line-parameter matrix:

```@example grammar_type_trees
print(Main.DocumentationTrees.type_tree(LineCableModels.LineParamsDomain))
```

## Fixed actions

One declarative action defines a fixed sequence selected by an abstract definition
type:

| Action | Definition root | Fixed sequence |
|:--|:--|:--|
| `ReportBuilder.report` | `AbstractReportDefinition` | `select → tabulate → illustrate → encode → write → ReportArtifact` |

```@example grammar_type_trees
using LineCableModels.ReportBuilder: AbstractReportDefinition

println("Report definitions")
print(Main.DocumentationTrees.type_tree(AbstractReportDefinition))
```

Each action has one method at its declared abstract root. Concrete definitions
implement stage methods. They do not specialize the public action itself.
Required stages are declared with RequiredInterfaces where the type family
admits that form. Optional stages have an explicit no-op at the abstract root.
Plotting is deliberately not such an action: the optional Makie extension
constructs native figures directly from published observations.

Observation uses a different structure. Result owners add `observe` methods.
Commons owns one `observables(source, requests::Tuple; ...)` publication
method. Standalone arrays and external result types do not need a shared
abstract observation type.

## Architecture checks

The architecture that the quality, core and integration tests protect requires:

- the action and its abstract root to belong to the same module.
- one public action method, with no more-specific definition methods.
- every fixed stage to remain visible in the declared action.
- every package-owned concrete definition to implement its required stages.
- every package-owned definition to remain below its declared root.
- one Commons-owned observation publication method with positional requests.
- wide scientific tables whose quantity and unit metadata remain attached to
  their columns.
- no calls to another package module's private functions or types, including
  through aliases, and no private wrappers whose methods only forward unchanged
  arguments to the same package-owned function. Owner-local numerical kernels
  and type-dispatch branches remain valid.

A new definition must satisfy these standards. Native interface and import checks,
tests through actual composed consumers and the structural guards below provide the
protection. Unchanged forwarding wrappers remain an advisory ReLint finding.

The quality selection also runs advisory source diagnostics. Fatou applies the
explicit rules in `fatou.toml`. ReLint identifies runtime-evaluation spellings,
empty or constant-result catches, and unchanged private forwarding methods.
These scans read all Julia files under `src/` and `ext/` without executing them.
They report candidates for review, including typed forwarding methods whose
dispatch purpose has not been determined. Macro expansion, arbitrary binding
ownership, and global architectural redundancy are outside their analysis.
Source findings do not fail the job. Tool setup errors, incomplete scans, and
failed integration controls do. The scans leave source files unchanged and
do not contribute executed-code coverage. See the
[test instructions](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/README.md#advisory-source-diagnostics)
for pinned tool installation, commands, and the retained CI report.

Integration tests also count actual blueprint lowering calls: one per selected
design point, shared across its formulation alternatives and frequency sweep.
They check material evaluation before local calculations, paired exterior
calculation before assembly and reduction at each frequency and explicit
phase-to-modal result transport without another geometry lowering. The visual
suite applies the ownership checks to the loaded Makie extensions and verifies
that material colors consume `Material` or `EarthLayer` objects directly.

InputValidation tests check owner dispatch, required interfaces, unchanged valid
inputs and rejection of damaged inputs through validation and computation. The
standards require checks to remain with the defining validator. Inspection of every method body requires a separate review beyond representative behavioral tests.

The native FEM environment uses pinned GetDP 3.5.0. Keep tests of its actual
execution, extraction, terminal identity, material and option transport and resume
behavior. Physical cross-backend accuracy and domain-convergence acceptance
belong to explicit research work under the testing requirements below. See the
[test commands](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/README.md) for the implementation checks and their
execution environments.

### Structural guards

`test/quality/architecture.jl` enforces one rule: each concept has one owner, one
name and one implementation, and each name has one meaning. The guards load the
package with the extensions that `explicit_imports.jl` loads. They enumerate every
package method and parse every Julia file under `src/` and `ext/`.

- Ownership (A1). A core method that extends a function owned by another package
  module mentions a type owned by its own module or a descendant. For a constructor,
  the function owner is the owner of the constructed type. For methods of the root
  module, functions and types of its submodules count as owned by another module.
  Extension methods are exempt.
- Placement (A2). No core method is defined in a file under `ext/`. The home of a
  module is the directory that holds its `<ModuleName>.jl` file. The nearest home
  around a method's file belongs to the method's module or one of its ancestors.
- Direction (A3). The submodules have the order Units, Commons, TextDisplay,
  InputValidation, PlotBuilder, Materials, Earth, DataModel, Engine, ModalAnalysis,
  ParametricBuilder, UQ, ReportBuilder, ImportExport, PSCAD. The lowered code of a
  submodule method references earlier submodules, the ancestors and descendants of
  its own module, and no later submodule. Each top-level submodule has a position in
  the order. The root module and the extensions are exempt.
- Names (A4). Functions and types owned by different package modules have distinct
  names. The per-family `Formula` types are exempt. A public package name that Base,
  DataFrames, LinearAlgebra, Statistics, Random, Dates or Logging also exports
  refers to the same object as the exported name.
- `validate` (A5). Each `validate` definition names its first positional argument.
  Each path through its body returns that argument or calls `throw`, `rethrow` or
  `error`. The guard reads source code. A path that it cannot classify counts as a
  violation.
- Reserved verbs (A6). No function name starts with `validate_`, `check_`,
  `require_`, `assert_`, `verify_` or `ensure_`, after any leading `_` characters.
  `validate` is the input-check verb.
- Symbol switches and probes (A7). Each source file has counts of `applicable`
  calls, `@eval` calls, and comparisons of a `kind` field against symbols, negated
  comparisons included. These counts can decrease and never increase.

`Commons` holds only what is defined once and used by several owners. An owner is
the root module, a top-level submodule or a package extension. Six guards keep
helpers from accumulating elsewhere.

- Commons admission (C1). Each function, type and constant defined in `Commons` is
  exported or declared `public`, has a docstring and is named in a test file under
  `test/unit/commons/`. A hook is a `Commons` function with a method defined outside
  `Commons`. The user API is the `Commons` names that the root module exports or
  declares `public`. Unless a name is user API, methods of at least two owners
  outside `Commons` use it. Unless it is a hook or user API, it has an entry in the
  reserved vocabulary of C2. A use is a reference in lowered code, a method
  signature, a method extension, a supertype or a field type. Uses by a public
  `Commons` definition count for the definitions it uses. A definition with one
  remaining owner moves to that owner. Algorithm-local helpers are nested functions
  inside the public definition.
- Reserved vocabulary (C2). `VOCABULARY` maps each `Commons` public name to the
  definition names it reserves. Outside `src/commons/`, no function, constant or
  assigned local takes a reserved name. A local assigned from a call to the
  reserving definition is exempt, as in `μ0 = vacuum_permeability(T)`.
- Literal fingerprints (C3). Outside `src/commons/consts.jl`, no numeric literal
  contains the digits `8854187` or `299792458`, and no statement that names `π`
  contains `1e-7` or a power of a base containing 10 with exponent `-7`.
- Clones (C4). Function bodies are tokenized with each identifier replaced by its
  role (call, field or name), and compared through runs of 8 tokens. Statements of
  the form `x || throw(...)` or `x && throw(...)`, with `throw`, `rethrow` or
  `error`, are left out. No function body under `src/` or `ext/` contains 75 % of the
  runs of a `Commons` public function body of at least 30 tokens. Hooks are not
  compared.
- Tiny helpers (C5). A tiny helper is a module-level function that is neither
  exported nor public, has one method, has at most three body statements, is
  referenced by exactly one method and is named in no test file. Files under
  `test/quality/` and `test/tools/` are not test files here. An exact forwarder,
  whose body is one call that passes its own arguments unchanged and in order,
  keywords included, counts whatever its references. Inline a new tiny helper,
  import the owner's definition, or test it directly.
- Root freeze (C6). The number of functions, types and constants defined by the root
  module does not grow. New shared definitions go to `Commons` under C1.

Each guard is a function of its inputs. The negative controls (A8 and C7) apply each
guard to a probe package with planted violations and to a clean probe package. A
guard reports exactly the planted violations and nothing in the clean probe.

`test/quality/architecture_baseline.toml` records the violations present when the
guards were introduced, with one table per guard. A table that no guard owns fails. Keys contain no line numbers. An
unlisted violation or a count above its entry fails. A listed entry without a live
violation, or a count below it, also fails. A change that removes a violation deletes
or lowers its entry. `test/tools/architecture_inventory.jl` prints the live inventory
in the baseline format. Fix a new violation in the source. Never add it to the
baseline.

Entries can be deleted or lowered, and none can be added or raised. Before the tests,
the quality CI job runs `test/tools/baseline_ratchet.jl` against the pull request
base or the previous push. An added entry or a raised count fails the job. The
ratchet first applies to the earlier keys the file renames that git detects, and the
module renames of renamed entry files `<Module>.jl` that declare their module. A key
renamed without a matching git rename counts as added. A table absent from the
earlier baseline belongs to a guard introduced since, and the ratchet lists it
without comparing its keys. Run
`julia test/tools/baseline_ratchet.jl HEAD` to compare local changes with the last
commit.

## Developer paths

- [Extension API](extensions.md) lists the equation methods and definition types omitted
  from the user API reference.
- [Computational engine](engine.md) covers formulations, options, supplemental
  calculation output, and external implementations.
- [Makie plotting](plotting.md) covers the small high-level API and native ownership.
- [Conventions](conventions.md) defines placement, dispatch, naming, and
  docstring rules.

## Testing requirements

### Release status and regressions

The package has no published stable API. Development changes update the
implementation, callers, and relevant tests together.

In this repository, a regression test protects against a bug detected in a
published stable release. It identifies the reported tracker issue and affected release, together with the protected behavior. Tests of current development behavior belong
to the behavior, mathematical implementation, integration, or architecture
suites, according to what they exercise.

Tests follow the current API and its intended behavior. Preserve independent
scientific expectations when updating callers. Architectural tests verify each module through its consumers and native interfaces.

### Verification coverage

The harness checks that the current implementation works and follows the
codebase standards. It exercises actual public workflows and their owned
kernels: calculations, dispatch, units, ordering, data transport, errors,
resource handling and side effects. Tests use current input builders and small,
distinguishable examples with an expected result or rejection.

Mathematical implementation tests remain appropriate. Checking an implemented
formula against a directly calculable result, a matrix reduction against a
constrained solve, or a derivative of a simple function verifies code correctness.
It does not certify the physical applicability of the model. Imports, shapes,
finiteness and round trips provide useful limited checks. They do not replace
value assertions where the implemented behavior has a decidable expectation.

Float32 support means that valid ordinary inputs containing Float32 values work
without type-induced crashes. It includes no additional numerical accuracy
promise. Do not impose a high-precision reference target on Float32 and then
change owned numerical code to meet it. No Float32-specific widening, compensated
arithmetic or precision infrastructure is justified by such a test. Preserve
existing dispatch, hooks, uncertainty and supported types. Keep API conversions
when an actual API or dependency interface requires them. A requested tolerance
does not create an accuracy guarantee.

Scientific validity, comparison of physical approximations, broad accuracy
claims, convergence research and scientific acceptance are outside the test
harness. Scientific acceptance belongs to the researcher's interpretation, not
the test runner. Such research is not a required CI job, and unfinished scientific
evidence is not a failing software check. Apply this distinction to passing and failing
experiments alike.

A fabricated failure is as unacceptable as a fabricated success. Verify an
assertion's premise before treating its outcome as a product defect. Correct or
remove a defective criterion with its reason. Do not preserve it merely because
it was fixed before execution. Do not change inputs, outputs or tolerances just
to obtain a pass. Preserve observed differences and distinguish assertion
failures from execution errors or unavailable verification.

### Architecture and coverage

Architectural tests protect current responsibilities through real behavior and
native method and interface checks: owner-local dispatch, fixed report stages,
observation and table owners, validated inputs, optional integrations and
caller-owned state. A conforming new leaf must work through the actual composed consumer.

The existing source-amended production line-coverage gate remains at 95%, with
its current `src/` and `ext/` inventory. Measure executed code, then cover actual
missing behavior in its existing test owner. Assertion counts, obsolete guards,
denominator changes and fabricated expectations cannot satisfy this objective.
Keep incomplete executions and real bugs visible. Coverage does not turn them
into successes.

## External interface methods

`test/quality/explicit_imports.jl` loads the numerical, XLSX and Cairo adapters
explicitly. All mechanical ownership and import checks remain active. Each
exception identifies an exact consumer, owner and name for a documented
upstream interface that lacks a Julia `public` annotation. Access to
package-private names remains prohibited.

The following external interfaces are used by the package:

| Access | Supported use |
| --- | --- |
| `CairoMakie.activate!` | Documented [backend activation](https://docs.makie.org/stable/explanations/backends/cairomakie.html). Cairo adapter only. |
| `Base.IOError` | Native I/O exception, including [filesystem errors](https://docs.julialang.org/en/v1/base/file/). Renderer export errors and FEM recovery from file and process errors. |
| `Base.unalias` | Documented native preventative-copy operation in Julia 1.12's `base/abstractarray.jl`. Used only by Engine's allocating Kron entry point to preserve source and destination aliasing. The workspace path uses separate preallocated buffers. |
| `Base.get_extension` | Identifies the loaded Cairo extension. Its public `CairoMakie` binding supplies the native save backend. |
| `Makie.automatic` | Documented [native attribute default](https://docs.makie.org/stable/api). Renderer only. |
| `Makie.current_backend` | Documented [backend-dependent API default](https://docs.makie.org/stable/api). Renderer only. |
| `Makie.get_ticks`, `Makie.get_tickvalues` | Documented [axis extension hooks](https://docs.makie.org/stable/reference/blocks/axis.html). Renderer only. |
| `Makie.pseudolog10` | Documented [axis scale](https://docs.makie.org/stable/reference/blocks/axis.html). Renderer only. |
| `Makie.inverse_transform` | Documented [custom axis scale interface](https://docs.makie.org/stable/reference/blocks/axis.html#xscale). The renderer applies axis margins to full uncertainty bounds in the selected scale, then maps them back without a second inverse-scale mapping. |
| `Makie.CategoricalConversion` | Documented [categorical axis conversion](https://docs.makie.org/stable/reference/generic/dimensional/). Renderer assembly axes only. |
| `Makie.defaultlimits` | The documented native scale-default hook named by Makie's Axis attributes. The renderer queries it only when an empty axis changes scale. The native scale then supplies its valid interval. |
| `propertynames`, `Makie.default_theme` | Axis attributes use `propertynames` and the native `palette` keyword. Scatter attributes use the exported [`default_theme`](https://docs.makie.org/v0.24/explanations/recipes) method. |
| `LegendElement.plots` | Associates legend glyphs with source plots through the documented [LegendElement extension interface](https://github.com/MakieOrg/Makie.jl/blob/v0.24.13/Makie/src/makielayout/types.jl). |
| `on`, `off`, `ObserverFunction.observable` | Manage visibility subscriptions using the documented [`ObserverFunction.observable`](https://juliagizmos.github.io/Observables.jl/stable/#Observables.ObserverFunction) notification target. |
| `Makie.fast_string_boundingboxes(Text)` | Public attribute subscriptions query the documented [text bounds](https://github.com/MakieOrg/Makie.jl/blob/v0.24.13/Makie/src/basic_recipes/text.jl), retaining marker-space extents. |
| `GridLayoutBase.remove_from_gridlayout!` | Retained as the [maintainer-prescribed nested-layout removal](https://discourse.julialang.org/t/makie-removing-gridlayouts/103935) workaround. It is not an exported stable API. Only this renderer call is admitted. Existing legend recreation and layout tests protect it. |

These integration methods are scoped to the supported Makie 0.24 family.
The layout workaround needs review when that compatibility range changes.
The owned legend, uncertainty visibility, series style and colorbar endpoint tests
exercise the actual affected paths. Behavioral compatibility requires execution of the affected paths. Negative gate controls reject unlisted renderer
internals, different consumers and package-owned private calls.
