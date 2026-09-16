# Grammar invariants

LineCableModels has three global calculation roles: a complete problem, a
formulation that selects how to calculate it, and a completed result. Package
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
tree. UQ is therefore not a separate calculation supertype. Its completed
collections sit below `AbstractUncertaintyResult`, while deterministic finite
collections sit below `AbstractParametricResult`.

`LineParamsDomain` is independent of this calculation grammar. It tags the
physical coordinate system of a completed line-parameter matrix:

```@example grammar_type_trees
print(Main.DocumentationTrees.type_tree(LineCableModels.LineParamsDomain))
```

## Fixed actions

One declarative action owns a fixed sequence selected by an abstract definition
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
Grammar owns one `observables(source, requests::Tuple; ...)` publication
method. Standalone arrays and external result types do not need a shared
abstract observation type.

## Architecture checks

The architecture that the quality, core and integration tests protect requires:

- the action and its abstract root to belong to the same module;
- one public action method, with no more-specific definition methods;
- every fixed stage to remain visible in the declared action;
- every package-owned concrete definition to implement its required stages;
- every package-owned definition to remain below its declared root;
- one Grammar-owned observation publication method with positional requests;
- wide scientific tables whose quantity and unit metadata remain attached to
  their columns;
- no calls to another package module's private functions or types, including
  through aliases, and no private wrappers whose methods only forward unchanged
  arguments to the same package-owned function. Owner-local numerical kernels
  and type-dispatch branches remain valid.

A new definition must satisfy these standards. Native interface/import checks
and tests through actual composed consumers provide the protection; the suite
does not claim a universal source analysis of every method or forwarding wrapper.

Integration tests also count actual blueprint lowering calls: one per selected
design point, shared across its formulation alternatives and frequency sweep.
They check local-before-earth assembly order at each frequency and explicit
phase-to-modal result transport without another geometry lowering. The visual
suite applies the ownership checks to the loaded Makie extensions and verifies
that material colors consume `Material` or `EarthLayer` objects directly.

InputValidation tests check owner dispatch, required interfaces, unchanged valid
inputs and rejection of damaged inputs through validation and computation. The
standards require checks to remain with the owning validator; representative
behavioral tests do not establish that every method body has been inspected.

The native FEM environment uses pinned GetDP 3.5.0. Keep tests of its actual
execution, extraction, terminal identity, material/option transport and resume
behavior. Physical cross-backend accuracy and domain-convergence acceptance
belong to explicit Gauntlet/research work under the testing policy below. See the
[test commands](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/README.md) for the implementation checks and their
execution environments.

## Developer paths

- [Extension API](extensions.md) lists the hooks and definition types omitted
  from the user API reference.
- [Computational engine](engine.md) covers formulations, options, supplemental
  calculation output, and external implementations.
- [Makie plotting](plotting.md) covers the small high-level API and native ownership.
- [Conventions](conventions.md) defines placement, dispatch, naming, and
  docstring rules.

## Testing policy

### Release status and regressions

The codebase is moving toward its first stable release. Everything currently
pushed to `main`, including the `0.1.0` tag, is unreleased. The intended `0.2.0`
candidate does not create a previous stable release or a compatibility promise.

In this repository, a regression test protects against a real bug detected in
released stable code. It identifies the reported tracker issue, affected stable
release and protected behavior. Current development defects are bugs to fix;
they are not regressions against an unreleased prototype. Ordinary behavior,
mathematical implementation, integration and architecture tests need no invented
issue. Existing useful tests remain under those purposes. Stable publication
does not retroactively turn WIP tests into bug-regression tests.

Deliberate API and architectural changes update the implementation, callers and
relevant tests together. Do not restore obsolete APIs or preserve prototype
outputs to satisfy tests. A test whose only purpose is to prove an old name
vanished or a new spelling appeared is not an architectural safeguard.

### What the harness checks

The harness checks that the current implementation works and follows the
codebase standards. It exercises actual public workflows and their owned
kernels: calculations, dispatch, units, ordering, data transport, errors,
resource handling and side effects. Tests use current input builders and small,
distinguishable examples with an expected result or rejection.

Mathematical implementation tests remain appropriate. Checking an implemented
formula against a directly calculable result, a matrix reduction against a
constrained solve, or a derivative of a simple function checks code correctness.
It does not certify the physical applicability of the model. Imports, shapes,
finiteness and round trips provide useful limited checks; they do not replace
value assertions where the implemented behavior has a decidable expectation.

Float32 support means that valid ordinary inputs containing Float32 values work
without type-induced crashes. It carries no additional numerical accuracy
promise. Do not impose a high-precision reference target on Float32 and then
change owned numerical code to meet it. No Float32-specific widening, compensated
arithmetic or precision infrastructure is justified by such a test. Preserve
existing dispatch, hook, uncertainty and type contracts; keep boundary conversions
when an actual API or dependency interface requires them. A requested tolerance
does not create an accuracy guarantee.

Scientific validity, comparison of physical approximations, broad accuracy
claims, convergence research and scientific acceptance are outside the test
harness. Gauntlet executes calculations and reports comparisons and timings;
scientific acceptance belongs to the researcher's interpretation, not its runner.
Such research is not a required CI job, and unfinished scientific evidence is not a failing code
contract. Apply this boundary to passing and failing experiments alike.

A fabricated failure is as unacceptable as a fabricated success. Verify an
assertion's premise before treating its outcome as a product defect. Correct or
remove a defective criterion with its reason; do not preserve it merely because
it was fixed before execution. Do not change inputs, outputs or tolerances just
to obtain a pass. Preserve observed differences and distinguish assertion
failures from execution errors or unavailable verification.

### Architecture and coverage

Architectural tests protect current responsibilities through real behavior and
native method/interface checks: owner-local dispatch, fixed report stages,
observation and table boundaries, validated inputs, optional integrations and
caller-owned state. A conforming new leaf must work through the actual composed
consumer. Unrelated old helper names, private storage layouts and incidental
source expressions are not substitutes for these checks. Keep the standards
in this guide and the conventions; do not invent a second architecture framework.

The existing source-amended production line-coverage gate remains at 95%, with
its current `src/` and `ext/` inventory. Measure executed code, then cover actual
missing behavior in its existing test owner. Assertion counts, obsolete guards,
denominator changes and fabricated expectations cannot satisfy this objective.
Keep incomplete executions and real bugs visible; coverage does not turn them
into successes.

The completed fixture reset and useful rendering, reporting and harness repairs
stand. Do not repeat the reset, restore legacy expected output, or turn completed
calibration into a recurring obligation. Prior completion does not justify a
test or numerical change whose requirement was unsupported, including the
accuracy-driven Float32 surface-evaluation change.

### Numerical snapshots after stable publication

Numerical snapshot testing is deferred until after the first stable publication.
The user will choose a small number of Gauntlet artifacts. Do not select, create,
refresh or approve numerical baselines, add a snapshot dependency, or activate a
snapshot CI job before then. The existing inactive provision is sufficient now;
its empty reference list is not a prerelease failure or missing approval task.

The future check compares each selected backend with its own retained output
across revisions. It records case/settings identity, frequencies, terminal order,
units/basis, relevant returned quantities and execution provenance. Compare
individual meaningful components with explicit continuity tolerances; an RMS
summary alone must not hide a local change. Do not snapshot private workspaces,
internal layouts or incidental paths. Keep a designated reference fixed until
an explicit user-approved update; no automatic refresh after a passing push.

These snapshots detect behavioral change, not scientific validity. They do not
replace architectural tests or create a new scientific approval process. Reuse
existing storage and test owners when the user activates this later work.

## External interface contracts

`test/quality/explicit_imports.jl` loads the numerical, XLSX and Cairo adapters
explicitly. All mechanical ownership/import checks remain active. Its exact
consumer/owner/name exceptions recognize documented upstream interfaces that
lack Julia `public` annotations; they grant no package-private access.

The inspected external accesses have these dispositions:

| Access | Contract and disposition |
| --- | --- |
| `CairoMakie.activate!` | Documented [backend activation](https://docs.makie.org/stable/explanations/backends/cairomakie.html); Cairo adapter only. |
| `Base.IOError` | Native I/O exception, including [filesystem errors](https://docs.julialang.org/en/v1/base/file/); renderer export error handling only. |
| `Base.require` | Removed from the renderer. `Base.get_extension` identifies the loaded Cairo extension, whose public `CairoMakie` binding supplies the native save backend. Export never loads packages. |
| `Makie.automatic` | Documented [native attribute default](https://docs.makie.org/stable/api); renderer only. |
| `Makie.current_backend` | Documented [backend-dependent API default](https://docs.makie.org/stable/api); renderer only. |
| `Makie.get_ticks`, `Makie.get_tickvalues` | Documented [axis extension hooks](https://docs.makie.org/stable/reference/blocks/axis.html); renderer only. |
| `Makie.pseudolog10` | Documented [axis scale](https://docs.makie.org/stable/reference/blocks/axis.html); renderer only. |
| `Makie.inverse_transform` | Documented [custom axis scale contract](https://docs.makie.org/stable/reference/blocks/axis.html#xscale); the renderer applies axis margins to full uncertainty bounds in the selected scale, then maps them back without a second inverse-scale registry. |
| `Makie.attribute_names` | Removed. Axis uses its `propertynames` interface plus the native `palette` keyword; Scatter uses the exported [`default_theme`](https://docs.makie.org/v0.24/explanations/recipes) contract. |
| `Makie.get_plots` | Removed. Legend glyphs use the `plots` vector required by the [LegendElement extension contract](https://github.com/MakieOrg/Makie.jl/blob/v0.24.13/Makie/src/makielayout/types.jl). This is the documented source association, not arbitrary private-field inspection. |
| `Makie.get_plot_visibilities` | Removed. Native `on`/`off` and the documented [`ObserverFunction.observable`](https://juliagizmos.github.io/Observables.jl/stable/#Observables.ObserverFunction) supply the notification target without changing visibility. |
| `Makie.fast_string_boundingboxes_obs` | Removed. Public attribute subscriptions query the documented [`fast_string_boundingboxes(Text)`](https://github.com/MakieOrg/Makie.jl/blob/v0.24.13/Makie/src/basic_recipes/text.jl) result, preserving marker-space extents without the internal observable helper. The exact documented query is allowed for the renderer. |
| `GridLayoutBase.remove_from_gridlayout!` | Retained as the [maintainer-prescribed nested-layout removal](https://discourse.julialang.org/t/makie-removing-gridlayouts/103935) workaround. It is not an exported stable API. Only this renderer call is admitted; existing legend recreation/layout tests protect it. |

These integration contracts are scoped to the supported Makie 0.24 family.
The layout workaround needs review when that compatibility range changes.
The owned legend, uncertainty visibility, series style and colorbar endpoint tests
exercise the actual affected paths. Source declarations alone do not establish
behavioral compatibility. Negative gate controls reject unlisted renderer
internals, different consumers and package-owned private calls.
