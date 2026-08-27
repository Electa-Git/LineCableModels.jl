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

## Declarative actions

Three actions own a fixed sequence selected by an abstract definition type:

| Action | Definition root | Fixed sequence |
|:--|:--|:--|
| `PlotBuilder.make_render` | `AbstractPlotDefinition` | `entitle → parse → resolve → fetch → finish` |
| `ReportBuilder.report` | `AbstractReportDefinition` | `entitle → select → tabulate → illustrate → encode → write → finish` |
| `ParametricBuilder.project` | `AbstractProjectionDefinition` | `entitle → select → derive → materialize → finish` |

```@example grammar_type_trees
using LineCableModels.PlotBuilder: AbstractPlotDefinition
using LineCableModels.ReportBuilder: AbstractReportDefinition
using LineCableModels.ParametricBuilder: AbstractProjectionDefinition

println("Plot definitions")
print(Main.DocumentationTrees.type_tree(AbstractPlotDefinition))
println("\nReport definitions")
print(Main.DocumentationTrees.type_tree(AbstractReportDefinition))
println("\nProjection definitions")
print(Main.DocumentationTrees.type_tree(AbstractProjectionDefinition))
```

Each action has one method at its declared abstract root. Concrete definitions
implement stage methods. They do not specialise the public action itself.
Required stages are declared with RequiredInterfaces where the type family
admits that form; PlotBuilder checks its definition-owned stages directly.

Observation uses a different structure. Result owners add `observe` methods.
Grammar owns one `observables(source, requests::Tuple; ...)` publication
method. Standalone arrays and external result types do not need a shared
abstract observation type.

## CI hard gates

The `Quality contracts` CI job rejects changes that violate these rules. Its
checks require:

- declared action metadata and a root present in the action signature;
- the action and its abstract root to belong to the same module;
- one public action method, with no more-specific definition methods;
- every fixed stage to remain visible in the declared action;
- every package-owned concrete definition to implement its required stages;
- every package-owned definition to remain below its declared root;
- one Grammar-owned observation publication method with positional requests;
- wide scientific tables whose quantity and unit metadata remain attached to
  their columns.

These checks run as test failures. A new definition must satisfy them before it
can enter the maintained type family.

## Developer paths

- [Extension API](extensions.md) lists the hooks and definition types omitted
  from the user API reference.
- [Computational engine](engine.md) covers formulations, options, supplemental
  calculation output, and external implementations.
- [PlotBuilder guide](plotbuilder.md) covers detached recipes and Makie drawing.
- [Conventions](conventions.md) defines placement, dispatch, naming, and
  docstring rules.
