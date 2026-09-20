# Conventions

## Documentation language

Use US English in docstrings and project documentation. Preserve API identifiers,
file paths, quoted titles, and proper names. State what a calculation does,
define its physical inputs and units, and distinguish implemented behavior from
approximations and literature results. Omit promotional claims and process
terminology when a concrete description suffices.

LineCableModels follows the SciML formatter style. Run:

```julia
using JuliaFormatter
format(".")
```

before committing Julia source.

## Native-first semantic economy

Use Julia's existing meanings before adding package vocabulary. Check, in
order, whether the operation is already expressed by Core, Base, a standard
library, a direct dependency, or an existing package generic.

A new name must identify a domain action, an invariant, a numerical method, an
external format, or real state. A name that only forwards arguments, reads one
field, repacks a tuple, or merges defaults does not add meaning.

| Prefer | Avoid |
|:--|:--|
| `Base.length`, iteration, indexing, `show`, and `showerror` | package copies of collection and display operations |
| constructors, `promote`, and `convert` | a second conversion vocabulary |
| another method on an owned generic | a synonym that forwards to that generic |
| dispatch on a scientific type | a symbol switch or dictionary that repeats dispatch |
| a local method beside its owner | a global `helpers` or `utils` bucket |

For example, a result that is a finite collection implements Julia's collection
methods:

```julia
Base.length(result::MyResult) = length(result.values)
Base.getindex(result::MyResult, index::Integer) = result.values[index]
Base.iterate(result::MyResult, state...) = iterate(result.values, state...)
```

It does not add `get_results`, `result_count`, and `iterate_results` as parallel
names.

Small methods are appropriate when each method owns a dispatch choice or an
invariant. Tiny forwarding helpers that only rename another operation are not.
Do not add speculative compatibility shims, runtime `eval`, exception-driven
feature tests, or lookup tables that duplicate Julia methods.

Mutating functions use `!` only when the operation can mutate an argument or
externally visible state.

## Dispatch-driven fixed actions

Use one public action when an operation has one required stage order. The
action method shows the complete sequence, and concrete definition types add
methods for the stages they own.

```julia
function process(definition::AbstractDefinition, source)
    selected = select(definition, source)
    product = build(definition, selected)
    return Product(product)
end
```

Good definition types are concrete and passive: their fields state scientific
choices or completed configuration. Runtime work belongs to stage methods.

```julia
struct FrequencyReport <: AbstractReportDefinition
    requests::Tuple
end

select(definition::FrequencyReport, source) =
    observables(source, definition.requests)
```

Do not replace type-directed stages with a generic `mode`, `style`, or options
dictionary interpreted by one large switch. Normalize public symbols once at
the owner's entry point when symbols are part of the public syntax, then use
the existing `Val` method family:

```julia
owned_action(selector::Symbol, args...; kwargs...) =
    owned_action(Val(selector), args...; kwargs...)

owned_action(::Val{:example}, args...; kwargs...) = ...
```

Use an explicit no-op method only when doing nothing is a valid stage result.
Reject unsupported definition/source pairs through required stage dispatch
before partial work. Introduce a mutable context only when several stages
genuinely share buffers, resources, or evolving state. CI checks the fixed
actions listed in [Grammar invariants](developers.md) directly; runtime
metadata that merely repeats their method definitions is not part of the
grammar.

## Ownership-centered recursive module layout

Place code first by the owner that defines when it changes, then by its precise
responsibility. Files, directories, and Julia modules solve different
problems:

- a file separates one responsibility within an owner;
- a directory groups several responsibilities that still belong to one owner;
- a submodule supplies a separate namespace, dependency set, or stated
  interface;
- a package is warranted only when the code has independent users and releases.

Grow code recursively:

```text
single responsibility in owner file
└─ several responsibilities in owner directory, same module
   └─ separate namespace or dependency set in child module
      └─ independent package only when independently consumed
```

A module entry file is an index. It contains the module description, explicit
imports, public names, includes in dependency order, and deliberate child
reexports. Constructors, algorithms, validation, plotting descriptions, and
format translations belong in focused files selected by the owner.

Place a method according to the reason it changes. A method that exposes an
Engine result through `observe` belongs with that result. A method that draws a
native figure with Makie belongs in the Makie extension. Scientific
observations and physical preview geometry remain with their owner; plot
request normalization, presentation groups, palettes, Makie blocks, layouts,
widgets, callbacks, and backend activation do not. A method that parses an
external file belongs with the format owner.

Optional dependencies remain in package extensions. Core source may define
package-neutral requests and completed values, but it does not import Makie,
XLSX, Measurements, or Distributions.

Prefer conceptual groupings such as:

```text
Owner
├─ types
├─ interfaces
├─ constructors
├─ action
├─ Base protocols
└─ owned optional translations
```

Avoid global mechanism-first trees such as `types/services/managers/handlers`,
one module per file or per type, empty directory scaffolds, and a common base
file that accumulates unrelated methods. Do not split a fixed call sequence
across files merely to make each stage visually separate.

## Scientific reads, tables, and plots

`observe` reads native scientific values. `ObservedResult` captures detached
products for one completed gridpoint. Ordinary collections lift that constructor.
Reports, tables, plots, exports, and saved report inspection consume observations.

```text
completed numerical result + completed comparisons + recorded timings
→ ObservedResult(gridpoint, quantities, errors, timings)
→ retained selection → tables / existing Makie renderers / persistence
```

The result owner implements `observation_quantity` and declares its requests.
`Grammar.observation_gridpoint` reads the description captured in completed result
storage. Completion captures actual physical inputs, formulation selections and
controls, coordinates, and original point identity. Basic descriptions are always
retained. Observation and presentation never reconstruct a problem from lazy axes.

The complete primary pairs are R/X, magnitude/angle of Z, or R/L for Z, and G/B,
magnitude/angle of Y, or G/C for Y. Defaults are R/X and G/B. The strict constructor
rejects an incomplete pair. Raw plotting conveniences use the same request
normalizer with `complete_pairs=true`; `plot(line; ydata=(R,))` retains R/X and G/B
and displays R. `plot(observed; ydata=(R,))` selects retained R only. Different
representations require a new explicit construction. No consumer derives an
absent quantity, repeats clipping, or acquires a raw source.

UQ acquisition belongs to the UQ owner: `ObservedResult(uq, point, requests)` joins
that point's primary values and requested statistics, samples, or histograms.
Product requests omit the point index. Mean/std estimates and precomputed
histogram/CDF/Q–Q coordinates remain ordinary records in `quantities`; sampling
information belongs to `gridpoint.sampling`. Uncertain primary values in an ordinary
collection use the primary owner and preserve their dependencies.

```julia
observed = ObservedResult(parameters, (R, L, G, C))
tables = ReportBuilder.tabulate(observed)  # tables.Z.R, tables.Z.L, tables.Y.G, tables.Y.C
plot(observed; ydata=(R,))
```

Every quantity has its own table. A full n×n matrix on m frequencies has m rows and
1+n² columns, in row-major coefficient order; both off-diagonals remain present.
Sparse and diagonal requests preserve original indices and matrix extent.
`DataFrame(observed)` is a diagnostic long view. `ObservedResult` is not a
Tables.jl table. XLSX writes one numeric values/std workbook per point and quantity;
native observation persistence preserves precision and uncertainty dependencies
across the complete candidate/reference archive.

Explicit benchmark comparison precedes construction:

```julia
completed = compare(reference, candidates, [R, L]; bands=(:all,))
points = observables(candidates; comparisons=completed, timings=recorded_timings)
reference_point = ObservedResult(reference)
artifact = report(BenchmarkTableDefinition(), points; reference=reference_point)
plot(artifact; ydata=(R,))
```

Comparison records join by original candidate identities. A reference remains
outside the candidate collection. External data needs explicit retained identity
for a benchmark; absent physical descriptions remain explicitly absent. Arithmetic
checks cover coordinates, dimensions, units, basis, and frequency agreement.
Scientific comparability is the caller's responsibility; no interpolation occurs.
`ReportArtifact.observed` and `.reference` contain only observations; `.tables`
contains the organized tables. Raw conveniences construct and delegate once.

`observation_groups` is the shared grouping owner. Eligibility requires the same
physical point, relevant selected formulas and controls, quantity/statistical
meaning, units, uncertainty interpretation, and coordinates. Exact numerical and
dependency agreement verifies that eligibility. It does not discover equivalence.
All group identities, original observations, tables, and files remain available.

## Numerical reporting and current behavior

An available scalar is engineering zero precisely when
`abs(nominal(value)) <= cutoff`. Nonfinite and unavailable values are separate.
Defaults per meter are R=1e-10 Ω/m, L=1e-15 H/m, G=1e-12 S/m, C=1e-16 F/m;
X and B use 2πf times their L and C cutoffs. Total quantities scale by a retained
physical length or require explicit cutoffs. L/C are unavailable at DC. No inferred
`eps`, matrix-norm, largest-coefficient, or uncertainty contribution is added.
Complex zero requires both Cartesian components to be zero. Polar products come
from the original complex values. Recentring preserves uncertainty dependencies
and every spread. Undefined first-order magnitude retains its components and zero
nominal magnitude with an explicit reason; its value and phase remain missing.

Comparison classifies original operands first. Any ineligible sample makes both
RMS metrics missing for that coefficient/band. Otherwise existing RMS mathematics
applies to original values. Small and zero errors remain valid; operand cutoffs
are never applied to errors. Actual cutoffs, units, selections, settings, and
missing reasons are retained.

Recorded timings remain associated with the original candidate and separate
reference in report tables. Equal elapsed times do not identify a shared event.
A measurement for a whole calculation retains that scope when several observed
points carry it; reporting does not invent per-point timings.

Intentional changes replace current internal behavior. Persistence and consumer
dependencies do not require generations of that behavior. Do not introduce a
renamed policy counter, code fingerprint, compatibility branch, or maintenance
instruction to recreate that obligation. File checksums and source revisions
retain their actual integrity and source-identification meanings.

| Removed owned representation | Scientific meaning and current owner |
| --- | --- |
| `ObservationPublication`, `publication_table`, parallel flattened columns | One `ObservedResult` per point; `ReportBuilder.tabulate` creates quantity tables on demand |
| Publication `contract` and column contracts | Explicit quantity, basis, units, coordinates, cutoffs, availability, and reasons in each quantity record |
| Publication `provenance` and raw source references | Captured inputs, formulations, original identities, sampling information, and completed measurements in the four observation sections |
| `getdp_provenance` | `getdp_selection` in FEM completion, retained details, recovery, and their tests |
| Resolution revision and policy-generation checks | Deleted; actual applied numerical settings are retained |
| `ReportArtifact.published` / `.table` | `.observed`, separate `.reference`, and `.tables`; no compatibility getters |
| Raw result-specific report/renderer preparation | Constructor conveniences delegate to the common observed workflow |

## Text display and table boundaries

Human inspection (`show`), scientific extraction (`observe`), detached acquisition
(`ObservedResult`/`observables`), and tabulation (`tabulate`) are separate actions.
Each public owned type defines `summary`, two-argument `show`, and text/plain
`show`. Display reads stored state only; it must not run builders, comparisons,
solvers, or lazy-grid materialization. Ordinary report display renders retained
tables and never constructs a figure implicitly. Existing Makie axes and the
plot window own drawing, controls, layout, and backend behavior.

## Docstrings

Docstrings use DocStringExtensions abbreviations so declarations remain aligned
with the implementation.

### Placement and content

- Place a docstring immediately before the documented module, type,
  constructor, function, or constant.
- Use triple double quotes, except for concise field and constant docstrings.
- Describe implemented behavior. Do not infer equations, units, or defaults
  from a name.
- State each fact once.
- Link related local bindings inline when the relationship helps the reader.

Use `@doc` for an inner constructor written inside a `struct`. An outer
constructor at module scope uses an ordinary preceding docstring.

### Physical quantities and equations

State the SI unit for every physical argument, return value, field, and
constant. In Julia docstring source, escape square brackets and LaTeX commands:

```julia
"Series resistance `\\[Ω/m\\]`."
"Relative permeability `\\[dimensionless\\]`."
```

Comments inside Julia examples use ordinary brackets, such as `# [m]`.

When code directly evaluates a physical law, approximation, or reduction that
matters to the method's meaning, include the equation and define its symbols:

````julia
"""
$(TYPEDSIGNATURES)

Return the series impedance:

```math
Z(f) = R + \\mathrm{j} 2 \\pi f L,
```

where ``R`` is resistance, ``L`` is inductance, and ``f`` is frequency.
"""
````

Accessors, forwarding methods, and bookkeeping functions do not need a
mathematical section unless they evaluate the documented expression.

### DocStringExtensions abbreviations

- `$(TYPEDSIGNATURES)` is the default opening for functions and constructors.
- `$(SIGNATURES)` is suitable when typed signatures obscure the public call.
- `$(FUNCTIONNAME)` keeps executable examples aligned with renames.
- `$(TYPEDEF)` inserts a type declaration.
- `$(TYPEDFIELDS)` inserts fields, declared types, and field docstrings.
- `$(FIELDS)` omits declared field types when they would distract from the
  public meaning.
- `$(METHODLIST)` is reserved for a multi-method interface whose purpose is to
  list implementations.
- `$(IMPORTS)` and `$(EXPORTS)` maintain a module inventory.

Do not repeat generated text by hand.

### Function structure

Use this section order, omitting sections that add no information:

1. description and any implemented equation;
2. `# Arguments`;
3. `# Keywords`;
4. `# Returns`;
5. `# Notes` for assumptions or limitations;
6. `# Errors` for deliberate exceptions;
7. `# Examples`.

````julia
"""
$(TYPEDSIGNATURES)

Describe the implemented operation.

# Arguments

- `value`: Physical input `\\[unit\\]`.

# Keywords

- `basis`: `:pul` or `:total`. Default: `:pul`.

# Returns

- Completed value in `\\[unit\\]`.

# Errors

- Throws `ArgumentError` when `basis` is unsupported.

# Examples

```jldoctest
result = $(FUNCTIONNAME)(1.0; basis=:pul) # [unit]
@assert isfinite(result)
# output
```
"""
````

List arguments in declaration order and document each returned tuple member.
Prefer `jldoctest` for a self-contained public example. Examples requiring an
external executable, graphical interaction, network access, or repository
fixtures belong in the developer guides.

### Type and module structure

Use `$(TYPEDEF)` and `$(TYPEDFIELDS)` for a type. Put a concise docstring above
each field:

````julia
"""
$(TYPEDEF)

Represent a cable section.

$(TYPEDFIELDS)
"""
struct CableSection{T <: Real}
    "Section thickness `\\[m\\]`."
    thickness::T
end
````

A module docstring begins with the indented module name, states its purpose,
then uses `$(IMPORTS)` and `$(EXPORTS)` when those lists aid the reader. A
physical constant uses a concise single-line docstring with its symbol and SI
unit.

## Repository practice

After the first stable publication, versions follow [Semantic Versioning](https://semver.org/).
The current 0.2.0 candidate establishes the initial intended contract; development API
renames are not regressions merely because an earlier spelling existed.

Commit subjects use scoped Conventional Commits, begin with a lowercase
description, and stay within 72 characters:

```text
fix(engine): reject unsupported formulation options
```

Every change includes tests at the closest relevant scope. Core tests do not
load optional packages. Rendering activation and dependency installation are distinct. CairoMakie is a current
package dependency; its rendering extensions activate when loaded. Rendering and other
extension paths also run in their dedicated test environments. Public examples should be executable and self-contained.

## Testing policy

The [developer testing policy](developers.md#testing-policy) is authoritative for
release status, regression terminology, test scope, architecture and the unchanged
95% coverage gate. Everything currently on `main`, including `0.1.0`, is unreleased.
The harness checks implementation correctness and architectural conformance;
scientific acceptance is outside it. User-selected numerical snapshots are
deferred until after the first stable publication. Do not duplicate that policy
as a separate validation or baseline-approval scheme.
