# Conventions

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
externally visible state. The suffix states behaviour; it is not decoration.

## Dispatch-driven fixed actions

Use one public action when an operation has one required stage order. The
action method shows the complete sequence, and concrete definition types add
methods for the stages they own.

```julia
function process(definition::AbstractDefinition, source)
    admitted = entitle(definition, source)
    selected = select(definition, admitted)
    product = build(definition, selected)
    return finish(definition, admitted, selected, product)
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
Reject unsupported definition/source pairs before partial work. Introduce a
mutable context only when several stages genuinely share buffers, resources,
or evolving state.

The sanctioned `@orchestrator` actions and their CI checks are listed in
[Grammar invariants](developers.md). Do not apply `@orchestrator` to an open
generic such as `observe`.

## Ownership-centred recursive module layout

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
detached plot page with Makie belongs in the Makie extension. A method that
parses an external file belongs with the format owner.

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

An owned result exposes numerical meaning through `observe`. A consumer asks
for explicit requests through `observables` before tabulation or drawing:

```text
scientific request
→ observe / observables
→ detached values + Quantity + UnitExpr
→ DataFrame, report, or plot
```

Apply these rules:

1. Construct scientific requests with `@observe` and retain their tuple form.
2. Define quantity identity, units, labels, and symbols in `Units`.
3. Let the result owner implement `observe` and declare supported requests.
4. Publish once before a table or plot consumes the values.
5. Keep scientific tables wide: coordinates identify rows and each observed
   quantity owns one column.
6. Select a drawing primitive with the Makie function itself.
7. Create managed axes through `axis!`, or attach an ordinary Makie axis with
   `register!`.
8. Add preview geometry beside the DataModel type that owns the cable part.

Do not inspect UQ storage fields from plotting or reporting code. `samples`,
`statistics`, and `histograms` select stored products; the scientific selector
still identifies the physical quantity.

## Docstrings

Docstrings use DocStringExtensions abbreviations so declarations remain aligned
with the implementation.

### Placement and content

- Place a docstring immediately before the documented module, type,
  constructor, function, or constant.
- Use triple double quotes, except for concise field and constant docstrings.
- Describe implemented behaviour. Do not infer equations, units, or defaults
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

Versions follow [Semantic Versioning](https://semver.org/). Public behaviour
remains compatible within a minor release. A deprecation includes a migration
path before removal.

Commit subjects use scoped Conventional Commits, begin with a lowercase
description, and stay within 72 characters:

```text
fix(engine): reject unsupported formulation options
```

Every change includes tests at the closest relevant scope. Core tests do not
load optional packages. CairoMakie and other optional paths run in their own
environments. Public examples should be executable and self-contained.
