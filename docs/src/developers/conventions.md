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

## Docstrings

Docstrings use DocStringExtensions abbreviations so declarations remain aligned
with the implementation.

### Placement and content

- Place a docstring immediately before the documented module, type,
  constructor, function, or constant.
- Use triple double quotes, except for concise field and constant docstrings.
- Describe implemented behavior. Do not infer equations, units, or defaults
  from a name.
- Describe the current API directly. Keep development history, removed APIs,
  and migration notes out of reference documentation and docstrings.
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

Use this section order, omitting sections that do not add information:

1. description and any implemented equation.
2. `# Arguments`.
3. `# Keywords`.
4. `# Returns`.
5. `# Notes` for assumptions or limitations.
6. `# Errors` for deliberate exceptions.
7. `# Examples`.

````julia
"""
$(TYPEDSIGNATURES)

Describe the implemented operation.

# Arguments

- `value`: Physical input `\\[unit\\]`.

# Keywords

- `length_prefix`: metric prefix of a per-length denominator. Default: `:kilo`.

# Returns

- Completed value in `\\[unit\\]`.

# Errors

- Throws `ArgumentError` when `length_prefix` is unsupported.

# Examples

```jldoctest
result = $(FUNCTIONNAME)(1.0; length_prefix=:kilo) # [unit]
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
then uses `$(IMPORTS)` and `$(EXPORTS)` when those lists aid the reader. For a physical constant, use a concise single-line docstring with its symbol and SI
unit.

## Prose checks

Vale and cspell, configured in `.github/prose/`, pass on every maintained file before a step
reports. A finding is fixed by rewriting the sentence, never by altering the spelling of a
word or by adding a disable directive. The dictionary gains a word, and a rule changes, only
by a decision of the maintainers.

