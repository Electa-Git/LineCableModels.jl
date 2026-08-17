# Docstrings

The following docstring standards apply throughout LineCableModels. Docstrings
use
[`DocStringExtensions.jl`](https://github.com/JuliaDocs/DocStringExtensions.jl)
abbreviations by default so signatures, type declarations, fields, and module
contents remain synchronized with the implementation.

## General principles

1. **Placement:** Place a docstring immediately before the struct,
   constructor, function, module, or constant that it describes.
2. **Delimiter:** Use triple double quotes (`"""`) except for individual struct
   fields and constants, which use single-line string docstrings.
3. **Implementation first:** Describe behavior implemented by the code. Do not
   infer behavior, units, or mathematics from a name alone.
4. **Generated structure:** Prefer the appropriate `DocStringExtensions`
   abbreviation over a handwritten signature, type declaration, field list, or
   module inventory.
5. **Conciseness:** Present each fact once, in the appropriate section.
6. **Tone:** Use precise scientific language. Avoid contractions,
   colloquialisms, and ambiguous wording.

Dedicated cross-reference sections are not maintained. Use an inline link such
as [`LineCableModels.DataModel.CableDesign`](@ref) where a relationship
materially helps the explanation.

## Physical unit formatting

Every argument, return value, struct field, and constant representing a
physical quantity states its SI unit.

1. Enclose units in double-backslash escaped square brackets, for example
   `\\[m\\]`, `\\[Hz\\]`, `\\[Ω\\]`, and `\\[H/m\\]`.
2. Mark a dimensionless physical quantity as `\\[dimensionless\\]`.
3. Do not add units to non-physical counters, flags, indices, symbols, or
   collection sizes.
4. Use standard SI symbols and the Unicode middle dot for multiplication, for
   example `\\[Ω·m\\]`.
5. Inside comments in Julia example blocks, use ordinary square brackets, such
   as `# [m]`; the comment is Julia source rather than docstring prose.
6. State the basis of distributed line quantities when it matters: for
   example, `\\[Ω/m\\]` for `:per_length` and `\\[Ω\\]` for `:total`.

## Mathematical formulation

Use LaTeX whenever an implementation directly evaluates a mathematical
expression, reduction, approximation, or physical law that matters to its
meaning. Put the expression in a `math` block within the method description, followed by the definitions of symbols that are not already unambiguous from `# Arguments`.

This requirement follows the implementation, not a function-name convention.
Simple accessors, wrappers, dispatch helpers, and bookkeeping functions generally do not need a
mathematical section.

Within a Julia docstring, escape LaTeX commands with a second backslash:

````julia
"""
$(TYPEDSIGNATURES)

Return the series impedance at a specified frequency, evaluated as:

```math
Z(f) = R + \\mathrm{j} 2 \\pi f L,
```

where ``R`` is the resistance, ``L`` is the inductance, and ``f`` is the
frequency.

# Arguments

- `resistance`: Series resistance `\\[Ω\\]`.
- `inductance`: Series inductance `\\[H\\]`.
- `frequency`: Frequency `\\[Hz\\]`.

# Returns

- Complex series impedance `\\[Ω\\]`.

# Examples

```jldoctest
z = $(FUNCTIONNAME)(0.1, 1e-3, 50.0) # [Ω]
@assert z == 0.1 + im * 2π * 50.0 * 1e-3
# output
```
"""
function series_impedance_at(resistance, inductance, frequency)
````

Preserve every correct existing mathematical description. If the
implementation and its documentation disagree, inspect and test the
implementation before changing either.

## `DocStringExtensions` abbreviations

The abbreviations are Julia string interpolations evaluated when the docstring
is attached. They are the preferred codebase-wide mechanism, not an optional
legacy style.

- `$(SIGNATURES)` inserts method signatures without argument type annotations.
  Use it only when the typed form would obscure the supported public call.
- `$(TYPEDSIGNATURES)` inserts method signatures with argument types. Use it at
  the start of function and constructor docstrings by default.
- `$(FUNCTIONNAME)` inserts the documented function name. Use it in examples
  that should survive a function rename.
- `$(TYPEDEF)` inserts a type declaration, including type parameters and its
  supertype.
- `$(TYPEDFIELDS)` inserts fields, declared types, and their field docstrings.
- `$(FIELDS)` inserts fields and their field docstrings without declared types.
  Prefer `$(TYPEDFIELDS)` unless hiding implementation types improves the public
  documentation.
- `$(METHODLIST)` inserts the package's sanitized method listing and source
  locations. Reserve it for documentation whose purpose is to enumerate a
  multi-method interface; it does not replace `$(TYPEDSIGNATURES)` in an
  ordinary method docstring.
- `$(IMPORTS)` inserts modules imported by a documented module.
- `$(EXPORTS)` inserts names exported by a documented module.

Use only abbreviations that add useful rendered information. Never repeat their
generated content manually. Each package module imports the required
`DocStringExtensions` abbreviations explicitly. The package's root-level
`METHODLIST` formatter preserves concise source locations while hiding CI runner and
workspace prefixes.

## Documentation templates

These templates are intentionally complete. Remove an optional section only
when it adds no information; do not replace generated declarations with
handwritten copies.

### Structs

Use `$(TYPEDEF)` for the declaration and `$(TYPEDFIELDS)` for the field list.
Place a single-line string docstring immediately above every field:

````julia
"""
$(TYPEDEF)

Represent a section of cable material.

$(TYPEDFIELDS)
"""
struct CableSection{T <: Real}
    "Section thickness `\\[m\\]`."
    thickness::T

    "Electrical resistivity `\\[Ω·m\\]`."
    resistivity::T

    "Relative permeability `\\[dimensionless\\]`."
    relative_permeability::T
end
````

Do not use block docstrings, inline comments, or a manually maintained field
list for field documentation.

### Outer constructors

Julia attaches an ordinary docstring to an outer constructor defined at module
scope. Do not add `@doc` to it:

````julia
"""
$(TYPEDSIGNATURES)

Construct a [`CableSection`](@ref).

# Arguments

- `thickness`: Section thickness `\\[m\\]`.
- `resistivity`: Electrical resistivity `\\[Ω·m\\]`.

# Keywords

- `relative_permeability`: Relative permeability `\\[dimensionless\\]`.
  Default: `1.0`.

# Returns

- A [`CableSection`](@ref) with a concrete shared numeric type.

# Errors

- Throws `DomainError` when `thickness` is not positive.

# Examples

```jldoctest
section = $(FUNCTIONNAME)(1e-3, 1.72e-8)
@assert section.thickness == 1e-3
# output
```
"""
function CableSection(
    thickness,
    resistivity;
    relative_permeability=1.0,
)
````

### Inner constructors

An inner constructor is a method written inside the `struct` body. Julia does
not attach a bare string there as constructor documentation, so an inner
constructor requires `@doc`:

````julia
"""
$(TYPEDEF)

Represent a cable section with validated dimensions.

$(TYPEDFIELDS)
"""
struct CableSection{T <: Real}
    "Section thickness `\\[m\\]`."
    thickness::T

    "Electrical resistivity `\\[Ω·m\\]`."
    resistivity::T

    @doc """
    $(TYPEDSIGNATURES)

    Construct a validated [`CableSection`](@ref).

    # Arguments

    - `thickness`: Section thickness `\\[m\\]`.
    - `resistivity`: Electrical resistivity `\\[Ω·m\\]`.

    # Returns

    - A validated [`CableSection`](@ref).

    # Errors

    - Throws `DomainError` when `thickness` is not positive.

    # Examples

    ```jldoctest
    section = $(FUNCTIONNAME)(1e-3, 1.72e-8)
    @assert section.thickness > 0
    # output
    ```
    """
    function CableSection(thickness::T, resistivity::T) where {T <: Real}
        thickness > zero(T) || throw(DomainError(thickness))
        return new{T}(thickness, resistivity)
    end
end
````

The distinction is structural: use `@doc` for a method inside a `struct`, and
an ordinary immediately preceding docstring for a method outside it.

### Functions and methods

Start with `$(TYPEDSIGNATURES)` and use the following section order:

1. Description, without a heading, including the mathematical formulation when applicable.
2. `# Arguments`.
3. `# Keywords`, when keyword arguments need documentation.
4. `# Returns`.
5. `# Notes`, when the implementation requires assumptions, limitations, or
   additional mathematical explanation.
6. `# Errors`, for deliberate exceptions callers should anticipate.
7. `# Examples`.

Separate the description and every section with exactly one blank line. List
arguments in declaration order, state keyword defaults, and document every
member of a returned tuple individually.

````julia
"""
$(TYPEDSIGNATURES)

Describe the function's implemented purpose concisely.

# Arguments

- `arg1`: Physical input `\\[unit\\]`.
- `arg2`: Dimensionless model parameter `\\[dimensionless\\]`.

# Keywords

- `basis`: Storage basis. Supported values are `:per_length` and `:total`.
  Default: `:per_length`.

# Returns

- Description of the return value and its physical unit `\\[unit\\]`.

# Notes

Include only implementation-relevant assumptions, limitations, or
mathematical explanations.

# Errors

- Throws `ArgumentError` when `basis` is unsupported.

# Examples

```jldoctest
result = $(FUNCTIONNAME)(1.0, 0.5; basis=:per_length) # [unit]
@assert isfinite(result)
# output
```
"""
function function_name(arg1, arg2; basis=:per_length)
````

Examples use meaningful values and exercise supported public syntax. Include
an expected result only when it can be stated accurately.

### Modules

Begin with the module name indented by four spaces. Describe its purpose and
main capabilities explicitly. `$(IMPORTS)` and `$(EXPORTS)` keep the inventory
synchronized with the module definition:

````julia
"""
    ModuleName

Describe the module's purpose within LineCableModels.

# Overview

- Summarize the principal capability.
- Summarize another principal capability.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module ModuleName
````

### Constants

Use a single-line docstring. Include the conventional symbol and SI unit when
the constant represents a physical quantity:

```julia
"Magnetic constant (vacuum permeability), μ₀ = 4π × 10⁻⁷ `\\[H/m\\]`."
const VACUUM_PERMEABILITY = 4π * 1e-7
```

## Executable examples

Prefer `jldoctest` for a self-contained public example. When textual output is
not the behavior under test, use assertions before the `# output` separator and
leave the expected output empty. This keeps the syntax executable without
making display formatting part of the public behavior.

An example that requires an external executable, repository fixture, network
access, or graphical interaction belongs in integration or developer
documentation instead of a docstring. Remove stale, placeholder, or
fixture-heavy examples rather than publishing syntax that cannot be verified.

Use `$(FUNCTIONNAME)` when the example invokes the documented function. Use an
explicit different function name when that distinction helps copied code remain
readable.

## Common mistakes to avoid

- Handwriting a signature, type declaration, field list, import list, or export
  list that an abbreviation can generate.
- Omitting a central formula that materially explains the implementation.
- Inventing physical meaning, symbols, units, defaults, or numerical results.
- Omitting `@doc` from an inner constructor defined inside a `struct`.
- Adding unnecessary `@doc` to an outer constructor defined at module scope.
- Using a block docstring or inline comment for a struct field.
- Using the wrong section order or retaining empty optional sections.
- Adding a dedicated cross-reference section.
- Escaping units inside Julia example comments, or failing to escape units and
  LaTeX commands in docstring prose.
- Repeating information already generated by a `DocStringExtensions`
  abbreviation.
- Writing an example that cannot run in the documentation environment.
