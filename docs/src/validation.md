# Validation module

## Contents

```@contents
Pages = ["validation.md"]
Depth = 3
```

The strict materialized constructors use a trait-driven validation sequence:

```julia
sanitize(Type, args, kwargs) → parse(Type, values) → apply rules → construct
```

This layer accepts numeric physical state. It does not interpret tuples,
previous/next layer objects, or radial wrapper types. Parameter variation belongs
to `Grid` and declarative construction belongs to `Gridspace` builders.

## Radial inputs

The typed physical constructors accept numeric inner and outer radii:

```julia
Tubular(r_in, r_ex, material)
```

Their narrow named convenience accepts one forward declaration:

```julia
Tubular(r_in, material; radius = r_ex)
Tubular(r_in, material; thickness = delta_r)
```

Exactly one of `radius` and `thickness` is required. Group `add!` methods retain
the useful stacking operation by supplying the group's current numeric outer
radius as the new layer's inner radius:

```julia
add!(group, Tubular, material; thickness = delta_r)
```

No layer object is accepted as a radius proxy by the physical constructors.

## Validation sequence

- `sanitize` checks positional arity, accepted keywords, defaults, and real
  numeric radii.
- `parse` normalizes component-specific fields. It is the identity by
  default.
- `_rules(T)` assembles the standard trait-driven rules and `extra_rules(T)`.
- `validate!` runs the complete sequence used by convenience constructors.

The principal traits are:

- `has_radii(::Type)` for annular radius rules;
- `has_temperature(::Type)` for temperature finiteness;
- `required_fields`, `keyword_fields`, and `keyword_defaults` for constructor
  shape;
- `coercive_fields` for numeric promotion;
- `extra_rules` for component-specific constraints.

## Rules

Rules are small `Rule` values applied to the normalized `NamedTuple`. Standard
rules include `Normalized`, `Finite`, `Nonneg`, `Positive`, `IntegerField`,
`Less`, `LessEq`, `IsA`, and `OneOf`. Logical violations throw
`ArgumentError`; non-finite numerical values throw `DomainError`.

For an annular component, the core trait setup is:

```julia
Validation.has_radii(::Type{Tubular}) = true
Validation.has_temperature(::Type{Tubular}) = true
Validation.required_fields(::Type{Tubular}) = (:r_in, :r_ex, :material_props)
Validation.keyword_fields(::Type{Tubular}) = (:temperature,)
Validation.keyword_defaults(::Type{Tubular}) = (T₀,)
Validation.extra_rules(::Type{Tubular}) = (IsA{Material}(:material_props),)
```

`has_radii` adds numeric normalization, finiteness, non-negativity, and the
strict `r_in < r_ex` ordering check.

## Component checklist

When adding a physical part:

1. implement the fully typed numeric core;
2. declare required, optional, default, and coercive fields;
3. opt into standard geometry and temperature bundles;
4. add only the component-specific rules;
5. normalize fields in `parse` when necessary;
6. route the convenience constructor through `validate!`.

Tests should cover arity, rejected strings/complex radii, non-finite and
inverted geometry, type promotion, and equivalence between the convenience and
typed-core construction paths.

## API reference

```@autodocs
Modules = [LineCableModels.Validation]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

## Index

```@index
Pages = ["validation.md"]
Order = [:module, :constant, :type, :function, :macro]
```
