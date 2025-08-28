# Validation module

## Contents
```@contents
Pages = ["validation.md"]
Depth = 3
```

---

This section documents the validation framework used by component constructors. The design is deterministic, non‑magical, and trait‑driven. The flow is:

```julia
sanitize(Type, args, kwargs)  →  parse(Type, nt)  →  apply rules  →  construct
```

The **typed cores** accept **numbers only**. All proxy handling happens in the **convenience constructors** via `validate!`.

---

## Architecture

### Pipeline

* **`sanitize(::Type{T}, args::Tuple, kwargs::NamedTuple)`**

  * Rejects wrong arities: exactly `length(required_fields(T))` positionals are expected; optionals must be passed as keywords listed in `keyword_fields(T)`.
  * Maps positionals to names using `required_fields(T)`; merges keyword arguments; rejects unknown keywords.
  * If `has_radii(T) == true`, checks admissibility of raw radius inputs with `is_radius_input(T, Val(:radius_in), x)` and `is_radius_input(T, Val(:radius_ext), x)`. The default accepts only real, non‑complex numbers; types may extend to allow proxies.
  * Returns a **raw** `NamedTuple`.

* **`parse(::Type{T}, nt)`**

  * Normalizes raw inputs to canonical representation (e.g., radius proxies → numeric radii) while preserving domain semantics (e.g., uncertainty reset rules).
  * Returns a **normalized** `NamedTuple`.

* **Rule application**

  * `_rules(T)` is generated from traits; evaluates over the normalized `NamedTuple`.
  * Standard bundles are injected when traits are `true` (e.g., for radii: `Normalized`, `Finite`, `Nonneg`, `Less`).
  * Per‑type extras come from `extra_rules(T)`.

* **`validate!(::Type{T}, args...; kwargs...)`**

  * Orchestrates the pipeline. Use this in all convenience constructors.

### Traits (configuration surface)

* `has_radii(::Type{T})::Bool` — enables the radii rule bundle and raw acceptance checks.
* `has_temperature(::Type{T})::Bool` — enables finiteness check on `:temperature`.
* `required_fields(::Type{T})::NTuple` — positional keys.
* `keyword_fields(::Type{T})::NTuple` — keyword argument keys.
* `coercive_fields(::Type{T})::NTuple` — values that participate in type promotion and will be coerced (default: `required_fields ∪ keyword_fields`).
* `is_radius_input(::Type{T}, Val(:field), x)::Bool` — raw admissibility predicate for radii inputs; extend to allow proxies by field.
* `extra_rules(::Type{T})::NTuple{K,Rule}` — additional constraints appended to the generated bundle.

**Import before extending**:

```julia
import ..Validation: has_radii, has_temperature, required_fields, keyword_fields,
                     coercive_fields, is_radius_input, parse, extra_rules
```

Failing to import will create shadow functions in your module; the engine will not see your methods.

---

## Rules

Rules are small value types `struct <: Rule` with an `_apply(::Rule, nt, ::Type{T})` method. All rule methods must:

* Read data from the **normalized** `NamedTuple` `nt`.
* Throw `ArgumentError` for logical violations; `DomainError` for numerical domain violations (non‑finite).
* Avoid allocations; use `@inline` where appropriate.

### Standard rules

* `Normalized(:field)` — field must be numeric post‑parse.
* `Finite(:field)` — `isfinite` must hold.
* `Nonneg(:field)` — value `≥ 0`.
* `Positive(:field)` — value `> 0`.
* `IntegerField(:field)` — value `isa Integer`.
* `Less(:a,:b)` — strict ordering `a < b`.
* `LessEq(:a,:b)` — non‑strict ordering `a ≤ b`.
* `IsA{M}(:field)` — type membership check.
* `OneOf(:field, set)` — membership in a finite set.

### Custom rule pattern

```julia
struct InRange{T} <: Validation.Rule
    name::Symbol; lo::T; hi::T
end

@inline function Validation._apply(r::InRange, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    (x isa Number && !(x isa Complex)) || throw(ArgumentError("[$(String(nameof(T)))] $(r.name) must be real"))
    (r.lo ≤ x ≤ r.hi) || throw(ArgumentError("[$(String(nameof(T)))] $(r.name) out of range $(r.lo):$(r.hi), got $(x)"))
end
```

Attach via `extra_rules(::Type{X}) = (InRange(:alpha, 0.0, 1.0), ...)`.

---

## Example implementation — `DataModel.Tubular`

**Typed core** (numbers only):

```julia
function Tubular(radius_in::T, radius_ext::T, material_props::Material{T}, temperature::T) where {T<:REALSCALAR}
```

Rationale: all proxy resolution must happen before reaching the typed core to avoid duplicate parsing and to keep type promotion deterministic.

**Trait configuration**:

```julia
const _REQ_TUBULAR = (:radius_in, :radius_ext, :material_props,)
const _OPT_TUBULAR = (:temperature,)
const _DEFS_TUBULAR = (T₀,)

Validation.has_radii(::Type{Tubular}) = true
Validation.has_temperature(::Type{Tubular}) = true
Validation.required_fields(::Type{Tubular}) = _REQ_TUBULAR
Validation.keyword_fields(::Type{Tubular}) = _OPT_TUBULAR
```

* `has_radii = true` enables the radii bundle: `Normalized(:radius_in)`, `Normalized(:radius_ext)`, `Finite`, `Nonneg`, `Less(:radius_in,:radius_ext)`.
* `has_temperature = true` adds `Finite(:temperature)`.
* `required_fields` defines the mandatory positional fields.
* `keyword_fields` defines the optional fields that will receive default values.

**Raw proxy acceptance**:

```julia
Validation.is_radius_input(::Type{Tubular}, ::Val{:radius_in}, x::AbstractCablePart) = true
Validation.is_radius_input(::Type{Tubular}, ::Val{:radius_ext}, x::Thickness) = true
Validation.is_radius_input(::Type{Tubular}, ::Val{:radius_ext}, x::Diameter) = true
```

The inner radius may take an existing cable part (its `radius_ext`); the outer radius may take a `Thickness` or `Diameter` wrapper.

**Extra rules**:

```julia
Validation.extra_rules(::Type{Tubular}) = (IsA{Material}(:material_props),)
```

**Parsing**:

```julia
Validation.parse(::Type{Tubular}, nt) = begin
    rin, rex = _normalize_radii(Tubular, nt.radius_in, nt.radius_ext)
    (; nt..., radius_in=rin, radius_ext=rex)
end
```

**Convenience constructor**:

```julia
@_ctor Tubular _REQ_TUBULAR _OPT_TUBULAR _DEFS_TUBULAR
```

This expands to a weakly‑typed method that calls `validate!`, promotes using `_promotion_T`, coerces via `_coerced_args`, and delegates to the numeric core.

**Failure modes intentionally trapped**:

* Wrong arity, missing keys → `sanitize` via `required_fields`.
* String/complex radii → `sanitize` via default `is_radius_input` (unless explicitly allowed).
* Forgotten parsing after allowing proxies → caught by `Normalized` rules.
* Geometry violations (`radius_in ≥ radius_ext`) → `Less(:radius_in,:radius_ext)`.

---

## Template for a new component

```julia
# 1) Numeric core (numbers only)
function NewPart(a::T, b::T, material::Material{T}, temperature::T) where {T<:REALSCALAR}
    # compute derived, then construct
end

# 2) Trait config
Validation.has_radii(::Type{NewPart}) = true
Validation.has_temperature(::Type{NewPart}) = true
Validation.required_fields(::Type{NewPart}) = (:a, :b, :material)
Validation.keyword_fields(::Type{NewPart})  = (:temperature,)

# 3) Raw acceptance (extend only what you intend to parse)
Validation.is_radius_input(::Type{NewPart}, ::Val{:radius_in},  x::AbstractCablePart) = true
Validation.is_radius_input(::Type{NewPart}, ::Val{:radius_ext}, x::Thickness)        = true

# 4) Extra rules
Validation.extra_rules(::Type{NewPart}) = (
    IsA{Material}(:material),
)

# 5) Parsing (proxy → numeric)
Validation.parse(::Type{NewPart}, nt) = begin
    a′, b′ = _normalize_radii(NewPart, nt.a, nt.b)
    (; nt..., a=a′, b=b′)
end

# 6) Convenience constructor — generated
@_ctor NewPart (:a, :b, :material) (:temperature,) (T₀,)
```

---

## Extending traits

Traits are just methods. Add traits only when behavior toggles; avoid proliferation.

### New feature flags

Example: a shielding flag with its own bundle.

```julia
Validation.has_shield(::Type) = false
Validation.has_shield(::Type{SomeType}) = true
```

Extend `_rules` inside `Validation` to splice the corresponding checks when `has_shield(T)` is true.

### Field‑specific admissibility

If only `:radius_in` should accept proxies, extend the field‑tagged predicate:

```julia
Validation.is_radius_input(::Type{X}, ::Val{:radius_in},  p::AbstractCablePart) = true
Validation.is_radius_input(::Type{X}, ::Val{:radius_ext}, ::AbstractCablePart)  = false
```

---

## Testing guidelines

* Test arity: `()`, `(1)`, `(1,2)` → `ArgumentError`.
* Test raw type rejections: strings, complex numbers.
* Test proxy acceptance: prior layer objects, `Thickness`, `Diameter` when allowed.
* Test parse correctness: outputs numeric and respect uncertainty rules.
* Test rule violations: negative radii, inverted radii, non‑finite values, invalid sets.
* Test constructor round‑trip: convenience path and numeric core produce equivalent instances after coercion.

## Usage notes

* Keep all proxy handling in `parse`. Do not call normalizers in constructors.
* Error messages must be terse and contextualized with the component type name.
* Prefer tuple returns and `NamedTuple` updates to avoid allocations.
* When adding rules, benchmark `_apply` implementations.

---

## API reference

```@autodocs
Modules = [LineCableModels.Validation]
Order = [:module, :constant, :type, :function, :macro]
Public = true
Private = true
```

---

## Index
```@index
Pages   = ["validation.md"]
Order   = [:module, :constant, :type, :function, :macro]
```
