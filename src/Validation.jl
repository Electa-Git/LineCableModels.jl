"""
	LineCableModels.Validation

The [`Validation`](@ref) module implements a trait‑driven, two‑stage input checking pipeline for component constructors in `LineCableModels`. Inputs are first *sanitized* (arity and shape checks on raw arguments), then *parsed* (proxy values normalized to numeric radii), and finally validated by a generated set of rules.

# Overview

- Centralized constructor input handling: `sanitize` → `parse` → rule application.
- Trait hooks configure per‑type behavior (`has_radii`, `has_temperature`, `required_fields`, `keyword_fields`, `coercive_fields`, etc.).
- Rules are small value objects (`Rule` subtypes) applied to a normalized `NamedTuple`.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Validation

# Export public API
export validate!, has_radii, has_temperature, extra_rules,
    sanitize, parse, is_radius_input, required_fields, keyword_fields, coercive_fields,
    Finite, Nonneg, Positive, IntegerField, Less, LessEq, IsA, Normalized, OneOf

# Load common dependencies
using ..LineCableModels
using ..Utils

include("utils/commondeps.jl")

"""
$(TYPEDEF)

Base abstract type for validation rules. All concrete rule types must subtype [`Rule`](@ref) and provide an `_apply(::Rule, nt, ::Type{T})` method that checks a field in the normalized `NamedTuple` `nt` for the component type `T`.

$(TYPEDFIELDS)
"""
abstract type Rule end

"""
$(TYPEDEF)

Rule that enforces finiteness of a numeric field.

$(TYPEDFIELDS)
"""
struct Finite <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a field to be non‑negative (`≥ 0`).

$(TYPEDFIELDS)
"""
struct Nonneg <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a field to be strictly positive (`> 0`).

$(TYPEDFIELDS)
"""
struct Positive <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a field to be of an integer type.

$(TYPEDFIELDS)
"""
struct IntegerField <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a strict ordering constraint `a < b` between two fields.

$(TYPEDFIELDS)
"""
struct Less <: Rule
    "Left‑hand field name."
    a::Symbol
    "Right‑hand field name."
    b::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a non‑strict ordering constraint `a ≤ b` between two fields.

$(TYPEDFIELDS)
"""
struct LessEq <: Rule
    "Left‑hand field name."
    a::Symbol
    "Right‑hand field name."
    b::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a field to be `isa M` for a specified type parameter `M`.

$(TYPEDFIELDS)
"""
struct IsA{M} <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces that a field has already been normalized to a numeric value during parsing. Intended to guard that `parse` has executed and removed proxies.

$(TYPEDFIELDS)
"""
struct Normalized <: Rule
    "Name of the field to check."
    name::Symbol
end

"""
$(TYPEDEF)

Rule that enforces a field to be `in` the set `S`.

$(TYPEDFIELDS)
"""
struct OneOf{S} <: Validation.Rule
    name::Symbol
    set::S
end


"""
$(TYPEDSIGNATURES)

Returns the simple (unqualified) name of type `T` as a `String`. Utility for constructing diagnostic messages.

# Arguments

- `::Type{T}`: Type whose name is requested \\[dimensionless\\].

# Returns

- `String` with the type name \\[dimensionless\\].

# Examples

```julia
name = $(FUNCTIONNAME)(Float64)  # "Float64"
```
"""
@inline _typename(::Type{T}) where {T} = String(nameof(T))

"""
$(TYPEDSIGNATURES)

Returns a compact textual representation of `x` for error messages.

# Arguments

- `x`: Value to represent \\[dimensionless\\].

# Returns

- `String` with a compact `repr` \\[dimensionless\\].

# Examples

```julia
s = $(FUNCTIONNAME)(:field)  # ":field"
```
"""
@inline _repr(x) = repr(x; context=:compact => true)

"""
$(TYPEDSIGNATURES)

Asserts that `x` is a real (non‑complex) number. Used by rule implementations before performing numeric comparisons.

# Arguments

- `field`: Field name used in diagnostics \\[dimensionless\\].
- `x`: Value to check \\[dimensionless\\].
- `::Type{T}`: Component type for contextualized messages \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.

# Errors

- `ArgumentError` if `x` is not `isa Number` or is a `Complex` value.

# Examples

```julia
$(FUNCTIONNAME)(:radius_in, 0.01, SomeType)  # ok
```
"""
@inline function _ensure_real(field::Symbol, x, ::Type{T}) where {T}
    if !(x isa Number) || x isa Complex
        throw(ArgumentError("[$(_typename(T))] $field must be a real number, got $(typeof(x)): $(_repr(x))"))
    end
end

"""
$(TYPEDSIGNATURES)

Applies [`Finite`](@ref) to ensure the target field is a finite real number.

# Arguments

- `r`: Rule instance \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Finite, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    _ensure_real(r.name, x, T)
    isfinite(x) || throw(DomainError("[$(_typename(T))] $(r.name) must be finite, got $x"))
end

"""
$(TYPEDSIGNATURES)

Applies [`Nonneg`](@ref) to ensure the target field is `≥ 0`.

# Arguments

- `r`: Rule instance \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Nonneg, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    _ensure_real(r.name, x, T)
    x >= 0 || throw(ArgumentError("[$(_typename(T))] $(r.name) must be ≥ 0, got $x"))
end

"""
$(TYPEDSIGNATURES)

Applies [`Positive`](@ref) to ensure the target field is `> 0`.

# Arguments

- `r`: Rule instance \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Positive, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    _ensure_real(r.name, x, T)
    x > 0 || throw(ArgumentError("[$(_typename(T))] $(r.name) must be > 0, got $x"))
end

"""
$(TYPEDSIGNATURES)

Applies [`IntegerField`](@ref) to ensure the target field is an `Integer`.

# Arguments

- `r`: Rule instance \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::IntegerField, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    x isa Integer || throw(ArgumentError("[$(_typename(T))] $(r.name) must be Integer, got $(typeof(x))"))
end

"""
$(TYPEDSIGNATURES)

Applies [`Less`](@ref) to ensure `nt[a] < nt[b]`.

# Arguments

- `r`: Rule instance with fields `a` and `b` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Less, nt, ::Type{T}) where {T}
    a = getfield(nt, r.a)
    b = getfield(nt, r.b)
    _ensure_real(r.a, a, T)
    _ensure_real(r.b, b, T)
    a < b || throw(ArgumentError("[$(_typename(T))] $(r.a) < $(r.b) violated (got $a ≥ $b)"))
end

"""
$(TYPEDSIGNATURES)

Applies [`LessEq`](@ref) to ensure `nt[a] ≤ nt[b]`.

# Arguments

- `r`: Rule instance with fields `a` and `b` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::LessEq, nt, ::Type{T}) where {T}
    a = getfield(nt, r.a)
    b = getfield(nt, r.b)
    _ensure_real(r.a, a, T)
    _ensure_real(r.b, b, T)
    a <= b || throw(ArgumentError("[$(_typename(T))] $(r.a) ≤ $(r.b) violated (got $a > $b)"))
end

"""
$(TYPEDSIGNATURES)

Applies [`IsA{M}`](@ref) to ensure a field is of type `M`.

# Arguments

- `r`: Rule instance parameterized by `M` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::IsA{M}, nt, ::Type{T}) where {T,M}
    x = getfield(nt, r.name)
    x isa M || throw(ArgumentError("[$(_typename(T))] $(r.name) must be $(M), got $(typeof(x))"))
end

"""
$(TYPEDSIGNATURES)

Applies [`Normalized`](@ref) to ensure the field has been converted to a numeric value during parsing.

# Arguments

- `r`: Rule instance \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Normalized, nt, ::Type{T}) where {T}
    x = getfield(nt, r.name)
    x isa Number || throw(ArgumentError("[$(_typename(T))] $(r.name) must be normalized Number; got $(typeof(x))"))
end

"""
$(TYPEDSIGNATURES)

Applies [`OneOf`](@ref) to ensure the target field is contained in a specified set.

# Arguments

- `r`: Rule instance with fields `name` and `set` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::OneOf{S}, nt, ::Type{T}) where {S,T}
    x = getfield(nt, r.name)
    (x in r.set) || throw(ArgumentError("[$(String(nameof(T)))] $(r.name) must be one of $(collect(r.set)); got $(x)"))
end

"""
$(TYPEDSIGNATURES)

Trait hook enabling the annular radii rule bundle on fields `:radius_in` and `:radius_ext` (normalized numbers required, finiteness, non‑negativity, and the ordering constraint `:radius_in` < `:radius_ext`). It does **not** indicate the mere existence of radii; it opts in to the annular/coaxial shell geometry checks.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- `Bool` flag.

# Examples

```julia
Validation.has_radii(Tubular)
```
"""
has_radii(::Type) = false       # Default = false/empty. Components extend these.

"""
$(TYPEDSIGNATURES)

Trait hook indicating whether a component type uses a `:temperature` field subject to finiteness checks.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- `Bool` flag.
"""
has_temperature(::Type) = false

"""
$(TYPEDSIGNATURES)

Trait hook providing additional rule instances for a component type. Used to append per‑type constraints after the standard bundles.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of [`Rule`](@ref) instances.
"""
extra_rules(::Type) = ()          # per-type extras

"""
$(TYPEDSIGNATURES)

Trait hook listing required fields that must be present after positional→named merge in `sanitize`.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of required field names.
"""
required_fields(::Type) = ()

"""
$(TYPEDSIGNATURES)

Trait hook listing optional keyword fields considered by `sanitize`.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of optional keyword field names.
"""
keyword_fields(::Type) = ()

"""
$(TYPEDSIGNATURES)

Trait hook listing **coercive** fields: values that participate in numeric promotion and will be converted to the promoted type by the convenience constructor. Defaults to all fields (`required_fields ∪ keyword_fields`). Types may override to exclude integers or categorical fields.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of field names that are coerced.
"""
coercive_fields(::Type{T}) where {T} = (required_fields(T)..., keyword_fields(T)...)

"""
$(TYPEDSIGNATURES)

Trait predicate that defines admissible *raw* radius inputs for a component type during `sanitize`. The default accepts real, non‑complex numbers only. Component code may extend this to allow proxies (e.g., `AbstractCablePart`, `Thickness`, `Diameter`).

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `x`: Candidate value \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance.

# Examples

```julia
Validation.is_radius_input(Tubular, 0.01)  # true by default
Validation.is_radius_input(Tubular, 1 + 0im)  # false (complex)
```

# See also

- [`sanitize`](@ref)
"""
is_radius_input(::Type{T}, x) where {T} = (x isa Number) && !(x isa Complex)

"""
$(TYPEDSIGNATURES)

Field‑aware acceptance predicate used by `sanitize` to distinguish inner vs. outer radius policies. The default forwards to [`is_radius_input(::Type{T}, x)`](@ref) when no field‑specific method is defined.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{F}`: Field tag; typically `Val(:radius_in)` or `Val(:radius_ext)` \\[dimensionless\\].
- `x`: Candidate value \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance.

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_in), 0.01)   # true
Validation.is_radius_input(Tubular, Val(:radius_ext), 0.01)  # true
```

# See also

- [`sanitize`](@ref)
- [`is_radius_input(::Type{T}, x)`](@ref)
"""
is_radius_input(::Type{T}, ::Val{F}, x) where {T,F} = is_radius_input(T, x)

"""
$(TYPEDSIGNATURES)

Default policy for **inner** radius raw inputs: accept real numbers.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{:radius_in}`: Field tag for the inner radius \\[dimensionless\\].
- `x::Number`: Candidate value \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance (`true` for real, non‑complex numbers).

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_in), 0.0)   # true
Validation.is_radius_input(Tubular, Val(:radius_in), 1+0im) # false
```
"""
is_radius_input(::Type{T}, ::Val{:radius_in}, x::Number) where {T} = (x isa Number) && !(x isa Complex)
is_radius_input(::Type{T}, ::Val{:radius_in}, ::Any) where {T} = false

"""
$(TYPEDSIGNATURES)

Default policy for **outer** radius raw inputs (annular shells): accept real numbers. Proxies are rejected at this stage to prevent zero‑thickness stacking.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `::Val{:radius_ext}`: Field tag for the outer radius \\[dimensionless\\].
- `x::Number`: Candidate value \\[dimensionless\\].

# Returns

- `Bool` indicating acceptance (`true` for real, non‑complex numbers).

# Examples

```julia
Validation.is_radius_input(Tubular, Val(:radius_ext), 0.02)  # true
```
"""
is_radius_input(::Type{T}, ::Val{:radius_ext}, x::Number) where {T} = (x isa Number) && !(x isa Complex)
is_radius_input(::Type{T}, ::Val{:radius_ext}, ::Any) where {T} = false

"""
$(TYPEDSIGNATURES)

Performs raw input checks and shapes the input into a `NamedTuple` without parsing proxies. Responsibilities: arity validation, positional→named mapping, required field presence, and raw acceptance of radius inputs via `is_radius_input` when `has_radii(T)` is true.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `args::Tuple`: Positional arguments as received by the convenience constructor \\[dimensionless\\].
- `kwargs::NamedTuple`: Keyword arguments \\[dimensionless\\].

# Returns

- `NamedTuple` with raw (unparsed) fields.

# Errors

- `ArgumentError` on invalid arity, excess positional arguments, missing required fields, or rejected raw radius inputs.

# Examples

```julia
nt = $(FUNCTIONNAME)(Tubular, (0.01, 0.02, material), (; temperature = 20.0,))
```
"""
function sanitize(::Type{T}, args::Tuple, kwargs::NamedTuple) where {T}

    # soviet-style arity check 
    req = required_fields(T)
    kw = keyword_fields(T)

    nreq = length(req)
    na = length(args)

    if na != nreq
        names = join(string.(req), ", ")
        throw(ArgumentError("[$(_typename(T))] expected exactly $nreq positional args ($names); got $na. Optionals must be keywords."))
    end

    # positional -> named (exactly nreq items)
    nt_pos = (; (req[i] => args[i] for i = 1:nreq)...)

    # Reject unknown keywords 
    for k in keys(kwargs)
        if !(k in kw) && !(k in req)
            throw(ArgumentError("[$(_typename(T))] unknown keyword '$k'. Allowed keywords: $(join(string.(kw), ", "))."))
        end
    end

    # merge (keywords override same-name positional if any — rare, but consistent)
    nt = merge(nt_pos, kwargs)

    # raw acceptance for radii
    if has_radii(T)
        haskey(nt, :radius_in) || throw(ArgumentError("[$(_typename(T))] missing 'radius_in'."))
        haskey(nt, :radius_ext) || throw(ArgumentError("[$(_typename(T))] missing 'radius_ext'."))
        is_radius_input(T, Val(:radius_in), nt.radius_in) ||
            throw(ArgumentError("[$(_typename(T))] radius_in not an accepted input: $(typeof(nt.radius_in))"))
        is_radius_input(T, Val(:radius_ext), nt.radius_ext) ||
            throw(ArgumentError("[$(_typename(T))] radius_ext not an accepted input: $(typeof(nt.radius_ext))"))
    end
    return nt
end

"""
$(TYPEDSIGNATURES)

Parses and normalizes raw inputs produced by [`sanitize`](@ref) into the canonical form expected by rules. Default is identity; component code overrides this to resolve proxy radii to numeric values while preserving uncertainty semantics.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `nt::NamedTuple`: Raw inputs from `sanitize` \\[dimensionless\\].

# Returns

- `NamedTuple` with normalized fields (e.g., numeric `:radius_in`, `:radius_ext`).
"""
parse(::Type, nt) = nt

"""
$(TYPEDSIGNATURES)

Generates, at compile time, the tuple of rules to apply for component type `T`. The result concatenates standard bundles driven by traits and any rules returned by `extra_rules(T)`.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of [`Rule`](@ref) instances to apply in order.
"""
@generated function _rules(::Type{T}) where {T}
    :((
        (has_radii(T) ? (Normalized(:radius_in), Normalized(:radius_ext),
            Finite(:radius_in), Nonneg(:radius_in),
            Finite(:radius_ext), Nonneg(:radius_ext),
            Less(:radius_in, :radius_ext)) : ())...,
        (has_temperature(T) ? (Finite(:temperature),) : ())...,
        extra_rules(T)...
    ))
end

"""
$(TYPEDSIGNATURES)

Runs the full validation pipeline for a component type: `sanitize` (arity and raw checks), `parse` (proxy normalization), then application of the generated rule set. Intended to be called from convenience constructors.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `args...`: Positional arguments \\[dimensionless\\].
- `kwargs...`: Keyword arguments \\[dimensionless\\].

# Returns

- `NamedTuple` containing normalized fields ready for construction.

# Errors

- `ArgumentError` from `sanitize` or rule checks; `DomainError` for finiteness violations.

# Examples

```julia
nt = $(FUNCTIONNAME)(Tubular, 0.01, 0.02, material; temperature = 20.0)
# use nt.radius_in, nt.radius_ext, nt.temperature thereafter
```
"""
function validate!(::Type{T}, args...; kwargs...) where {T}
    # One validate! to rule them all

    nt0 = sanitize(T, args, (; kwargs...))
    nt1 = parse(T, nt0)
    # if has_radii: Normalized ensures numbers post-parse; if not numbers, rules will throw
    rules = _rules(T)
    @inbounds for i in eachindex(rules)
        _apply(rules[i], nt1, T)
    end
    return nt1
end

end # module Validation
