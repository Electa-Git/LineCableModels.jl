"""
	LineCableModels.Validation

The [`Validation`](@ref) module implements a trait-driven, three-phase input checking pipeline for component constructors in `LineCableModels`. Inputs are first *sanitized* (arity and shape checks on raw arguments), then *parsed* (proxy values normalized to numeric radii), and finally validated by a generated set of rules.

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
	sanitize, parse, is_radius_input, required_fields, keyword_fields, keyword_defaults,
	coercive_fields, Finite, Nonneg, Positive, IntegerField, Less, LessEq, IsA, Normalized,
	OneOf, GreaterEq, Greater, PhysicalFillLimit, Satisfies

# Module-specific dependencies
using ..Commons

include("rules.jl")
include("applyrules.jl")

"""
$(TYPEDSIGNATURES)

Trait hook enabling the annular radii rule bundle on fields `:radius_in` and `:radius_ext` (normalized numbers required, finiteness, non‑negativity, and the ordering constraint `:radius_in` < `:radius_ext`). It does **not** indicate the mere existence of radii; it opts in to the annular/coaxial shell geometry checks.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- `Bool` flag.

# Examples

```julia
Validation.has_radii(Tubular)  # true
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

Trait hook supplying **default values** for optional keyword fields.

Return either:
- a `NamedTuple` mapping defaults by name (e.g., `(temperature = T₀,)`), or
- a plain `Tuple` of defaults aligned with `keyword_fields(T)` (same order + length).

Defaults are applied in `sanitize` **after** positional→named merge and before rule application.
User-provided keywords always override these defaults.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- `NamedTuple` or `Tuple` with default keyword values.
"""
keyword_defaults(::Type) = ()

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
is_radius_input(::Type{T}, ::Val{F}, x) where {T, F} = is_radius_input(T, x)

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
is_radius_input(::Type{T}, ::Val{:radius_in}, x::Number) where {T} =
	(x isa Number) && !(x isa Complex)
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
is_radius_input(::Type{T}, ::Val{:radius_ext}, x::Number) where {T} =
	(x isa Number) && !(x isa Complex)
is_radius_input(::Type{T}, ::Val{:radius_ext}, ::Any) where {T} = false

"""
$(TYPEDSIGNATURES)

Internal helper that canonicalizes `keyword_defaults(T)` into a `NamedTuple`
keyed by `keyword_fields(T)`. Accepts:

- `()` → returns an empty `NamedTuple()`.
- `NamedTuple` → returned unchanged.
- `Tuple` → zipped **by position** with `keyword_fields(T)`; lengths must match.

This function does **not** merge user-provided keywords; callers should perform
`merge(_kwdefaults_nt(T), kwargs)` so that user values take precedence.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- `NamedTuple` mapping optional keyword field names to their default values.

# Errors

- If `keyword_defaults(T)` returns a `Tuple` whose length differs from
  `length(keyword_fields(T))`.
- If `keyword_defaults(T)` returns a value that is neither `()` nor a
  `NamedTuple` nor a `Tuple`.

# Examples

```julia
# Suppose:
#   keyword_fields(::Type{X})   = (:temperature, :lay_direction)
#   keyword_defaults(::Type{X}) = (T₀, 1)

$(FUNCTIONNAME)(X)  # => (temperature = T₀, lay_direction = 1)

# If defaults are already a NamedTuple:
#   keyword_defaults(::Type{Y}) = (temperature = 25.0,)
$(FUNCTIONNAME)(Y)  # => (temperature = 25.0,)
````

# See also

* [`keyword_fields`](@ref)
* [`keyword_defaults`](@ref)
* [`sanitize`](@ref)
  """
@inline function _kwdefaults_nt(::Type{T}) where {T}
	defs = keyword_defaults(T)
	defs === () && return NamedTuple()
	if defs isa NamedTuple
		return defs
	elseif defs isa Tuple
		keys = keyword_fields(T)
		length(keys) == length(defs) ||
			Base.error(
				"[$(String(nameof(T)))] keyword_defaults length $(length(defs)) ≠ keyword_fields length $(length(keys))",
			)
		return NamedTuple{keys}(defs)
	else
		Base.error(
			"[$(String(nameof(T)))] keyword_defaults must be NamedTuple or Tuple; got $(typeof(defs))",
		)
	end
end

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

# Notes

Required arguments must be positional; optional arguments must be passed as keywords. Positional arity must equal `length(required_fields(T))`.

# Examples

```julia
nt = $(FUNCTIONNAME)(Tubular, (0.01, 0.02, material), (; temperature = 20.0,))
```
"""
function sanitize(::Type{T}, args::Tuple, kwargs::NamedTuple) where {T}
	# -- hard arity on required positionals --
	req = required_fields(T)
	kw = keyword_fields(T)
	nreq = length(req)
	na = length(args)
	if na != nreq
		names = join(string.(req), ", ")
		throw(
			ArgumentError(
				"[$(_typename(T))] expected exactly $nreq positional args ($names); got $na. Optionals must be keywords.",
			),
		)
	end

	# positional -> named
	nt_pos = (; (req[i] => args[i] for i ∈ 1:nreq)...)

	# reject unknown keywords (strict)
	for k in keys(kwargs)
		if !(k in kw) && !(k in req)
			throw(
				ArgumentError(
					"[$(_typename(T))] unknown keyword '$k'. Allowed keywords: $(join(string.(kw), ", ")).",
				),
			)
		end
	end

	# user kw override positionals (if any same names)
	nt = merge(nt_pos, kwargs)

	# backfill missing optional keywords with trait defaults ---
	# defaults first, then user-provided values win
	nt = merge(_kwdefaults_nt(T), nt)

	# radii raw acceptance (unchanged)
	if has_radii(T)
		haskey(nt, :radius_in) ||
			throw(ArgumentError("[$(_typename(T))] missing 'radius_in'."))
		haskey(nt, :radius_ext) ||
			throw(ArgumentError("[$(_typename(T))] missing 'radius_ext'."))
		is_radius_input(T, Val(:radius_in), nt.radius_in) ||
			throw(
				ArgumentError(
					"[$(_typename(T))] radius_in not an accepted input: $(typeof(nt.radius_in))",
				),
			)
		is_radius_input(T, Val(:radius_ext), nt.radius_ext) ||
			throw(
				ArgumentError(
					"[$(_typename(T))] radius_ext not an accepted input: $(typeof(nt.radius_ext))",
				),
			)
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

Generates (at compile time, via a `@generated` function) the tuple of rules to apply for component type `T`. The result concatenates standard bundles driven by traits and any rules returned by `extra_rules(T)`.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Tuple of [`Rule`](@ref) instances to apply in order.
"""
@generated function _rules(::Type{T}) where {T}
	:((
		(
			has_radii(T) ?
			(Normalized(:radius_in), Normalized(:radius_ext),
				Finite(:radius_in), Nonneg(:radius_in),
				Finite(:radius_ext), Nonneg(:radius_ext),
				Less(:radius_in, :radius_ext)) : ()
		)...,
		(has_temperature(T) ? (Finite(:temperature),) : ())...,
		extra_rules(T)...,
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

# See also
- [`sanitize`](@ref)
- [`parse`](@ref)
- [`coercive_fields`](@ref)
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
