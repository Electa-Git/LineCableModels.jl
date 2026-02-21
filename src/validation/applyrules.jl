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
@inline _repr(x) = repr(x; context = :compact => true)

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
		throw(
			ArgumentError(
				"[$(_typename(T))] $field must be a real number, got $(typeof(x)): $(_repr(x))",
			),
		)
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
	x isa Integer || throw(
		ArgumentError("[$(_typename(T))] $(r.name) must be Integer, got $(typeof(x))"),
	)
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
	a < b ||
		throw(ArgumentError("[$(_typename(T))] $(r.a) < $(r.b) violated (got $a ≥ $b)"))
end

"""
$(TYPEDSIGNATURES)

Applies [`Greater`](@ref) to ensure `nt[a] > nt[b]`.

# Arguments

- `r`: Rule instance with fields `a` and `b` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::Greater, nt, ::Type{T}) where {T}
	a = getfield(nt, r.a)
	b = getfield(nt, r.b)
	_ensure_real(r.a, a, T)
	_ensure_real(r.b, b, T)
	a > b ||
		throw(ArgumentError("[$(_typename(T))] $(r.a) > $(r.b) violated (got $a ≤ $b)"))
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
	a <= b ||
		throw(ArgumentError("[$(_typename(T))] $(r.a) ≤ $(r.b) violated (got $a > $b)"))
end

"""
$(TYPEDSIGNATURES)

Applies [`GreaterEq`](@ref) to ensure `nt[a] ≥ nt[b]`.

# Arguments

- `r`: Rule instance with fields `a` and `b` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws on failure.
"""
@inline function _apply(r::GreaterEq, nt, ::Type{T}) where {T}
	a = getfield(nt, r.a)
	b = getfield(nt, r.b)
	_ensure_real(r.a, a, T)
	_ensure_real(r.b, b, T)
	a >= b ||
		throw(ArgumentError("[$(_typename(T))] $(r.a) ≥ $(r.b) violated (got $a < $b)"))
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
@inline function _apply(r::IsA{M}, nt, ::Type{T}) where {T, M}
	x = getfield(nt, r.name)
	x isa M ||
		throw(ArgumentError("[$(_typename(T))] $(r.name) must be $(M), got $(typeof(x))"))
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
	x isa Number || throw(
		ArgumentError(
			"[$(_typename(T))] $(r.name) must be normalized Number; got $(typeof(x))",
		),
	)
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
@inline function _apply(r::OneOf{S}, nt, ::Type{T}) where {S, T}
	x = getfield(nt, r.name)
	(x in r.set) || throw(
		ArgumentError(
			"[$(String(nameof(T)))] $(r.name) must be one of $(collect(r.set)); got $(x)",
		),
	)
end

"""
	maxfill(::Type{T}, args...)

Calculates the maximum physical number of strands that can fit for component `T`.
Custom shapes must overload this method.
"""
function maxfill end

"""
$(TYPEDSIGNATURES)

Applies [`PhysicalFillLimit`](@ref) to ensure the element count does not exceed the capacity computed by the extensible [`maxfill`](@ref) interface.

# Arguments

- `rule`: Rule instance with fields `n_field` and `geometry_fields` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws an `ArgumentError` on failure.
"""
@inline function Validation._apply(rule::PhysicalFillLimit, nt, ::Type{T}) where T
	n = getfield(nt, rule.n_field)

	# We only run the math if the types are sound. If they are garbage, we quietly 
	# return and let the dedicated type rules (like IntegerField or _ensure_real) throw.
	if n isa Integer
		geom_args = getfield.(Ref(nt), rule.geometry_fields)

		if all(x -> x isa Real, geom_args)
			limit = maxfill(T, geom_args...)

			n <= limit || throw(
				ArgumentError(
					"[$(_typename(T))] $(rule.n_field) $(_repr(n)) exceeds the physical maximum " *
					"limit ($limit) given geometry dimensions $(_repr(geom_args))",
				),
			)
		end
	end
	return nothing
end

"""
$(TYPEDSIGNATURES)

Applies [`Satisfies`](@ref) to evaluate an arbitrary predicate function against the specified fields.

# Arguments

- `rule`: Rule instance with `fields`, `predicate`, and `error_msg` \\[dimensionless\\].
- `nt`: Normalized `NamedTuple` of inputs \\[dimensionless\\].
- `::Type{T}`: Component type \\[dimensionless\\].

# Returns

- Nothing. Throws an `ArgumentError` on failure.

# Examples

```julia
Validation.extra_rules(::Type{SomeWeirdType}) = (
	# ... basic type rules ...
	
	Satisfies(
		(:width, :lay_angle, :overlap_pct),
		(w, a, pct) -> pct < 1.0 && pct >= 0.0 && w * cos(a) > 0,
		"Overlap percentage must be between 0 and 1, and effective width must be positive. Because I said so."
	)
)
```
"""
@inline function Validation._apply(rule::Satisfies, nt, ::Type{T}) where T
	args = getfield.(Ref(nt), rule.fields)

	# Evaluate the arbitrary predicate. We assume the predicate is robust enough
	# or that prior type-enforcing rules have already sanitized the inputs.
	if !rule.predicate(args...)
		throw(
			ArgumentError(
				"[$(_typename(T))] Validation failed for fields $(_repr(rule.fields)): " *
				"$(rule.error_msg) (Got values: $(_repr(args)))",
			),
		)
	end

	return nothing
end
