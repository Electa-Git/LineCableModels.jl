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

Rule that enforces a strict ordering constraint `a > b` between two fields.

$(TYPEDFIELDS)
"""
struct Greater <: Rule
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

Rule that enforces a non‑strict ordering constraint `a ≥ b` between two fields.

$(TYPEDFIELDS)
"""
struct GreaterEq <: Rule
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
struct OneOf{S} <: Rule
	name::Symbol
	set::S
end

"""
$(TYPEDEF)

Rule that enforces that the number of discrete elements (e.g., wires or strands) does not exceed the theoretical physical maximum for a given geometry.

$(TYPEDFIELDS)
"""
struct PhysicalFillLimit <: Rule
	"Symbol representing the field containing the element count (e.g., `:num_wires`)."
	n_field::Symbol
	"Tuple of symbols representing the geometric fields required to compute the limit."
	geometry_fields::Tuple{Vararg{Symbol}}
end

"""
$(TYPEDEF)

A highly flexible rule that enforces an arbitrary user-defined predicate across one or more fields.
Useful for complex, cross-field physics or one-off geometrical constraints without polluting the ruleset with bespoke types.

$(TYPEDFIELDS)
"""
struct Satisfies <: Rule
	"Tuple of symbols representing the fields to be evaluated."
	fields::Tuple{Vararg{Symbol}}
	"A function (often anonymous) that accepts the extracted fields as arguments and returns a boolean."
	predicate::Function
	"The diagnostic message appended to the error if the predicate returns `false`."
	error_msg::String
end
