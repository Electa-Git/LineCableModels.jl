"""
    LineCableModels.Validation

Provide nonmutating, declarative validation rules for project-owned values.
Concrete types specialize [`rules`](@ref); [`validate`](@ref) applies those
rules in declaration order and returns its input unchanged.
"""
module Validation

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: maxfill

export Rule, rules, validate
export Finite, Nonnegative, Positive, IntegerField, Less, LessEqual, Greater,
       GreaterEqual, IsA, OneOf, PhysicalFillLimit, OwnerRule

"Return the declarative validation rules owned by a type."
rules(::Type) = ()

include("rules.jl")
include("applyrules.jl")

"""
$(TYPEDSIGNATURES)

Apply the ordered rules owned by the concrete type of `value` and return
`value` unchanged.
"""
validate(value) = check(typeof(value), value)

end
