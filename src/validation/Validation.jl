"""
    LineCableModels.Validation

Provide nonmutating, declarative validation rules for project-owned values.
Concrete types specialize [`rules`](@ref) and may add a root [`validate`](@ref)
method for checks that are clearer as direct dispatch.
"""
module Validation

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: maxfill, validate

export Rule, rules, check
export Finite, Nonnegative, Positive, IntegerField, Less, LessEqual, Greater,
       GreaterEqual, IsA, OneOf, PhysicalFillLimit

"Return the declarative validation rules owned by a type."
rules(::Type) = ()

include("rules.jl")
include("applyrules.jl")

end
