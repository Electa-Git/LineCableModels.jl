"""
    LineCableModels.Validation

Define nonmutating validation rules for LineCableModels values.
Concrete types specialise [`rules`](@ref). [`validate`](@ref) applies those
rules in declaration order and returns its input unchanged.
"""
module Validation

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
export Rule, rules, validate
export Finite, Nonnegative, Positive, IntegerField, Less, LessEqual, Greater,
       GreaterEqual, IsA, OneOf, OwnerRule

include("interfaces.jl")
include("rules.jl")
include("validate.jl")

end
