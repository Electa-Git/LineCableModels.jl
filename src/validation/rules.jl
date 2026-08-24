"Base type for declarative validation rules."
abstract type Rule end

struct Finite <: Rule
    field::Symbol
end

struct Nonnegative <: Rule
    field::Symbol
end

struct Positive <: Rule
    field::Symbol
end

struct IntegerField <: Rule
    field::Symbol
end

struct Less <: Rule
    left::Symbol
    right::Symbol
end

struct LessEqual <: Rule
    left::Symbol
    right::Symbol
end

struct Greater <: Rule
    left::Symbol
    right::Symbol
end

struct GreaterEqual <: Rule
    left::Symbol
    right::Symbol
end

struct IsA{T} <: Rule
    field::Symbol
end

struct OneOf{T} <: Rule
    field::Symbol
    values::T
end

"""
$(TYPEDEF)

Run a named cross-field check declared beside its owning type.

$(TYPEDFIELDS)
"""
struct OwnerRule{F} <: Rule
    "Stable rule name used in diagnostics and tests."
    name::Symbol
    "Nonmutating callable that checks one complete value."
    check::F
end
