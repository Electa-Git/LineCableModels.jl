"""
$(TYPEDEF)

Store cable constants per unit length.

The fields `R`, `L`, and `C` are stored in Ω/m, H/m, and F/m respectively.
Display conversions belong to [`LineCableModels.Units`](@ref) and presentation
adapters.

$(TYPEDFIELDS)
"""
struct CableConstants{T <: Number} <: AbstractCoreResult
    "Series resistance per unit length \\[Ω/m\\]."
    R::T
    "Series inductance per unit length \\[H/m\\]."
    L::T
    "Shunt capacitance per unit length \\[F/m\\]."
    C::T
end

function Base.:(==)(left::CableConstants, right::CableConstants)
    left.R == right.R && left.L == right.L && left.C == right.C
end

function CableConstants(R::Real, L::Real, C::Real)
    values = promote(R, L, C)
    return CableConstants{typeof(first(values))}(values...)
end

observe(constants::CableConstants, ::typeof(R)) = constants.R
observe(constants::CableConstants, ::typeof(L)) = constants.L
observe(constants::CableConstants, ::typeof(C)) = constants.C

R(constants::CableConstants) = observe(constants, R)
L(constants::CableConstants) = observe(constants, L)
C(constants::CableConstants) = observe(constants, C)
basis(::CableConstants) = :pul
resistance(constants::CableConstants) = R(constants)
inductance(constants::CableConstants) = L(constants)
capacitance(constants::CableConstants) = C(constants)

observables(::Type{<:CableConstants}) = (R, L, C)
