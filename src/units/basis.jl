"""
$(TYPEDEF)

The length basis of line quantities: [`PerUnitLength`](@ref), written `:pul`, or
[`TotalLength`](@ref), written `:total`. `LengthBasis(name)` gives the basis that a
Symbol writes and rejects any other Symbol.
"""
abstract type LengthBasis end

"""
$(TYPEDEF)

Values per unit length, such as Ω/m. Written `:pul`.
"""
struct PerUnitLength <: LengthBasis end

"""
$(TYPEDEF)

Values for the total line length, such as Ω. Written `:total`.
"""
struct TotalLength <: LengthBasis end

LengthBasis(basis::LengthBasis) = basis
LengthBasis(name::Symbol) = LengthBasis(Val(name))
LengthBasis(::Val{:pul}) = PerUnitLength()
LengthBasis(::Val{:total}) = TotalLength()
LengthBasis(::Val{name}) where {name} =
    throw(ArgumentError("basis must be :pul or :total; got $(repr(name))"))
Base.Symbol(::PerUnitLength) = :pul
Base.Symbol(::TotalLength) = :total
