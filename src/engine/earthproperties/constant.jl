defaults(::Type{Formula{:constant}}) = (;)

"""
$(TYPEDSIGNATURES)

Return the earth `layer` at the frequency `f` by the formula `Formula{:constant}`, whose equation
is in the section "Formulas" of [`EarthProperties`](@ref).

# Arguments

- `Formula{:constant}`: frequency-independent properties law.
- `layer`: layer of the earth or of the air.
- `f`: frequency \\[Hz\\].

# Returns

- `layer`, whose properties hold at every frequency.
"""
function (::Formula{:constant})(layer::EarthLayer{T}, f::T) where {T <: Real}
    return layer
end
