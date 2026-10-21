"""
$(TYPEDSIGNATURES)

Return the vacuum permeability ``\\mu_0 = 4\\pi \\times 10^{-7}`` \\[H/m\\] in the scalar
type `T`.

The factors are evaluated in `T`. The value has the precision of `T`.
"""
vacuum_permeability(::Type{T}) where {T <: Real} = one(T) * 4 * (one(T) * π) * (one(T) * 10)^(-7)
