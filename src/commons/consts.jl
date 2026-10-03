"""
$(TYPEDSIGNATURES)

Return the vacuum permittivity ``\\varepsilon_0 = 8.8541878128 \\times 10^{-12}``
\\[F/m\\] in the scalar type `T`.

The integer mantissa and the power of ten are evaluated in `T`. The value is exact
for rational `T` and has the precision of `T` otherwise.
"""
vacuum_permittivity(::Type{T}) where {T <: Real} = one(T) * 88541878128 * (one(T) * 10)^(-22)

"""
$(TYPEDSIGNATURES)

Return the vacuum permeability ``\\mu_0 = 4\\pi \\times 10^{-7}`` \\[H/m\\] in the scalar
type `T`.

The factors are evaluated in `T`. The value has the precision of `T`.
"""
vacuum_permeability(::Type{T}) where {T <: Real} = one(T) * 4 * (one(T) * π) * (one(T) * 10)^(-7)
