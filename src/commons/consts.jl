"""
$(TYPEDSIGNATURES)

Return the vacuum permittivity ``\\varepsilon_0 = 8.8541878128 \\times 10^{-12}``
\\[F/m\\] in the scalar type `T`.

The integer mantissa and the power of ten are evaluated in `T`. The value is exact
for rational `T` and has the precision of `T` otherwise.
"""
vacuum_permittivity(::Type{T}) where {T <: Real} = one(T) * 88541878128 * (one(T) * 10)^(-22)
