defaults(::Type{Formula{:ideal}}) = (;)

# The formula gives the potential coefficients alone.
validate(formula::Formula{:ideal}, ::typeof(potential!)) = formula

"""
$(TYPEDSIGNATURES)

Return an entry of the earth-return potential coefficients \\[m/F\\] of `Formula{:ideal}`. Its
equations are in [`EarthReturn`](@ref).
"""
function (f::Formula{:ideal})(::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    ε₁ = vacuum_permittivity * media[1].eps_r
    input.y[i] > 0 || return zero(ε₁)
    return Λ(f, Val(:self), zero(s), input.radius[i], input.y[i]) / (2π * ε₁)
end

"""
$(TYPEDSIGNATURES)

Return an entry of the earth-return potential coefficients \\[m/F\\] of `Formula{:ideal}`. Its
equations are in [`EarthReturn`](@ref).
"""
function (f::Formula{:ideal})(::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    ε₁ = vacuum_permittivity * media[1].eps_r
    input.y[i] > 0 && input.y[j] > 0 || return zero(ε₁)
    D = hypot(input.x[i] - input.x[j], input.y[i] - input.y[j])
    D′ = hypot(input.x[i] - input.x[j], input.y[i] + input.y[j])
    return Λ(f, Val(:mutual), zero(s), input.radius[i], D, D′) / (2π * ε₁)
end
