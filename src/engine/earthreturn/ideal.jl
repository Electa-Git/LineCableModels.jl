defaults(::Type{Formula{:ideal}}) = (;)

# The formula gives the potential coefficients alone.
validate(formula::Formula{:ideal}, ::typeof(potential!)) = formula

function (f::Formula{:ideal})(::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    ε₁ = vacuum_permittivity * media[1].eps_r
    input.y[i] > 0 || return zero(complex(ε₁))
    return Λ(f, Val(:self), zero(s), input.radius[i], input.y[i]) / (2π * ε₁)
end

function (f::Formula{:ideal})(::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    ε₁ = vacuum_permittivity * media[1].eps_r
    input.y[i] > 0 && input.y[j] > 0 || return zero(complex(ε₁))
    D = hypot(input.x[i] - input.x[j], input.y[i] - input.y[j])
    D′ = hypot(input.x[i] - input.x[j], input.y[i] + input.y[j])
    return Λ(f, Val(:mutual), zero(s), input.radius[i], D, D′) / (2π * ε₁)
end
