defaults(::Type{Formula{:ideal}}) = (;)
domain(::Formula{:ideal}) = (Γ = 0.0 .. 0.0,)

# The formula gives the potential coefficients alone.
validate(formula::Formula{:ideal}, ::typeof(potential!)) = formula
validate(formula::Formula{:ideal}, ::Int, ::Int) = formula

function (f::Formula{:ideal})(::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    ε₁ = ε₀ * media[1].eps_r
    input.earth_layer[i] == 1 || return zero(complex(ε₁))
    return Λ(f, Val(:self), zero(s), input.radius[i], input.y[i]) / (2π * ε₁)
end

function (f::Formula{:ideal})(::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    ε₁ = ε₀ * media[1].eps_r
    input.earth_layer[i] == 1 && input.earth_layer[j] == 1 || return zero(complex(ε₁))
    D = hypot(input.x[i] - input.x[j], input.y[i] - input.y[j])
    D′ = hypot(input.x[i] - input.x[j], input.y[i] + input.y[j])
    return Λ(f, Val(:mutual), zero(s), input.radius[i], D, D′) / (2π * ε₁)
end
