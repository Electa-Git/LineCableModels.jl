defaults(::Type{Formula{:wedepohl1973}}) = (;)

validate(formula::Formula{:wedepohl1973}, ::typeof(impedance!)) = formula
function validate(formula::Formula{:wedepohl1973}, t::Int, s::Int)
    t > 1 && s > 1 || invoke(validate, Tuple{Formula, Int, Int}, formula, t, s)
    return formula
end

"""
$(TYPEDSIGNATURES)

Warn once when the highest frequency of the scan leaves the bound of the paper
on the buried cables. The bound reads the static resistivity and permeability of the
deepest layer.

# Arguments

- `problem`: cable geometry, earth model and frequencies.
- `formula`: earth-return formula `Formula{:wedepohl1973}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:wedepohl1973})
    system = problem.system
    layer = last(problem.earth.layers)
    f = last(problem.frequencies)
    m = sqrt(2π * f * μ₀ * layer.mu_r / layer.rho)
    buried = findall(placement -> placement.y < 0, system.placements)
    reach = maximum(
        (i == j ? support(system.designs[i].outer) :
         hypot(system.placements[i].x - system.placements[j].x,
             system.placements[i].y - system.placements[j].y)
        for i in buried, j in buried); init = zero(m))
    m * reach < 0.25 ||
        @warn "the formula :wedepohl1973 holds for |m D| and |m r| below 0.25; on $(system.id) at $f Hz the largest is $(m * reach)"
    return problem
end

function (f::Formula{:wedepohl1973})(
        ::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    r_i = input.radius[i]
    h_i = abs(input.y[i])
    μ₂ = μ₀ * media[2].mu_r
    m = sqrt(s * μ₂ / media[2].rho)
    return s * μ₂ / (2π) *
           (-log(exp(MathConstants.γ) * m * r_i / 2) + 1 / 2 - 4m * h_i / 3)
end

function (f::Formula{:wedepohl1973})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    h_i = abs(input.y[i])
    h_j = abs(input.y[j])
    y = abs(input.x[i] - input.x[j])
    H = h_i + h_j
    D = hypot(y, h_i - h_j)
    μ₂ = μ₀ * media[2].mu_r
    m = sqrt(s * μ₂ / media[2].rho)
    return s * μ₂ / (2π) *
           (-log(exp(MathConstants.γ) * m * D / 2) + 1 / 2 - 2m * H / 3)
end
