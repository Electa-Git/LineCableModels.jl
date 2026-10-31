defaults(::Type{Formula{:wedepohl1973}}) = (;)

function (::Formula{:wedepohl1973})(
        r_in::T, r_ex::T, ρ::T, μ_r::T, s::Complex{T}) where {T <: Real}
    μ = vacuum_permeability * μ_r
    σ = inv(ρ)
    m = sqrt(s * μ * σ)
    if iszero(r_in)
        outer = (ρ * m / (2π * r_ex)) * coth(0.777 * m * r_ex) +
                0.356 * ρ / (π * r_ex^2)
        return (inner = zero(outer), outer, transfer = zero(outer))
    end
    Δ = r_ex - r_in
    inner = (ρ * m / (2π * r_in)) * coth(m * Δ) -
            ρ / (2π * r_in * (r_in + r_ex))
    outer = (ρ * m / (2π * r_ex)) * coth(m * Δ) +
            ρ / (2π * r_ex * (r_in + r_ex))
    transfer = (ρ * m / (π * (r_in + r_ex))) * csch(m * Δ)
    return (; inner, outer, transfer)
end

"""
$(TYPEDSIGNATURES)

Warn once for each conductor whose tube leaves the bound of the paper.

# Arguments

- `problem`: cable geometry and earth model.
- `formula`: internal-impedance formula `Formula{:wedepohl1973}` that consumes the problem.

# Returns

- `problem`, unchanged.
"""
function validate(problem::LineParametersProblem, formula::Formula{:wedepohl1973})
    for (k, design) in enumerate(problem.system.designs)
        for tube in flatten(design).conductors
            ratio = (tube.r_ex - tube.r_in) / (tube.r_ex + tube.r_in)
            iszero(tube.r_in) || ratio < 1 / 8 ||
                @warn "the formula :wedepohl1973 holds for (r_ex - r_in) / (r_ex + r_in) < 1/8; conductor $(tube.terminal) of the cable $k of $(problem.system.id) has $ratio"
        end
    end
    return problem
end
