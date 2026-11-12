defaults(::Type{Formula{:linear}}) = (;)
domain(::Formula{:linear}) = (Δθ = 0.0 .. 150.0,)

# The factor 1 + α (θ − T0) of each conductor at the operating temperature is positive and
# finite.
function validate(problem::LineParametersProblem, formula::Formula{:linear})
    factors = [1 +
               region.source.material.alpha *
               (problem.temperature - region.source.material.T0)
               for region in problem.system.regions
               if region.source.material isa Material{<:Real, :conductor}]
    all(factor -> isfinite(factor) && factor > 0, factors) || throw(DomainError(factors,
        "the factor 1 + α (θ − T0) of the linear law must be positive and finite"))
    return problem
end

function (::Formula{:linear})(material::Material{T}, θ::T) where {T <: Real}
    difference = θ - material.T0
    factor = 1 + material.alpha * difference
    return isinf(material.rho) ? material.rho : material.rho * factor
end
