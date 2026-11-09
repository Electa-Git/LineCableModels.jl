defaults(::Type{Formula{:linear}}) = (;)
domain(::Formula{:linear}) = (Δθ = 0.0 .. 150.0,)

function (::Formula{:linear})(material::Material{T}, θ::T) where {T <: Real}
    difference = θ - material.T0
    factor = 1 + material.alpha * difference
    isfinite(factor) && factor > 0 || throw(DomainError(factor,
        "the factor 1 + α (θ − T0) of the linear law must be positive and finite"))
    return isinf(material.rho) ? material.rho : material.rho * factor
end
