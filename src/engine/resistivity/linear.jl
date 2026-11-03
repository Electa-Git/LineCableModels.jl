defaults(::Type{Formula{:linear}}) = (;)

function (::Formula{:linear})(material::Material{T}, θ::T) where {T <: Real}
    difference = θ - material.T0
    abs(difference) < 150 || throw(DomainError(θ,
        "the linear resistivity law holds within 150 K of Material.T0 = $(material.T0) °C"))
    factor = 1 + material.alpha * difference
    isfinite(factor) && factor > 0 || throw(DomainError(factor,
        "the factor 1 + α (θ − T0) of the linear law must be positive and finite"))
    return isinf(material.rho) ? material.rho : material.rho * factor
end
