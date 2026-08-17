"""
$(TYPEDEF)

Select frequency-invariant earth resistivity, permittivity, and permeability.
"""
struct CPEarth <: AbstractEarthPropertiesFormulation end

description(::CPEarth) = "Constant properties (CP)"

"""
$(TYPEDSIGNATURES)

Evaluate static earth layers at every requested frequency. The returned
matrices have one row per layer and one column per frequency. Permittivity and
permeability are absolute SI values.
"""
function evaluate(
        ::CPEarth,
        model::EarthModel{T},
        frequencies::AbstractVector{T}
) where {T <: Real}
    nlayer = length(model.layers)
    nfrequency = length(frequencies)
    rho = Matrix{T}(undef, nlayer, nfrequency)
    epsilon = Matrix{T}(undef, nlayer, nfrequency)
    mu = Matrix{T}(undef, nlayer, nfrequency)
    unit = one(first(frequencies))
    epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    @inbounds for row in eachindex(model.layers)
        layer = model.layers[row]
        for column in eachindex(frequencies)
            rho[row, column] = layer.rho
            epsilon[row, column] = epsilon0 * layer.eps_r
            mu[row, column] = mu0 * layer.mu_r
        end
    end
    return (; rho, epsilon, mu)
end
