function routes(::Val{:Pollaczek})
    (
        self = pollaczek_self,
        mutual = pollaczek_mutual,
        Γ = pollaczek_gamma
    )
end

function assumptions(::Val{:Pollaczek})
    (
        air = _vacuum,
        earth = _vacuum,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pollaczek}) = Val(:zero)
description(::Formula{:Pollaczek}) = "Pollaczek"

function pollaczek_gamma(jω, permeability, permittivity)
    value = zero(jω)
    return (Γ = value, squared = value)
end

function (formula::Formula{:Pollaczek})(
        resistivity::AbstractVector{T},
        permittivity::AbstractVector{T},
        permeability::AbstractVector{T},
        jω::Complex{T},
        Γ,
        segments = nothing
) where {T <: Real}
    _check(resistivity, permittivity, permeability)
    values = formula.assumptions
    source = 2
    other = 1
    source_conductivity = conductivity(resistivity[source])
    other_conductivity = conductivity(resistivity[other])
    mu_source = _permeability(permeability, source, values.permeability)
    mu_other = _permeability(permeability, other, values.permeability)
    gamma_source = _wave(
        values,
        source,
        jω,
        mu_source,
        source_conductivity,
        permittivity
    )
    gamma_other = _wave(
        values,
        other,
        jω,
        mu_other,
        other_conductivity,
        permittivity
    )
    gamma_source_squared = gamma_source^2
    gamma_other_squared = gamma_other^2
    longitudinal = _longitudinal(
        formula,
        Γ,
        jω,
        mu_source,
        permittivity[source]
    )
    R = typeof(float(nominal(one(T))))
    tolerance = max(R(1e-8), eps(R))
    state = (;
        formula,
        jω,
        Γ = longitudinal.Γ,
        gamma_source,
        gamma_other,
        gamma_source_squared,
        gamma_other_squared,
        gamma_squared = longitudinal.squared,
        source_permeability = mu_source,
        other_permeability = mu_other,
        source_conductivity,
        source_permittivity = permittivity[source],
        source_layer = source,
        target_layer = source,
        tolerance,
        segments
    )
    return Functor{:Pollaczek, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

pollaczek_self(functor, pair) = _admittance(functor, pair)
pollaczek_mutual(functor, pair) = _admittance(functor, pair)

@inline function (functor::Functor{:Pollaczek})(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor{:Pollaczek})(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

:Pollaczek
