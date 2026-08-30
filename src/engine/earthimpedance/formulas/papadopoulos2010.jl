function routes(::Val{:Papadopoulos2010})
    (
        self = papadopoulos_self,
        mutual = papadopoulos_mutual,
        Γ = papadopoulos_gamma
    )
end

function assumptions(::Val{:Papadopoulos2010})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2010}) = Val(:explicit)
description(::Formula{:Papadopoulos2010}) = "Papadopoulos"

function papadopoulos_gamma(jω, permeability, permittivity)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2010})(
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
    mu_source = _permeability(permeability, source, values.permeability)
    mu_other = _permeability(permeability, other, values.permeability)
    gamma_source = _wave(
        values,
        source,
        jω,
        mu_source,
        conductivity(resistivity[source]),
        permittivity
    )
    gamma_other = _wave(
        values,
        other,
        jω,
        mu_other,
        conductivity(resistivity[other]),
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
        gamma_source_squared,
        gamma_other_squared,
        gamma_squared = longitudinal.squared,
        source_permeability = mu_source,
        other_permeability = mu_other,
        source_layer = source,
        target_layer = source,
        tolerance,
        segments
    )
    return Functor{:Papadopoulos2010, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

papadopoulos_self(functor, pair) = _impedance(functor, pair)
papadopoulos_mutual(functor, pair) = _impedance(functor, pair)

@inline function (functor::Functor{:Papadopoulos2010})(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor{:Papadopoulos2010})(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

:Papadopoulos2010
