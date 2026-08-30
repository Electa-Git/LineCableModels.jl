routes(::Val{:Carson}) = (
    self = carson_self,
    mutual = carson_mutual,
    Γ = carson_gamma
)

function assumptions(::Val{:Carson})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Carson}) = Val(:zero)
description(::Formula{:Carson}) = "Carson"

function carson_gamma(jω, permeability, permittivity)
    value = zero(jω)
    return (Γ = value, squared = value)
end

function (formula::Formula{:Carson})(
        resistivity::AbstractVector{T},
        permittivity::AbstractVector{T},
        permeability::AbstractVector{T},
        jω::Complex{T},
        Γ,
        segments = nothing
) where {T <: Real}
    _check(resistivity, permittivity, permeability)
    values = formula.assumptions
    source = 1
    other = 2
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
    return Functor{:Carson, typeof(formula.routes), typeof(state)}(
        formula.routes,
        state
    )
end

carson_self(functor, pair) = _impedance(functor, pair)
carson_mutual(functor, pair) = _impedance(functor, pair)

@inline function (functor::Functor{:Carson})(::Val{:self}, pair)
    return functor.routes.self(functor, pair)
end

@inline function (functor::Functor{:Carson})(::Val{:mutual}, pair)
    return functor.routes.mutual(functor, pair)
end

:Carson
