@inline function _full(jω, permeability, conductivity, permittivity)
    return sqrt(jω * permeability * (conductivity + jω * permittivity))
end

@inline function _vacuum(jω, permeability, conductivity, permittivity)
    return jω * sqrt(permeability * vacuum_permittivity(permittivity))
end

@inline _material(permeability) = permeability

@inline function _complex_result(::Complex{T}, value)::Complex{T} where {T <: Real}
    return value
end

@inline function _placement(pair)
    source_air = pair.layers[1] == 1
    target_air = pair.layers[2] == 1
    source_air && target_air && return Val(:overhead)
    !source_air && !target_air && return Val(:underground)
    return Val(:mixed)
end

@inline function _geometry(pair)
    h_i = abs(pair.heights[1])
    h_j = abs(pair.heights[2])
    y_ij = pair.separation
    H = h_i + h_j
    d_ij = hypot(y_ij, h_i - h_j)
    D_ij = hypot(y_ij, H)
    return (; h_i, h_j, y_ij, H, d_ij, D_ij)
end

function _quadrature(state, integrand)
    R = typeof(state.tolerance)
    state.segments === nothing || empty!(state.segments)
    radial_limit = sqrt(floatmax(R)) / 2
    guarded = lambda -> lambda > radial_limit ? zero(state.jω) : integrand(lambda)
    value, _ = quadgk(
        guarded,
        zero(R),
        R(Inf);
        rtol = state.tolerance,
        segbuf = state.segments,
        norm = z -> abs(complex(nominal(real(z)), nominal(imag(z))))
    )
    return _complex_result(state.jω, value)
end

@inline _quadrature(integrand::F, state::NamedTuple) where {F <: Function} = _quadrature(state, integrand)

@inline function _tolerance(::Type{T}) where {T}
    R = typeof(float(nominal(one(T))))
    return max(R(1e-8), eps(R))
end

function _homogeneous_functor(
        ::Val{ID},
        formula::Formula{ID},
        resistivity::AbstractVector{T},
        permittivity::AbstractVector{T},
        permeability::AbstractVector{T},
        jω::Complex{T},
        Γ,
        segments
) where {ID, T <: Real}
    validate(formula, resistivity, permittivity, permeability, nothing)
    values = formula.assumptions
    μ = (
        _permeability(permeability, 1, values.permeability),
        _permeability(permeability, 2, values.permeability)
    )
    σ = (conductivity(resistivity[1]), conductivity(resistivity[2]))
    γ = (
        _wave(values, 1, jω, μ[1], σ[1], permittivity),
        _wave(values, 2, jω, μ[2], σ[2], permittivity)
    )
    longitudinal = _longitudinal(formula, Γ, jω, μ[2], permittivity[2])
    state = (;
        formula,
        jω,
        Γ = longitudinal.Γ,
        gamma_squared = longitudinal.squared,
        rho = (resistivity[1], resistivity[2]),
        epsilon = (permittivity[1], permittivity[2]),
        mu = μ,
        sigma = σ,
        gamma = γ,
        gamma_medium_squared = (γ[1]^2, γ[2]^2),
        tolerance = _tolerance(T),
        segments
    )
    return Functor{ID, typeof(formula.routes), typeof(state)}(formula.routes, state)
end

function _stratified_functor(
        ::Val{ID},
        formula::Formula{ID},
        resistivity::AbstractVector{T},
        permittivity::AbstractVector{T},
        permeability::AbstractVector{T},
        jω::Complex{T},
        Γ,
        segments,
        thickness::AbstractVector{T}
) where {ID, T <: Real}
    validate(formula, resistivity, permittivity, permeability, thickness)
    values = formula.assumptions
    μ = map(eachindex(permeability)) do layer
        _permeability(permeability, layer, values.permeability)
    end
    σ = conductivity.(resistivity)
    γ = map(eachindex(resistivity)) do layer
        _wave(values, layer, jω, μ[layer], σ[layer], permittivity)
    end
    longitudinal = _longitudinal(formula, Γ, jω, μ[2], permittivity[2])
    state = (;
        formula,
        jω,
        Γ = longitudinal.Γ,
        gamma_squared = longitudinal.squared,
        rho = resistivity,
        epsilon = permittivity,
        mu = μ,
        sigma = σ,
        gamma = γ,
        gamma_medium_squared = γ .^ 2,
        thickness,
        tolerance = _tolerance(T),
        segments
    )
    return Functor{ID, typeof(formula.routes), typeof(state)}(formula.routes, state)
end

@inline function _permeability(values, layer, transform)
    value = values[layer]
    return layer == 1 ? value : transform(value)
end

@inline function _wave(
        values,
        layer,
        jω,
        permeability,
        conductivity_value,
        permittivity
)
    evaluator = layer == 1 ? values.air : values.earth
    return evaluator(jω, permeability, conductivity_value, permittivity[layer])
end
