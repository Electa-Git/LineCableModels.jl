@inline function _full(jω, permeability, conductivity, permittivity)
    return sqrt(jω * permeability * (conductivity + jω * permittivity))
end

@inline function _vacuum(jω, permeability, conductivity, permittivity)
    return jω * sqrt(permeability * vacuum_permittivity(permittivity))
end

@inline _material(permeability) = permeability

function _check(resistivity, permittivity, permeability)
    length(resistivity) == length(permittivity) == length(permeability) ||
        throw(DimensionMismatch("earth-property vectors must have equal lengths"))
    length(resistivity) >= 2 || throw(DimensionMismatch(
        "an earth-admittance formula requires air and at least one earth layer"
    ))
    return nothing
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

"""
$(TYPEDEF)

Store pair-specific geometry for one improper earth-admittance integral.

$(TYPEDFIELDS)
"""
struct Integrand{F, T}
    "Formula-owned frequency functor."
    functor::F
    "Sum of conductor depths or heights \\[m\\]."
    height_sum::T
    "Horizontal conductor separation \\[m\\]."
    separation::T
end

@inline function (integrand::Integrand)(lambda::Real)
    state = integrand.functor.state
    source_attenuation = sqrt(
        lambda * lambda + state.gamma_source_squared + state.gamma_squared
    )
    other_attenuation = sqrt(
        lambda * lambda + state.gamma_other_squared + state.gamma_squared
    )
    common = source_attenuation * state.other_permeability +
             other_attenuation * state.source_permeability
    decay = exp(-source_attenuation * integrand.height_sum)
    magnetic = state.other_permeability * decay / common
    coupling_numerator = state.other_permeability *
                         state.source_permeability *
                         source_attenuation *
                         (state.gamma_source_squared -
                          state.gamma_other_squared) * decay
    coupling_denominator = common * (
        source_attenuation * state.gamma_other_squared *
        state.source_permeability +
        other_attenuation * state.gamma_source_squared *
        state.other_permeability
    )
    coupling = coupling_numerator / coupling_denominator
    return (magnetic + coupling) * cos(integrand.separation * lambda)
end

function _integral(functor::Functor, height_sum, separation)
    integrand = Integrand(functor, height_sum, separation)
    state = functor.state
    R = typeof(state.tolerance)
    value, _ = quadgk(
        integrand,
        zero(R),
        R(Inf);
        rtol = state.tolerance,
        segbuf = state.segments,
        norm = z -> abs(complex(nominal(real(z)), nominal(imag(z))))
    )
    return 2value
end

function _pair(functor::Functor, pair)
    state = functor.state
    pair.layers[1] == state.source_layer || throw(ArgumentError(
        "source conductor is in layer $(pair.layers[1]) but formula :$(formula_id(state.formula)) expects layer $(state.source_layer)"
    ))
    pair.layers[2] == state.target_layer || throw(ArgumentError(
        "target conductor is in layer $(pair.layers[2]) but formula :$(formula_id(state.formula)) expects layer $(state.target_layer)"
    ))
    return nothing
end

@inline function _admittance(functor::Functor, pair)
    _pair(functor, pair)
    state = functor.state
    T = typeof(pair.separation)
    height_i = abs(pair.heights[1])
    height_j = abs(pair.heights[2])
    direct_distance = hypot(pair.separation, height_i - height_j)
    image_distance = hypot(pair.separation, height_i + height_j)
    perfect_ground = bessel_difference(
        state.gamma_source,
        direct_distance,
        image_distance
    )
    complex_conductivity = state.source_conductivity +
                           state.jω * state.source_permittivity
    if state.source_layer == 1 && iszero(state.gamma_squared) &&
       isapprox(nominal(real(state.gamma_other)), 0; atol = 1.0e-6)
        return state.jω /
               (2 * (one(T) * π) * complex_conductivity) * perfect_ground
    end
    if state.source_layer == 2 && state.target_layer == 2 &&
       iszero(state.gamma_squared) &&
       isapprox(nominal(real(state.gamma_source)), 0; atol = 1.0e-6)
        return state.jω /
               (2 * (one(T) * π) * complex_conductivity) * perfect_ground
    end
    correction = _integral(
        functor,
        height_i + height_j,
        pair.separation
    )
    return state.jω /
           (2 * (one(T) * π) * complex_conductivity) *
           (perfect_ground + correction)
end
