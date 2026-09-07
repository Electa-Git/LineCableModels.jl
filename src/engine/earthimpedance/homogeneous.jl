@inline function _full(jω, permeability, conductivity, permittivity)
    return sqrt(jω * permeability * (conductivity + jω * permittivity))
end

@inline function _lossless(jω, permeability, conductivity, permittivity)
    return jω * sqrt(permeability * permittivity)
end

@inline function _conductive(jω, permeability, conductivity, permittivity)
    return sqrt(jω * permeability * conductivity)
end

@inline function _require_horizontal_separation(pair)
    iszero(pair.separation) && throw(DomainError(
        pair.separation,
        "this mutual closed form requires nonzero horizontal cable separation"
    ))
    return nothing
end

@inline _material(permeability) = permeability

function _image_plane_impedance(functor, pair, depth, interaction::Val)
    _require(pair, Val(:overhead))
    state, geometry = functor.state, _geometry(pair)
    πT = one(geometry.y_ij) * π
    if interaction === Val(:self)
        ratio = (geometry.H + 2depth) / geometry.y_ij
    else
        image = sqrt((geometry.H + 2depth)^2 + geometry.y_ij^2)
        ratio = image / geometry.d_ij
    end
    return state.jω * state.mu[1] / (2πT) * log(ratio)
end

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

function _require(pair, expected)
    typeof(_placement(pair)) === typeof(expected) || throw(ArgumentError(
        "earth-impedance formula is incompatible with this conductor placement"
    ))
    return nothing
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
    if hasproperty(state.formula.assumptions,:quadrature) &&
            state.formula.assumptions.quadrature===:double_exponential
        points=[sqrt(-real(g)) for g in state.gamma_medium_squared if
            iszero(imag(g)) && real(g)<0]
        return _complex_result(state.jω,double_exponential(guarded,R;
            rtol=state.tolerance,breakpoints=points).value)
    end
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

# Resolve the near-zero spectral interval for large physical heights. Without
# this substitution, all initial quadrature nodes can underflow to zero.
function _height_quadrature(integrand::F, state, height) where {F <: Function}
    scale=max(height,one(height))
    return _quadrature(state,t->integrand(t/scale)/scale)
end

# The selected material state determines whether soil displacement current
# is retained. The overhead self correction has zero lateral separation;
# the finite conductor radius belongs only to the ideal-ground logarithm.
function _carson_integral(state,height,lateral)
    return _height_quadrature(state,height) do λ
        exp(-height*λ)*cos(lateral*λ)/
            (λ+sqrt(λ^2+state.gamma_medium_squared[2]))
    end
end

function _homogeneous_overhead_coefficient(state,pair)
    geometry=_geometry(pair)
    self=pair.row==pair.column
    lateral=self ? zero(geometry.H) : geometry.y_ij
    ideal=self ? log(geometry.H/geometry.y_ij) : log(geometry.D_ij/geometry.d_ij)
    integral=_carson_integral(state,geometry.H,lateral)
    return state.jω*state.mu[1]/(2*(one(geometry.H)*π))*(ideal+2integral)
end

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
    _check(resistivity, permittivity, permeability)
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
    _check(resistivity, permittivity, permeability)
    length(thickness) == length(resistivity) || throw(DimensionMismatch(
        "earth-layer thickness and material vectors must align"
    ))
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

function _check(resistivity, permittivity, permeability)
    length(resistivity) == length(permittivity) == length(permeability) ||
        throw(DimensionMismatch("earth-property vectors must have equal lengths"))
    length(resistivity) >= 2 || throw(DimensionMismatch(
        "an earth-impedance formula requires air and at least one earth layer"
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

Store pair-specific geometry for one improper earth-return integral.

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
    decay = exp(-source_attenuation * integrand.height_sum)
    denominator = source_attenuation * state.other_permeability +
                  other_attenuation * state.source_permeability
    return state.other_permeability * decay / denominator *
           cos(integrand.separation * lambda)
end

function _integral(functor::Functor, height_sum, separation)
    integrand = Integrand(functor, height_sum, separation)
    state = functor.state
    state.segments === nothing || empty!(state.segments)
    R = typeof(state.tolerance)
    value, _ = quadgk(
        integrand,
        zero(R),
        R(Inf);
        rtol = state.tolerance,
        segbuf = state.segments,
        norm = z -> abs(complex(nominal(real(z)), nominal(imag(z))))
    )
    return _complex_result(state.jω, 2value)
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

@inline function _impedance(functor::Functor, pair)
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
    correction = _integral(
        functor,
        height_i + height_j,
        pair.separation
    )
    return state.jω * state.source_permeability /
           (2 * (one(T) * π)) * (perfect_ground + correction)
end

# K₂(z)-2/z², with its divergent leading term removed before evaluation.
function _sunde_bessel_remainder(z)
    T=typeof(real(z))
    logarithm=log(z/2)+T(Base.MathConstants.eulergamma)
    result=-one(z)/2; power=one(z); harmonic=zero(T)
    for n in 1:10000
        power*=z^2/(4n^2); harmonic+=inv(T(n))
        term=power*(T(n)/(n+1)*(harmonic-logarithm)-inv(T(2)*(n+1)^2))
        result+=term
        abs(term)<=eps(T)*max(one(T),abs(result)) && return result
    end
    throw(ErrorException("Sunde Bessel remainder did not converge"))
end

# [1-(1+hz)exp(-hz)]/z², including its finite z=0 limit.
function _sunde_exponential_remainder(z,h)
    T=typeof(real(z)); term=one(z)*h^2/2; result=term
    for n in 2:10000
        term*=(-h*z)*n/((n+1)*(n-1))
        result+=term
        abs(term)<=eps(T)*max(one(T),abs(result)) && return result
    end
    throw(ErrorException("Sunde exponential remainder did not converge"))
end

function _sunde_image_correction(γ,H,r)
    D=hypot(H,r); z=γ*D; h=H/D
    angular=(H^2-r^2)/D^2
    combined=if abs(z)<one(H)/2
        _sunde_bessel_remainder(z)+2 * _sunde_exponential_remainder(z,h)
    else
        special_besselk(2,z)-2exp(-γ*H)*(1+γ*H)/z^2
    end
    return angular*combined
end

function _sunde_image_terms(γ,d,H,r)
    return special_besselk(0,γ*d)+_sunde_image_correction(γ,H,r)
end

# Surface field ratio q/mu, matched upward through each physical layer.
# Callers prescribe a longitudinal reference for which the air root is lambda.
function _layered_input_step(q,weight,d,load)
    if iszero(q)
        transmission=inv(1+load*d/weight)
        return (input=load*transmission,transmission)
    end
    local_ratio=q*weight; t=tanh(q*d); decay=exp(-q*d)
    input=local_ratio*(load+local_ratio*t)/(local_ratio+load*t)
    transmission=2local_ratio*decay/
        (local_ratio+load+(local_ratio-load)*decay^2)
    return (;input,transmission)
end

function _layered_overhead_kernel(lambda,state)
    last=length(state.rho)
    input=sqrt(lambda^2+state.gamma_medium_squared[last]+state.gamma_squared)/state.mu[last]
    for layer in (last-1):-1:2
        q=sqrt(lambda^2+state.gamma_medium_squared[layer]+state.gamma_squared)
        input=_layered_input_step(q,inv(state.mu[layer]),state.thickness[layer],input).input
    end
    return inv(lambda+state.mu[1]*input)
end

function _layered_overhead_coefficient(state,pair)
    geometry=_geometry(pair)
    self=pair.row==pair.column
    lateral=self ? zero(geometry.H) : geometry.y_ij
    ideal=self ? log(geometry.H/geometry.y_ij) : log(geometry.D_ij/geometry.d_ij)
    integral=_height_quadrature(state,geometry.H) do lambda
        _layered_overhead_kernel(lambda,state)*
            exp(-lambda*geometry.H)*cos(lambda*lateral)
    end
    return state.jω*state.mu[1]/(2*(one(geometry.H)*π))*(ideal+2integral)
end
