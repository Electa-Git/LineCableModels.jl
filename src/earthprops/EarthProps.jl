"""
    LineCableModels.EarthProps

Static physical descriptions of homogeneous and layered earth. Analysis
frequencies and frequency-dependent property formulations belong to
`LineCableModels.Engine`.
"""
module EarthProps

export EarthLayer, EarthModel

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using DataFrames
import ..LineCableModels: add!, validate

"""
$(TYPEDEF)

Store the static physical properties of one earth or air layer.

$(TYPEDFIELDS)
"""
struct EarthLayer{T <: Real}
    "Electrical resistivity \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Layer thickness \\[m\\]."
    thickness::T

    function EarthLayer{T}(rho::T, eps_r::T, mu_r::T, thickness::T) where {T <: Real}
        return validate(new{T}(rho, eps_r, mu_r, thickness))
    end
end

Base.eltype(::EarthLayer{T}) where {T} = T
Base.eltype(::Type{EarthLayer{T}}) where {T} = T

function validate(layer::EarthLayer)
    isnan(layer.rho) && throw(DomainError(layer.rho, "resistivity cannot be NaN"))
    layer.rho > zero(layer.rho) ||
        throw(DomainError(layer.rho, "resistivity must be positive"))
    isfinite(layer.eps_r) && layer.eps_r > zero(layer.eps_r) ||
        throw(DomainError(layer.eps_r, "relative permittivity must be positive and finite"))
    isfinite(layer.mu_r) && layer.mu_r > zero(layer.mu_r) ||
        throw(DomainError(layer.mu_r, "relative permeability must be positive and finite"))
    layer.thickness > zero(layer.thickness) ||
        throw(DomainError(layer.thickness, "thickness must be positive"))
    return layer
end

"""
$(TYPEDSIGNATURES)

Construct one static earth layer after floating and promoting its properties.
"""
function EarthLayer(rho::Real, eps_r::Real, mu_r::Real, thickness::Real)
    values = promote(float(rho), float(eps_r), float(mu_r), float(thickness))
    T = typeof(first(values))
    return EarthLayer{T}(values...)
end

function EarthLayer(rho::Real, eps_r::Real, mu_r::Real)
    values = promote(float(rho), float(eps_r), float(mu_r))
    T = typeof(first(values))
    return EarthLayer{T}(values..., T(Inf))
end

function Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer) where {T <: Real}
    return EarthLayer{T}(
        convert(T, layer.rho), convert(T, layer.eps_r),
        convert(T, layer.mu_r), convert(T, layer.thickness)
    )
end

Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer{T}) where {T <: Real} = layer

"""
$(TYPEDEF)

Store a static layered earth description. The first layer is always air.

$(TYPEDFIELDS)
"""
mutable struct EarthModel{T <: Real}
    "Whether earth interfaces are vertical rather than horizontal."
    vertical_layers::Bool
    "Static layers beginning with semi-infinite air."
    layers::Vector{EarthLayer{T}}

    function EarthModel{T}(
            vertical_layers::Bool,
            layers::Vector{EarthLayer{T}}
    ) where {T <: Real}
        return validate(new{T}(vertical_layers, layers))
    end
end

Base.eltype(::EarthModel{T}) where {T} = T
Base.eltype(::Type{EarthModel{T}}) where {T} = T

function _valid_layer_sequence(vertical::Bool, layers)
    length(layers) >= 2 || throw(ArgumentError("an earth model requires air and earth"))
    air = first(layers)
    isinf(air.rho) && isinf(air.thickness) || throw(ArgumentError(
        "the first layer must be semi-infinite air",
    ))
    vertical && !isinf(layers[2].thickness) &&
        throw(ArgumentError(
            "the first vertical earth layer must be semi-infinite",
        ))
    for index in 2:(length(layers) - 1)
        if isinf(layers[index].thickness)
            vertical && index == 2 && continue
            throw(ArgumentError("an infinite layer must be the final layer"))
        end
    end
    return layers
end

function validate(model::EarthModel)
    foreach(validate, model.layers)
    _valid_layer_sequence(model.vertical_layers, model.layers)
    return model
end

"""
$(TYPEDSIGNATURES)

Construct a static earth model. Frequencies are supplied later by a
`LineParametersProblem`.
"""
function EarthModel(
        rho::Real,
        eps_r::Real = one(float(rho)),
        mu_r::Real = one(float(rho));
        thickness::Real = oftype(float(rho), Inf),
        vertical_layers::Bool = false,
        air_layer::Union{Nothing, EarthLayer} = nothing
)
    T = promote_type(
        typeof(float(rho)), typeof(float(eps_r)), typeof(float(mu_r)),
        typeof(float(thickness)),
        air_layer === nothing ? typeof(float(rho)) : eltype(air_layer)
    )
    air = air_layer === nothing ?
          EarthLayer{T}(T(Inf), convert(T, 1), convert(T, 1), T(Inf)) :
          convert(EarthLayer{T}, air_layer)
    earth = EarthLayer(
        convert(T, rho), convert(T, eps_r), convert(T, mu_r),
        convert(T, thickness)
    )
    return validate(EarthModel{T}(
        vertical_layers, EarthLayer{T}[air, convert(EarthLayer{T}, earth)]
    ))
end

function Base.convert(::Type{EarthModel{T}}, model::EarthModel) where {T <: Real}
    return validate(EarthModel{T}(
        model.vertical_layers,
        EarthLayer{T}[convert(EarthLayer{T}, layer) for layer in model.layers]
    ))
end

Base.convert(::Type{EarthModel{T}}, model::EarthModel{T}) where {T <: Real} = model

function add!(model::EarthModel{T}, layer::EarthLayer{T}) where {T}
    candidate = EarthLayer{T}[model.layers; layer]
    _valid_layer_sequence(model.vertical_layers, candidate)
    push!(model.layers, layer)
    return validate(model)
end

function add!(model::EarthModel{T}, layer::EarthLayer{U}) where {T, U}
    throw(ArgumentError(
        "cannot add EarthLayer{$U} to EarthModel{$T}; explicitly convert the " *
        "complete model before mutation",
    ))
end

function add!(
        model::EarthModel{T},
        rho::Real,
        eps_r::Real,
        mu_r::Real;
        thickness::Real = T(Inf)
) where {T}
    if !(rho isa T && eps_r isa T && mu_r isa T && thickness isa T)
        throw(ArgumentError(
            "earth-layer scalar inputs must all be $T; construct EarthLayer{$T} " *
            "explicitly or convert the complete model before mutation",
        ))
    end
    return add!(model, EarthLayer{T}(rho, eps_r, mu_r, thickness))
end

include("dataframe.jl")
include("base.jl")

end
