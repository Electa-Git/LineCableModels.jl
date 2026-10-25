"""
    LineCableModels.Earth

Describe the earth below the air as an ordered stack of layers, semi-infinite air first.
"""
module Earth

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: validate

export EarthLayer, EarthModel, homogeneous

"""
$(TYPEDEF)

Store the static properties of one earth or air layer.

$(TYPEDFIELDS)
"""
struct EarthLayer{T <: Real}
    "Electrical resistivity \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Thickness, infinite for a semi-infinite layer \\[m\\]."
    thickness::T
    function EarthLayer{T}(rho::T, eps_r::T, mu_r::T, thickness::T) where {T <: Real}
        return validate(new{T}(rho, eps_r, mu_r, thickness))
    end
end

"""
$(TYPEDSIGNATURES)

Construct one layer after floating and promoting its properties. The arguments are the fields
of [`EarthLayer`](@ref).

# Errors

- Throws `DomainError` when a property is not positive, or when `eps_r` or `mu_r` is not
  finite.
"""
function EarthLayer(rho::Real, eps_r::Real, mu_r::Real, thickness::Real)
    values = map(float, promote(rho, eps_r, mu_r, thickness))
    return EarthLayer{eltype(values)}(values...)
end

"""
$(TYPEDSIGNATURES)

Construct one semi-infinite layer: its thickness is infinite in the floating type that its
properties promote to.
"""
function EarthLayer(rho::Real, eps_r::Real, mu_r::Real)
    values = map(float, promote(rho, eps_r, mu_r))
    T = eltype(values)
    return EarthLayer{T}(values..., T(Inf))
end

function Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer) where {T <: Real}
    return EarthLayer{T}(convert(T, layer.rho), convert(T, layer.eps_r),
        convert(T, layer.mu_r), convert(T, layer.thickness))
end

function validate(layer::EarthLayer)
    isnan(layer.rho) && throw(DomainError(layer.rho, "EarthLayer.rho must not be NaN"))
    layer.rho > zero(layer.rho) ||
        throw(DomainError(layer.rho, "EarthLayer.rho must be positive"))
    isfinite(layer.eps_r) && layer.eps_r > zero(layer.eps_r) ||
        throw(DomainError(layer.eps_r, "EarthLayer.eps_r must be positive and finite"))
    isfinite(layer.mu_r) && layer.mu_r > zero(layer.mu_r) ||
        throw(DomainError(layer.mu_r, "EarthLayer.mu_r must be positive and finite"))
    layer.thickness > zero(layer.thickness) ||
        throw(DomainError(layer.thickness, "EarthLayer.thickness must be positive"))
    return layer
end

"""
$(TYPEDEF)

Store an earth description: semi-infinite air above the earth layers. This release admits one
semi-infinite earth layer.

$(TYPEDFIELDS)
"""
struct EarthModel{T <: Real, N}
    "Layers from the top down, semi-infinite air first."
    layers::NTuple{N, EarthLayer{T}}
    function EarthModel{T, N}(layers::NTuple{N, EarthLayer{T}}) where {T <: Real, N}
        return validate(new{T, N}(layers))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an earth model from its layers, semi-infinite air first.

# Errors

- Throws `ArgumentError` unless the layers are semi-infinite air and one semi-infinite earth
  layer.
"""
function EarthModel(layers::Tuple{EarthLayer{T}, Vararg{EarthLayer{T}}}) where {T <: Real}
    return EarthModel{T, length(layers)}(layers)
end

Base.eltype(::EarthModel{T}) where {T} = T

function Base.convert(::Type{EarthModel{T}}, model::EarthModel) where {T <: Real}
    return EarthModel(map(layer -> convert(EarthLayer{T}, layer), model.layers))
end

function validate(model::EarthModel)
    length(model.layers) == 2 || throw(ArgumentError(
        "EarthModel.layers must be air and one earth layer; received " *
        "$(length(model.layers)) layers"))
    air, earth = model.layers
    isinf(air.rho) && isinf(air.thickness) ||
        throw(ArgumentError("EarthModel.layers[1] must be semi-infinite air"))
    isinf(earth.thickness) ||
        throw(ArgumentError("EarthModel.layers[2] must be semi-infinite"))
    return model
end

"""
$(TYPEDSIGNATURES)

Build the [`EarthModel`](@ref) of semi-infinite air, of relative permittivity and permeability
1, above one semi-infinite earth layer.

# Keywords

- `rho`: electrical resistivity of the earth \\[Ω·m\\].
- `eps_r`: relative permittivity of the earth \\[dimensionless\\].
- `mu_r`: relative permeability of the earth \\[dimensionless\\].

# Returns

- An `EarthModel{T, 2}`, where `T` is the floating type that the earth's properties promote to.
"""
function homogeneous(; rho::Real, eps_r::Real = 1, mu_r::Real = 1)
    earth = EarthLayer(rho, eps_r, mu_r)
    T = eltype(earth.rho)
    return EarthModel((EarthLayer(T(Inf), one(T), one(T)), earth))
end

end
