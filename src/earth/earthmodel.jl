"""
$(TYPEDEF)

Store an immutable static layered-earth description. The first layer is always
air, and `layers` is a read-only ordered tuple.

$(TYPEDFIELDS)
"""
struct EarthModel{T <: Real, N} <: AbstractEarthModel
    "Whether earth interfaces are vertical rather than horizontal."
    vertical_layers::Bool
    "Read-only static layers beginning with semi-infinite air."
    layers::NTuple{N, EarthLayer{T}}

    function EarthModel{T, N}(
            vertical_layers::Bool,
            layers::NTuple{N, EarthLayer{T}}
    ) where {T <: Real, N}
        return validate(new{T, N}(vertical_layers, layers))
    end
end

function EarthModel{T}(
        vertical_layers::Bool,
        layers::NTuple{N, EarthLayer{T}}
) where {T <: Real, N}
    return EarthModel{T, N}(vertical_layers, layers)
end

Base.eltype(::EarthModel{T}) where {T} = T
Base.eltype(::Type{<:EarthModel{T}}) where {T} = T

function validate(model::EarthModel)
    length(model.layers) >= 2 || throw(ArgumentError(
        "EarthModel.layers must contain air and at least one earth layer; " *
        "received $(length(model.layers)) layers"
    ))
    air = first(model.layers)
    isinf(air.rho) && isinf(air.thickness) || throw(ArgumentError(
        "EarthModel.layers[1] must be semi-infinite air",
    ))
    model.vertical_layers && !isinf(model.layers[2].thickness) &&
        throw(ArgumentError(
            "EarthModel.layers[2].thickness must be infinite for vertical layers",
        ))
    for index in 2:(length(model.layers) - 1)
        if isinf(model.layers[index].thickness)
            model.vertical_layers && index == 2 && continue
            throw(ArgumentError(
                "EarthModel.layers[$index].thickness must be finite because " *
                "only the final horizontal layer may be semi-infinite"
            ))
        end
    end
    return model
end

"""
$(TYPEDSIGNATURES)

Construct a complete immutable static earth model. Frequencies are supplied
later by a `LineParametersProblem`.
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
    return EarthModel{T}(
        vertical_layers, (air, convert(EarthLayer{T}, earth))
    )
end

function _earth_model(layers, vertical_layers::Bool, air_layer::Union{Nothing, EarthLayer})
    declared = _declared_layers(layers)
    T = promote_type(
        eltype.(declared)...,
        air_layer === nothing ? eltype(first(declared)) : eltype(air_layer)
    )
    air = air_layer === nothing ?
          EarthLayer{T}(T(Inf), convert(T, 1), convert(T, 1), T(Inf)) :
          convert(EarthLayer{T}, air_layer)
    earth = map(layer -> convert(EarthLayer{T}, layer), declared)
    return EarthModel{T}(vertical_layers, (air, earth...))
end
function _earth_model(layers, vertical_layers::Bool, air_layer)
    throw(ArgumentError("air_layer must be nothing or a completed EarthLayer"))
end
function _earth_model(layers, vertical_layers, air_layer)
    throw(ArgumentError("vertical_layers must be true or false"))
end

# Declared earth layers are one layer or a nonempty collection of layers.
_declared_layers(layer::EarthLayer) = (layer,)
function _declared_layers(layers::Union{Tuple, AbstractVector})
    foreach(_declared_layer, layers)
    declared = Tuple(layers)
    isempty(declared) && throw(ArgumentError(
        "an earth model requires at least one earth layer"
    ))
    return declared
end
function _declared_layers(layers)
    throw(ArgumentError(
        "earth layers must be an EarthLayer or a nonempty layer collection"
    ))
end
_declared_layer(layer::EarthLayer) = layer
function _declared_layer(layer)
    throw(ArgumentError("earth layers must contain completed EarthLayer objects"))
end

"""
$(TYPEDSIGNATURES)

Build a complete immutable static earth model from one or more completed earth
layers. A semi-infinite air layer is prepended unless `air_layer` is supplied.

# Arguments

- `layers`: one `EarthLayer` or an ordered collection. Horizontal models are
  ordered from the earth surface downward.

# Keywords

- `vertical_layers=false`: whether earth interfaces are vertical.
- `air_layer=nothing`: optional explicit semi-infinite air layer.
- `combine=:product`: gridspace composition rule.

# Returns

- An `EarthModel`, or a `Gridspace{EarthModel}` when an explicit finite source
  is supplied.
"""
function build(
        ::Type{EarthModel},
        layers;
        vertical_layers = false,
        air_layer = nothing,
        combine::Symbol = :product
)
    values = (layers, vertical_layers, air_layer)
    return parameterize(EarthModel, _earth_model, values; combine)
end

function Base.convert(
        ::Type{EarthModel{T}}, model::EarthModel{U, N}
) where {T <: Real, U <: Real, N}
    return EarthModel{T}(
        model.vertical_layers,
        ntuple(index -> convert(EarthLayer{T}, model.layers[index]), N)
    )
end

Base.convert(::Type{EarthModel{T}}, model::EarthModel{T}) where {T <: Real} = model

struct HomogeneousBuilder <: Builder{EarthModel} end

"""
```julia
homogeneous(; rho, eps_r = nothing, mu_r = nothing, thickness = nothing,
    vertical_layers = false, air_layer = nothing, combine = :product)
```

Declare a homogeneous earth model through [`layer`](@ref) and
`build(EarthModel, ...)`. Semi-infinite air is implicit.

# Keywords

- `rho`: electrical resistivity \\[Ω·m\\].
- `eps_r=nothing`: relative permittivity \\[dimensionless\\]. `nothing`
  selects unity in the resistivity scalar type.
- `mu_r=nothing`: relative permeability \\[dimensionless\\]. `nothing`
  selects unity in the resistivity scalar type.
- `thickness=nothing`: earth-layer thickness \\[m\\]. `nothing` selects a
  semi-infinite earth layer.
- `vertical_layers=false`: whether earth interfaces are vertical.
- `air_layer=nothing`: optional explicit semi-infinite air layer.
- `combine=:product`: composition of the sources when an argument or a keyword is a `Grid`
  or a `Gridspace`, `:product` or `:zip`.

# Returns

- An `EarthModel`, or a `Gridspace{EarthModel}` when a keyword varies.
"""
const homogeneous = HomogeneousBuilder()

function build(::HomogeneousBuilder; rho, eps_r = nothing, mu_r = nothing,
        thickness = nothing, vertical_layers = false, air_layer = nothing)
    return build(EarthModel, build(layer; rho, eps_r, mu_r, thickness); vertical_layers,
        air_layer)
end
