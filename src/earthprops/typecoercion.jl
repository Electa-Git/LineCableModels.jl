"""
$(TYPEDSIGNATURES)

Converts an `EarthModel{S}` to `EarthModel{T}` by reconstructing the model with
all layers coerced to the new scalar type `T`. Layer conversion is delegated to
[`coerce_to_T(::EarthLayer, ::Type)`](@ref), and non-numeric metadata are
forwarded unchanged.

# Arguments

- `model`: Source Earth model \\[dimensionless\\].
- `::Type{T}`: Target element type for numeric fields \\[dimensionless\\].

# Returns

- `EarthModel{T}` rebuilt with each layer and numeric payload converted to `T`.

# Examples

```julia
m64 = $(FUNCTIONNAME)(model, Float64)
mM  = $(FUNCTIONNAME)(model, Measurement{Float64})
```

# See also

- [`coerce_to_T`](@ref)
- [`resolve_T`](@ref)
"""
function coerce_to_T(model::EarthModel, ::Type{T}) where {T}
    # 1. Coerce all existing layers recursively
    new_layers = [coerce_to_T(layer, T) for layer in model.layers]

    # 2. Use the inner constructor to build the new, promoted model
    return EarthModel{T}(
        model.freq_dependence,
        model.vertical_layers,
        new_layers
    )
end

"""
$(TYPEDSIGNATURES)

Converts an `EarthLayer{S}` to `EarthLayer{T}` by coercing each stored field to
the target element type `T` and rebuilding the layer via its inner constructor.
Scalar and array fields are converted using the generic [`coerce_to_T`](@ref)
machinery.

# Arguments

- `layer`: Source Earth layer \\[dimensionless\\].
- `::Type{T}`: Target element type for numeric fields \\[dimensionless\\].

# Returns

- `EarthLayer{T}` with all numeric state converted to `T`.

# Examples

```julia
ℓ64 = $(FUNCTIONNAME)(layer, Float64)
ℓM  = $(FUNCTIONNAME)(layer, Measurement{Float64})
```

# See also

- [`coerce_to_T`](@ref)
- [`resolve_T`](@ref)
"""
function coerce_to_T(layer::EarthLayer, ::Type{T}) where {T}
    # Reconstruct the layer using the correct internal constructor.
    # The existing coerce_to_T methods for scalars and arrays will be
    # dispatched automatically for each field.
    return EarthLayer{T}(
        coerce_to_T(layer.base_rho_g, T),
        coerce_to_T(layer.base_epsr_g, T),
        coerce_to_T(layer.base_mur_g, T),
        coerce_to_T(layer.t, T),
        coerce_to_T(layer.rho_g, T),
        coerce_to_T(layer.eps_g, T),
        coerce_to_T(layer.mu_g, T)
    )
end

