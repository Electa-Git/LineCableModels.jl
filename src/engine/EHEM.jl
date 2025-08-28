"""
	LineCableModels.Engine.EHEM

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module EHEM

# Export public API
export EnforceLayer

# Load common dependencies
using ...LineCableModels
include("../utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using ...Utils
import ...LineCableModels: get_description
import ..Engine: AbstractEHEMFormulation

"""
$(TYPEDEF)

An EHEM formulation that creates a homogeneous earth model by enforcing the properties of a single, specified layer from a multi-layer model.

# Attributes
$(TYPEDFIELDS)
"""
struct EnforceLayer <: AbstractEHEMFormulation
    "Index of the earth layer to enforce. `-1` selects the bottommost layer."
    layer::Int

    @doc """
    $(TYPEDSIGNATURES)

    Constructs an `EnforceLayer` instance.

    # Arguments
    - `layer::Int`: The index of the layer to enforce.
        - `-1` (default): Enforces the properties of the bottommost earth layer.
        - `2`: Enforces the properties of the topmost earth layer (the one directly below the air).
        - `> 2`: Enforces the properties of a specific layer by its index.
    """
    function EnforceLayer(; layer::Int=-1)
        @assert (layer == -1 || layer >= 2) "Invalid layer index. Must be -1 (bottommost) or >= 2."
        new(layer)
    end
end

function get_description(f::EnforceLayer)
    if f.layer == -1
        return "Assume bottom layer"
    elseif f.layer == 2
        return "Assume top earth layer"
    else
        return "Assume layer $(f.layer)"
    end
end

"""
$(TYPEDSIGNATURES)

Functor implementation for `EnforceLayer`.

Takes a multi-layer `EarthModel` and returns a new two-layer model (air + one effective earth layer) based on the properties of the layer specified in the `EnforceLayer` instance.

# Returns
- A `Vector{EarthLayer}` containing two layers: the original air layer and the selected earth layer.
"""
function (f::EnforceLayer)(model::EarthModel, frequencies::Vector{<:REALSCALAR}, T::DataType)
    num_layers = length(model.layers)

    # Determine the index of the layer to select
    layer_idx = f.layer == -1 ? num_layers : f.layer

    # Validate the chosen index
    if !(2 <= layer_idx <= num_layers)
        error("Invalid layer index: $layer_idx. The model only has $num_layers layers (including air). Valid earth layer indices are from 2 to $num_layers.")
    end

    # The air layer is always the first layer in the original model
    air_layer = model.layers[1]

    # The enforced earth layer is the one at the selected index
    enforced_layer = model.layers[layer_idx]

    # Create a NamedTuple for the air layer with type-promoted property vectors
    air_data = (
        rho_g=T.(air_layer.rho_g),
        eps_g=T.(air_layer.eps_g),
        mu_g=T.(air_layer.mu_g)
    )

    # Create a NamedTuple for the enforced earth layer
    earth_data = (
        rho_g=T.(enforced_layer.rho_g),
        eps_g=T.(enforced_layer.eps_g),
        mu_g=T.(enforced_layer.mu_g)
    )

    # Return a new vector containing only these two layers
    return [air_data, earth_data]
end

# """
# $(TYPEDEF)

# Represents a homogeneous earth model defined using the properties of a specific earth layer, with atttribute:

# $(TYPEDFIELDS)
# """
# struct EnforceLayer <: AbstractEHEMFormulation
#     "Index of the enforced earth layer."
#     layer::Int

#     @doc """
#     $(TYPEDSIGNATURES)

#     Constructs an [`EnforceLayer`](@ref) instance, meaning that the properties of the specified layer are extended to the entire semi-infinite earth.

#     # Arguments

#     - `layer`: Integer specifying the layer to be enforced as the equivalent homogeneous layer.
#     - `-1` selects the bottommost layer.
#     - `2` selects the topmost earth layer (excluding air).
#     - Any integer `≥ 2` selects a specific layer.

#     # Returns

#     - An [`EnforceLayer`](@ref) instance with the specified layer index.

#     # Examples

#     ```julia
#     # Enforce the bottommost layer
#     bottom_layer = $(FUNCTIONNAME)(-1)
#     println(get_description(bottom_layer)) # Output: "Bottom layer"

#     # Enforce the topmost earth layer
#     first_layer = $(FUNCTIONNAME)(2)
#     println(get_description(first_layer)) # Output: "Top layer"

#     # Enforce a specific layer
#     specific_layer = $(FUNCTIONNAME)(3)
#     println(get_description(specific_layer)) # Output: "Layer 3"
#     ```
#     """
#     function EnforceLayer(layer::Int)
#         @assert (layer == -1 || layer >= 2) "Invalid earth layer choice."
#         return new(layer)
#     end
# end

# function get_description(formulation::EnforceLayer)
#     if formulation.layer == -1
#         return "Bottom layer"
#     elseif formulation.layer == 2
#         return "Top layer"
#     else
#         return "Layer $(formulation.layer)"
#     end
# end

# """
# $(TYPEDSIGNATURES)

# Computes the effective homogeneous earth model (EHEM) properties for an [`EarthModel`](@ref), overriding the layered model with the properties of the layer defined in [`EnforceLayer`](@ref).

# # Arguments

# - `model`: Instance of [`EarthModel`](@ref) for which effective properties are computed.
# - `frequencies`: Vector of frequency values \\[Hz\\].
# - `formulation`: Instance of [`AbstractEHEMFormulation`](@ref) specifying how the effective properties should be determined.

# # Returns

# - Modifies `model` in place by updating `rho_eff`, `eps_eff`, and `mu_eff` with the corresponding values.

# # Notes

# - If `formulation` is an [`EnforceLayer`](@ref), the effective properties (`rho_eff`, `eps_eff`, `mu_eff`) are **directly assigned** from the specified layer.
# - If `layer_idx = -1`, the **last** layer in `model.layers` is used.
# - If `layer_idx < 2` or `layer_idx > length(model.layers)`, an error is raised.

# # Examples

# ```julia
# frequencies = [1e3, 1e4, 1e5]
# earth_model = EarthModel(frequencies, 100, 10, 1, t=5)
# earth_model.AbstractEHEMFormulation = EnforceLayer(-1)  # Enforce the last layer as the effective
# $(FUNCTIONNAME)(earth_model, frequencies, EnforceLayer(-1))

# println(earth_model.rho_eff) # Should match the last layer rho_g = 100
# println(earth_model.eps_eff) # Should match the last layer eps_g = 10*ε₀
# println(earth_model.mu_eff)  # Should match the last layer mu_g = 1*μ₀
# ```

# # See also

# - [`EarthModel`](@ref)
# - [`AbstractEHEMFormulation`](@ref)
# """
# function calc_ehem_properties(
#     model::EarthModel,
#     frequencies::Vector{<:T},
#     formulation::AbstractEHEMFormulation=EnforceLayer(-1),
# ) where {T<:REALSCALAR}
#     layer_idx = formulation.layer

#     if layer_idx == -1
#         layer_idx = length(model.layers)  # Use last layer
#     elseif layer_idx < 2 || layer_idx > length(model.layers)
#         error(
#             "Invalid layer index: $layer_idx. Must be between 2 and $(length(model.layers))",
#         )
#     end

#     layer = model.layers[layer_idx]

#     # Just reference the selected layer
#     rho_eff = layer.rho_g
#     eps_eff = layer.eps_g
#     mu_eff = layer.mu_g

#     return (rho_g=rho_eff, eps_g=eps_eff, mu_g=mu_eff)
# end

end # module EHEM