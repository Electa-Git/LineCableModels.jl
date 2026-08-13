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
	function EnforceLayer(; layer::Int = -1)
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

Builds a 2-layer (air + one enforced earth layer) data pack as three matrices
ρ, ε, μ of size (2 × n_freq), already converted to `T`.

# Returns
- `(ρ, ε, μ) :: (Matrix{T}, Matrix{T}, Matrix{T})`
  with row 1 = air, row 2 = enforced earth layer.
"""
function (f::EnforceLayer)(
	model::EarthModel,
	freq::AbstractVector{<:REALSCALAR},
	::Type{T},
) where {T <: REALSCALAR}

	nL = length(model.layers)
	nF = length(freq)

	layer_idx = f.layer == -1 ? nL : f.layer
	(2 <= layer_idx <= nL) || error(
		"Invalid layer index: $layer_idx. Model has $nL layers (including air). " *
		"Valid earth layer indices are 2:$nL.",
	)

	Lair = model.layers[1]
	Lsel = model.layers[layer_idx]

	ρ = Matrix{T}(undef, 2, nF)
	ε = similar(ρ)
	μ = similar(ρ)

	@inbounds for j in 1:nF
		ρ[1, j] = T(Lair.rho_g[j])
		ε[1, j] = T(Lair.eps_g[j])
		μ[1, j] = T(Lair.mu_g[j])

		ρ[2, j] = T(Lsel.rho_g[j])
		ε[2, j] = T(Lsel.eps_g[j])
		μ[2, j] = T(Lsel.mu_g[j])
	end

	return ρ, ε, μ
end

# """
# $(TYPEDSIGNATURES)

# Functor implementation for `EnforceLayer`.

# Takes a multi-layer `EarthModel` and returns a new two-layer model (air + one effective earth layer) based on the properties of the layer specified in the `EnforceLayer` instance.

# # Returns
# - A `Vector{EarthLayer}` containing two layers: the original air layer and the selected earth layer.
# """
# function (f::EnforceLayer)(
# 	model::EarthModel,
# 	freq::Vector{<:REALSCALAR},
# 	T::DataType,
# )
# 	num_layers = length(model.layers)

# 	# Determine the index of the layer to select
# 	layer_idx = f.layer == -1 ? num_layers : f.layer

# 	# Validate the chosen index
# 	if !(2 <= layer_idx <= num_layers)
# 		Base.error(
# 			"Invalid layer index: $layer_idx. The model only has $num_layers layers (including air). Valid earth layer indices are from 2 to $num_layers.",
# 		)
# 	end

# 	# The air layer is always the first layer in the original model
# 	air_layer = model.layers[1]

# 	# The enforced earth layer is the one at the selected index
# 	enforced_layer = model.layers[layer_idx]

# 	# Create a NamedTuple for the air layer with type-promoted property vectors
# 	air_data = (
# 		rho_g = T.(air_layer.rho_g),
# 		eps_g = T.(air_layer.eps_g),
# 		mu_g = T.(air_layer.mu_g),
# 	)

# 	# Create a NamedTuple for the enforced earth layer
# 	earth_data = (
# 		rho_g = T.(enforced_layer.rho_g),
# 		eps_g = T.(enforced_layer.eps_g),
# 		mu_g = T.(enforced_layer.mu_g),
# 	)

# 	# Return a new vector containing only these two layers
# 	return [air_data, earth_data]
# end