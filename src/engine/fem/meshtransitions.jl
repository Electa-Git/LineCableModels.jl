"""
$(TYPEDEF)

Defines a mesh transition region for improved mesh quality in earth/air regions around cable systems.

$(TYPEDFIELDS)
"""
struct MeshTransition
	"Center coordinates (x, y) [m]"
	center::Tuple{Float64, Float64}
	"Minimum radius (must be ≥ bounding radius of cables) [m]"
	r_min::Float64
	"Maximum radius [m]"
	r_max::Float64
	"Minimum mesh size factor at r_min [m]"
	mesh_factor_min::Float64
	"Maximum mesh size factor at r_max [m]"
	mesh_factor_max::Float64
	"Number of transition regions [dimensionless]"
	n_regions::Int
	"Earth layer index (1=air, 2+=earth layers from top to bottom, nothing=auto-detect)"
	earth_layer::Union{Int, Nothing}

	function MeshTransition(
		center,
		r_min,
		r_max,
		mesh_factor_min,
		mesh_factor_max,
		n_regions,
		earth_layer,
	)
		# Basic validation
		r_min >= 0 || Base.error("r_min must be greater than or equal to 0")
		r_max > r_min || Base.error("r_max must be greater than r_min")
		mesh_factor_min > 0 || Base.error("mesh_factor_min must be positive")
		mesh_factor_max <= 1 ||
			Base.error("mesh_factor_max must be smaller than or equal to 1")
		mesh_factor_max > mesh_factor_min ||
			Base.error("mesh_factor_max must be > mesh_factor_min")
		n_regions >= 1 || Base.error("n_regions must be at least 1")

		# Validate earth_layer if provided
		if !isnothing(earth_layer)
			earth_layer >= 1 ||
				Base.error("earth_layer must be >= 1 (1=air, 2+=earth layers)")
		end

		new(center, r_min, r_max, mesh_factor_min, mesh_factor_max, n_regions, earth_layer)
	end
end

# Convenience constructor
function MeshTransition(
	cable_system::LineCableSystem,
	cable_indices::Vector{Int};
	r_min::Number,
	r_length::Number,
	mesh_factor_min::Number,
	mesh_factor_max::Number,
	n_regions::Int = 3,
	earth_layer::Union{Int, Nothing} = nothing,
)
	(r_min, r_length, mesh_factor_min, mesh_factor_max) =
		to_nominal.((r_min, r_length, mesh_factor_min, mesh_factor_max))

	# Validate cable indices
	all(1 <= idx <= length(cable_system.cables) for idx in cable_indices) ||
		Base.error("Cable indices out of bounds")

	isempty(cable_indices) && Base.error("Cable indices cannot be empty")

	# Get centroid and bounding radius
	cx, cy, bounding_radius, _ =
		to_nominal.(get_system_centroid(cable_system, cable_indices))

	# Calculate parameters
	if r_min < bounding_radius
		@warn "r_min ($r_min m) is smaller than bounding radius ($bounding_radius m). Adjusting r_min to match."
		r_min = bounding_radius
	end

	r_max = r_min + r_length

	# Auto-detect layer if not specified
	if isnothing(earth_layer)
		# Simple detection: y >= 0 is air (layer 1), y < 0 is first earth layer (layer 2)
		earth_layer = cy >= 0 ? 1 : 2
		@debug "Auto-detected earth_layer=$earth_layer for transition at ($cx, $cy)"
	end

	# Validate no surface crossing for underground transitions
	if earth_layer > 1 && cy + r_max > 0
		Base.error(
			"Transition region would cross earth surface (y=0). Reduce r_length or use separate transition regions.",
		)
	end

	return MeshTransition(
		(cx, cy),
		r_min,
		r_max,
		mesh_factor_min,
		mesh_factor_max,
		n_regions,
		earth_layer,
	)
end
