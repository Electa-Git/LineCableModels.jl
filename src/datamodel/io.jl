# Full inheritance over composition madness
function Base.getproperty(part::AbstractCablePart, sym::Symbol)
	# Fast path: Is it a real field on the top-level struct? (r_in, gmr, shape)
	if hasfield(typeof(part), sym)
		return getfield(part, sym) # MUST use getfield here to prevent infinite recursion
	end

	# Fallback path: Route it to the shape payload (if the struct has a shape)
	if hasfield(typeof(part), :shape)
		shape_payload = getfield(part, :shape)
		# We use hasproperty on the shape just in case the shape itself has custom routing
		if hasproperty(shape_payload, sym)
			return getproperty(shape_payload, sym)
		end
	end

	# If it doesn't exist anywhere, fallback to standard getfield to throw the normal error
	return getfield(part, sym)
end

function Base.hasproperty(part::AbstractCablePart, sym::Symbol)
	# 1. Does it exist at the top level?
	if hasfield(typeof(part), sym)
		return true
	end

	# 2. Does it exist in the shape payload?
	if hasfield(typeof(part), :shape)
		return hasproperty(getfield(part, :shape), sym)
	end

	return false
end

function Base.propertynames(part::AbstractCablePart, private::Bool = false)
	top_fields = fieldnames(typeof(part))

	if hasfield(typeof(part), :shape)
		shape_payload = getfield(part, :shape)
		shape_fields = propertynames(shape_payload, private)

		# Merge and deduplicate the fields
		return Tuple(unique((top_fields..., shape_fields...)))
	end

	return top_fields
end

"""
$(TYPEDSIGNATURES)

Defines the display representation of an [`AbstractCablePart`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `part`: The [`AbstractCablePart`](@ref) instance to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the object.
"""
function Base.show(io::IO, ::MIME"text/plain", part::T) where {T <: AbstractCablePart}
	# Start output with type name
	print(io, "$(nameof(T)): [")

	# Use _print_fields to display all relevant fields
	_print_fields(
		io,
		part,
		[
			:r_in,
			:r_ex,
			:cross_section,
			:resistance,
			:gmr,
			:shunt_capacitance,
			:shunt_conductance,
		],
	)

	println(io, "]")

	# Display material properties if available
	if hasproperty(part, :material_props)
		print(io, "└─ Material properties: [")
		_print_fields(
			io,
			part.material_props,
			[:rho, :eps_r, :mu_r, :alpha],
		)
		println(io, "]")
	end
end


"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`ConductorGroup`](@ref) or [`InsulatorGroup`](@ref)objects for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `group`: The [`ConductorGroup`](@ref) or [`InsulatorGroup`](@ref) instance to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the object.
"""
function Base.show(io::IO, ::MIME"text/plain", group::Union{ConductorGroup, InsulatorGroup})

	print(io, "$(length(group.layers))-element $(nameof(typeof(group))): [")
	_print_fields(
		io,
		group,
		[
			:r_in,
			:r_ex,
			:cross_section,
			:resistance,
			:gmr,
			:shunt_capacitance,
			:shunt_conductance,
		],
	)
	println(io, "]")

	# Tree-like layer representation

	for (i, layer) in enumerate(group.layers)
		# Determine prefix based on whether it's the last layer
		prefix = i == length(group.layers) ? "└─" : "├─"
		# Print layer information with only selected fields
		print(io, prefix, "$(nameof(typeof(layer))): [")
		_print_fields(
			io,
			layer,
			[
				:r_in,
				:r_ex,
				:cross_section,
				:resistance,
				:gmr,
				:shunt_capacitance,
				:shunt_conductance,
			],
		)
		println(io, "]")
	end
end

"""
$(TYPEDSIGNATURES)

Print the specified fields of an object in a compact format.

# Arguments
- `io`: The output stream.
- `obj`: The object whose fields will be displayed.
- `fields_to_show`: Vector of field names (as Symbols) to display.
- `sigdigits`: Number of significant digits for rounding numeric values.

# Returns

- Number of fields that were actually displayed.
"""
function _print_fields(io::IO, obj, fields_to_show::Vector{Symbol}; sigdigits::Int = 4)
	displayed_fields = 0
	for field in fields_to_show
		if hasproperty(obj, field)
			value = getproperty(obj, field)
			# Skip NaN values
			if value isa Number && isnan(value)
				continue
			end
			# Add comma if not the first item
			if displayed_fields > 0
				print(io, ", ")
			end
			# Format numbers with rounding
			if value isa Number
				print(io, "$field=$(round(value, sigdigits=sigdigits))")
			else
				print(io, "$field=$value")
			end
			displayed_fields += 1
		end
	end
	return displayed_fields
end
