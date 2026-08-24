# Base property and display protocols owned by DataModel.
# RectStrands exposes its shape fields as part of its geometry.
function Base.getproperty(part::RectStrands, sym::Symbol)
    # Read fields from RectStrands itself before consulting its shape.
    if hasfield(typeof(part), sym)
        return getfield(part, sym)
    end

    # Route other properties to the shape payload.
    if hasfield(typeof(part), :shape)
        shape_payload = getfield(part, :shape)
        # Honour custom property access implemented by the shape.
        if hasproperty(shape_payload, sym)
            return getproperty(shape_payload, sym)
        end
    end

    # Use the standard missing-field error for unknown properties.
    return getfield(part, sym)
end

function Base.hasproperty(part::RectStrands, sym::Symbol)
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

function Base.propertynames(part::RectStrands, private::Bool = false)
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

Write the plain-text representation of a cable part to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `part`: Cable part to display.

# Returns

- `nothing` after writing to `io`.
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
            :shunt_conductance
        ]
    )

    println(io, "]")

    # Display material properties if available
    if hasproperty(part, :material_props)
        print(io, "└─ Material properties: [")
        _print_fields(
            io,
            part.material_props,
            [:rho, :eps_r, :mu_r, :alpha]
        )
        println(io, "]")
    end
end

"""
$(TYPEDSIGNATURES)

Write the plain-text representation of a conductor or insulator group to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `group`: Conductor or insulator group to display.

# Returns

- `nothing` after writing to `io`.
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
            :shunt_conductance
        ]
    )
    println(io, "]")

    # Display each layer as a tree entry.

    for (i, layer) in enumerate(group.layers)
        # Select the prefix from the layer position.
        prefix = i == length(group.layers) ? "└─" : "├─"
        # Print selected layer fields.
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
                :shunt_conductance
            ]
        )
        println(io, "]")
    end
end

"""
$(TYPEDSIGNATURES)

Write the available fields in `fields_to_show` as comma-separated `name=value`
pairs.

# Arguments

- `io`: Output stream.
- `obj`: Object to inspect.
- `fields_to_show`: Field names in display order.

# Keywords

- `sigdigits`: Significant digits used for numeric values. Default: `4`.

# Returns

- Number of fields written to `io`.
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
