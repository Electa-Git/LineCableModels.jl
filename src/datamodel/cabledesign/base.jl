import Base: get, show, delete!, length, setindex!, iterate, keys, values, haskey, getindex, eltype

eltype(::CableDesign{T}) where {T} = T
eltype(::Type{CableDesign{T}}) where {T} = T


"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`CableDesign`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `design`: The [`CableDesign`](@ref) object to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the object.
"""
function show(io::IO, ::MIME"text/plain", design::CableDesign)
    # Print header with cable ID and count of components
    print(io, "$(length(design.components))-element CableDesign \"$(design.cable_id)\"")

    # Add nominal values if available
    nominal_values = []
    if design.nominal_data.resistance !== nothing
        push!(
            nominal_values,
            "resistance=$(round(design.nominal_data.resistance, sigdigits=4))",
        )
    end
    if design.nominal_data.inductance !== nothing
        push!(
            nominal_values,
            "inductance=$(round(design.nominal_data.inductance, sigdigits=4))",
        )
    end
    if design.nominal_data.capacitance !== nothing
        push!(
            nominal_values,
            "capacitance=$(round(design.nominal_data.capacitance, sigdigits=4))",
        )
    end

    if !isempty(nominal_values)
        print(io, ", with nominal values: [", join(nominal_values, ", "), "]")
    end
    println(io)

    # For each component, display its properties in a tree structure
    for (i, component) in enumerate(design.components)
        # Determine if this is the last component
        is_last_component = i == length(design.components)

        # Determine component prefix and continuation line
        component_prefix = is_last_component ? "└─" : "├─"
        continuation_line = is_last_component ? "   " : "│  "

        # Print component name and header
        println(io, component_prefix, " Component \"", component.id, "\":")

        # Display conductor group with combined properties
        print(io, continuation_line, "├─ ConductorGroup: [")

        # Combine properties from conductor_group and conductor_props
        conductor_props = [
            "radius_in" => component.conductor_group.radius_in,
            "radius_ext" => component.conductor_group.radius_ext,
            "rho" => component.conductor_props.rho,
            "eps_r" => component.conductor_props.eps_r,
            "mu_r" => component.conductor_props.mu_r,
            "alpha" => component.conductor_props.alpha,
        ]

        # Display combined conductor properties
        displayed_fields = 0
        for (field, value) in conductor_props
            if !(value isa Number && isnan(value))
                if displayed_fields > 0
                    print(io, ", ")
                end
                print(io, "$field=$(round(value, sigdigits=4))")
                displayed_fields += 1
            end
        end
        println(io, "]")

        # Display insulator group with combined properties
        print(io, continuation_line, "└─ InsulatorGroup: [")

        # Combine properties from insulator_group and insulator_props
        insulator_props = [
            "radius_in" => component.insulator_group.radius_in,
            "radius_ext" => component.insulator_group.radius_ext,
            "rho" => component.insulator_props.rho,
            "eps_r" => component.insulator_props.eps_r,
            "mu_r" => component.insulator_props.mu_r,
            "alpha" => component.insulator_props.alpha,
        ]

        # Display combined insulator properties
        displayed_fields = 0
        for (field, value) in insulator_props
            if !(value isa Number && isnan(value))
                if displayed_fields > 0
                    print(io, ", ")
                end
                print(io, "$field=$(round(value, sigdigits=4))")
                displayed_fields += 1
            end
        end
        println(io, "]")
    end
end