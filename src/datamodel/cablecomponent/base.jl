import Base: get, show, delete!, length, setindex!, iterate, keys, values, haskey, getindex, eltype

eltype(::CableComponent{T}) where {T} = T
eltype(::Type{CableComponent{T}}) where {T} = T

"""
$(TYPEDSIGNATURES)

Defines the display representation of a [`CableComponent`](@ref) object for REPL or text output.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `component`: The [`CableComponent`](@ref) object to be displayed.

# Returns

- Nothing. Modifies `io` by writing text representation of the object.
"""
function show(io::IO, ::MIME"text/plain", component::CableComponent)
    # Calculate total number of parts across both groups
    total_parts =
        length(component.conductor_group.layers) + length(component.insulator_group.layers)

    # Print header
    println(io, "$(total_parts)-element CableComponent \"$(component.id)\":")

    # Display conductor group parts in a tree structure
    print(io, "├─ $(length(component.conductor_group.layers))-element ConductorGroup: [")
    _print_fields(
        io,
        component.conductor_group,
        [:radius_in, :radius_ext, :cross_section, :resistance, :gmr],
    )
    println(io, "]")
    print(io, "│  ", "├─", " Effective properties: [")
    _print_fields(io, component.conductor_props, [:rho, :eps_r, :mu_r, :alpha])
    println(io, "]")

    for (i, part) in enumerate(component.conductor_group.layers)

        prefix = i == length(component.conductor_group.layers) ? "└───" : "├───"

        # Print part information with proper indentation
        print(io, "│  ", prefix, " $(nameof(typeof(part))): [")

        # Print each field with proper formatting
        _print_fields(
            io,
            part,
            [:radius_in, :radius_ext, :cross_section, :resistance, :gmr],
        )

        println(io, "]")
    end

    # Display insulator group parts
    if !isempty(component.insulator_group.layers)
        print(
            io,
            "└─ $(length(component.insulator_group.layers))-element InsulatorGroup: [",
        )
        _print_fields(
            io,
            component.insulator_group,
            [
                :radius_in,
                :radius_ext,
                :cross_section,
                :shunt_capacitance,
                :shunt_conductance,
            ],
        )
        println(io, "]")
        print(io, "   ", "├─", " Effective properties: [")
        _print_fields(io, component.insulator_props, [:rho, :eps_r, :mu_r, :alpha])
        println(io, "]")
        for (i, part) in enumerate(component.insulator_group.layers)
            # Determine prefix based on whether it's the last part
            prefix = i == length(component.insulator_group.layers) ? "└───" : "├───"

            # Print part information with proper indentation
            print(io, "   ", prefix, " $(nameof(typeof(part))): [")

            # Print each field with proper formatting
            _print_fields(
                io,
                part,
                [
                    :radius_in,
                    :radius_ext,
                    :cross_section,
                    :shunt_capacitance,
                    :shunt_conductance,
                ],
            )

            println(io, "]")
        end
    end
end