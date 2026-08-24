
Base.eltype(::CableComponent{T}) where {T} = T
Base.eltype(::Type{CableComponent{T}}) where {T} = T

"""
$(TYPEDSIGNATURES)

Write the plain-text representation of a cable component to `io`.

# Arguments

- `io`: Output stream.
- `::MIME"text/plain"`: MIME type for plain text output.
- `component`: Cable component to display.

# Returns

- `nothing` after writing to `io`.
"""
function Base.show(io::IO, ::MIME"text/plain", component::CableComponent)
    # Calculate total number of parts across both groups
    total_parts = length(component.conductor_group.layers) +
                  length(component.insulator_group.layers)

    # Print header
    println(io, "$(total_parts)-element CableComponent \"$(component.id)\":")

    # Display conductor group parts in a tree structure
    print(io, "├─ $(length(component.conductor_group.layers))-element ConductorGroup: [")
    _print_fields(
        io,
        component.conductor_group,
        [:r_in, :r_ex, :cross_section, :resistance, :gmr]
    )
    println(io, "]")
    print(io, "│  ", "├─", " Effective properties: [")
    _print_fields(io, component.conductor_props, [:rho, :eps_r, :mu_r, :alpha])
    println(io, "]")

    for (i, part) in enumerate(component.conductor_group.layers)
        prefix = i == length(component.conductor_group.layers) ? "└───" : "├───"

        # Indent the cable-part entry.
        print(io, "│  ", prefix, " $(nameof(typeof(part))): [")

        # Print selected cable-part fields.
        _print_fields(
            io,
            part,
            [:r_in, :r_ex, :cross_section, :resistance, :gmr]
        )

        println(io, "]")
    end

    # Display insulator group parts
    if !isempty(component.insulator_group.layers)
        print(
            io,
            "└─ $(length(component.insulator_group.layers))-element InsulatorGroup: ["
        )
        _print_fields(
            io,
            component.insulator_group,
            [
                :r_in,
                :r_ex,
                :cross_section,
                :shunt_capacitance,
                :shunt_conductance
            ]
        )
        println(io, "]")
        print(io, "   ", "├─", " Effective properties: [")
        _print_fields(io, component.insulator_props, [:rho, :eps_r, :mu_r, :alpha])
        println(io, "]")
        for (i, part) in enumerate(component.insulator_group.layers)
            # Select the prefix from the part position.
            prefix = i == length(component.insulator_group.layers) ? "└───" : "├───"

            # Indent the cable-part entry.
            print(io, "   ", prefix, " $(nameof(typeof(part))): [")

            # Print selected cable-part fields.
            _print_fields(
                io,
                part,
                [
                    :r_in,
                    :r_ex,
                    :cross_section,
                    :shunt_capacitance,
                    :shunt_conductance
                ]
            )

            println(io, "]")
        end
    end
end
