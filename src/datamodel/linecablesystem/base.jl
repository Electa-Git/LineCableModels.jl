
Base.eltype(::CablePosition{T}) where {T} = T
Base.eltype(::Type{CablePosition{T}}) where {T} = T
Base.eltype(::LineCableSystem{T}) where {T} = T
Base.eltype(::Type{LineCableSystem{T}}) where {T} = T

function Base.show(io::IO, ::MIME"text/plain", system::LineCableSystem)
    # Print top level info
    println(
        io,
        "LineCableSystem \"$(system.system_id)\": [line_length=$(system.line_length), num_cables=$(system.num_cables), num_phases=$(system.num_phases)]",
    )

    # Print cable definitions
    println(io, "└─ $(length(system.cables))-element CablePosition:")

    # Display each cable definition
    for (i, cable_position) in enumerate(system.cables)
        # Cable prefix
        prefix = i == length(system.cables) ? "   └─" : "   ├─"

        # Format connections as a string
        components = [comp.id for comp in cable_position.design_data.components]
        conn_str = join(
            ["$(comp)→$(phase)" for (comp, phase) in zip(components, cable_position.conn)],
            ", ",
        )

        # Print cable info
        println(
            io,
            "$(prefix) CableDesign \"$(cable_position.design_data.cable_id)\": [horz=$(round(cable_position.horz, sigdigits=4)), vert=$(round(cable_position.vert, sigdigits=4)), conn=($(conn_str))]",
        )
    end
end