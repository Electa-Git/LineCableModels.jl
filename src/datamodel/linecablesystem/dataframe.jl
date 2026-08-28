import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a line-cable system to one row per positioned cable.

# Arguments

- `system`: Materialised line-cable system.

# Returns

- A `DataFrame` with `cable_id`, horizontal position `horz` \\[m\\], vertical
  position `vert` \\[m\\], and the component-to-phase mapping
  `phase_mapping`.

"""
function DataFrame(system::LineCableSystem)::DataFrame
    cable_ids = String[]
    horz_coords = Number[]
    vert_coords = Number[]
    mappings = String[]

    for (design, position, connections) in zip(
            system.designs,
            system.positions,
            system.connections
    )
        push!(cable_ids, design.cable_id)
        push!(horz_coords, position.x)
        push!(vert_coords, position.y)

        mapping_str = join(
            ["$(name): $(phase)"
             for (name, phase) in zip(design.terminal_order, connections)],
            ", "
        )
        push!(mappings, mapping_str)
    end
    data = DataFrame(
        cable_id = cable_ids,
        horz = horz_coords,
        vert = vert_coords,
        phase_mapping = mappings
    )
    return data
end
