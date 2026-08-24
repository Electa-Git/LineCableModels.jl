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

    for cable_position in system.cables
        push!(cable_ids, cable_position.design_data.cable_id)
        push!(horz_coords, cable_position.horz)
        push!(vert_coords, cable_position.vert)

        component_names = [comp.id for comp in cable_position.design_data.components]
        mapping_str = join(
            ["$(name): $(phase)"
             for (name, phase) in zip(component_names, cable_position.conn)],
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
