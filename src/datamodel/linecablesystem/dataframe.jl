import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Generates a summary DataFrame for cable positions and phase mappings within a [`LineCableSystem`](@ref).

# Arguments

- `system`: A [`LineCableSystem`](@ref) object containing the cable definitions and their configurations.

# Returns

- A `DataFrame` containing:
  - `cable_id`: Identifier of each cable design.
  - `horz`: Horizontal coordinate of each cable \\[m\\].
  - `vert`: Vertical coordinate of each cable \\[m\\].
  - `phase_mapping`: Human-readable string representation mapping each cable component to its assigned phase.

# Examples

```julia
df = $(FUNCTIONNAME)(cable_system)
println(df)
# Output:
# │ cable_id   │ horz │ vert  │ phase_mapping           │
# │------------│------│-------│-------------------------│
# │ "Cable1"   │ 0.0  │ -0.5  │ core: 1, sheath: 0      │
# │ "Cable2"   │ 0.35 │ -1.25 │ core: 2, sheath: 0      │
```

# See also

- [`LineCableSystem`](@ref)
- [`CablePosition`](@ref)
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
            ["$(name): $(phase)" for (name, phase) in zip(component_names, cable_position.conn)],
            ", ",
        )
        push!(mappings, mapping_str)
    end
    data = DataFrame(
        cable_id=cable_ids,
        horz=horz_coords,
        vert=vert_coords,
        phase_mapping=mappings
    )
    return data
end