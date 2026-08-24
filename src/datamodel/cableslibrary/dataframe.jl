import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a cable-design library to one row per stored design.

# Arguments

- `library`: Cable-design library.

# Returns

- A `DataFrame` with `cable_id`, `nominal_data`, and a comma-separated
  `components` column.

"""
function DataFrame(library::CablesLibrary)::DataFrame
    ids = keys(library)
    nominal_data = [string(design.nominal_data) for design in values(library)]
    components = [join([comp.id for comp in design.components], ", ")
                  for
                  design in values(library)]
    df = DataFrame(
        cable_id = collect(ids),
        nominal_data = nominal_data,
        components = components
    )
    return (df)
end
