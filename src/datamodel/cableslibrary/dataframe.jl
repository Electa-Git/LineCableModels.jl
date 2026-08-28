import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a cable-design library to one row per stored design.

# Arguments

- `library`: Cable-design library.

# Returns

- A `DataFrame` with `cable_id`, `catalogue`, and a comma-separated
  `terminals` column.

"""
function DataFrame(library::CablesLibrary)::DataFrame
    ids = keys(library)
    catalogue_data = [string(catalogue(library, id)) for id in ids]
    terminals = [join(string.(library[id].terminal_order), ", ") for id in ids]
    df = DataFrame(
        cable_id = collect(ids),
        catalogue = catalogue_data,
        terminals = terminals
    )
    return (df)
end
