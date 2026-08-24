import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Convert a material library to one row per stored material.

# Arguments

- `library`: Material library.

# Returns

- A `DataFrame` with `name`, `rho`, `eps_r`, `mu_r`, `T0`, and `alpha`
  columns.

"""
function DataFrame(library::MaterialsLibrary)::DataFrame
    rows = [(
                name = name,
                rho = m.rho,
                eps_r = m.eps_r,
                mu_r = m.mu_r,
                T0 = m.T0,
                alpha = m.alpha
            )
            for (name, m) in library]
    data = DataFrame(rows)
    return data
end
