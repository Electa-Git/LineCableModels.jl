import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Lists the contents of a [`MaterialsLibrary`](@ref) as a `DataFrame`.

# Arguments

- `library`: Instance of [`MaterialsLibrary`](@ref) to be displayed.

# Returns

- A `DataFrame` containing the material properties.

# Examples

```julia
library = MaterialsLibrary()
df = $(FUNCTIONNAME)(library)
```

# See also

- [`LineCableModels.ImportExport.save`](@ref)
"""
function DataFrame(library::MaterialsLibrary)::DataFrame
    rows = [
        (
            name=name,
            rho=m.rho,
            eps_r=m.eps_r,
            mu_r=m.mu_r,
            T0=m.T0,
            alpha=m.alpha,
        )
        for (name, m) in library
    ]
    data = DataFrame(rows)
    return data
end