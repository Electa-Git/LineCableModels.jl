import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Lists the cable designs in a [`CablesLibrary`](@ref) object as a `DataFrame`.

# Arguments

- `library`: An instance of [`CablesLibrary`](@ref) whose cable designs are to be displayed.

# Returns

- A `DataFrame` object with the following columns:
  - `cable_id`: The unique identifier for each cable design.
  - `nominal_data`: A string representation of the nominal data for each cable design.
  - `components`: A comma-separated string listing the components of each cable design.

# Examples

```julia
library = CablesLibrary()
design1 = CableDesign("example1", nominal_data=NominalData(...), components=Dict("A"=>...))
design2 = CableDesign("example2", nominal_data=NominalData(...), components=Dict("C"=>...))
add!(library, design1)
add!(library, design2)

# Display the library as a DataFrame
df = $(FUNCTIONNAME)(library)
first(df, 5)  # Show the first 5 rows of the DataFrame
```

# See also

- [`CablesLibrary`](@ref)
- [`CableDesign`](@ref)
- [`add!`](@ref)
"""
function DataFrame(library::CablesLibrary)::DataFrame
    ids = keys(library)
    nominal_data = [string(design.nominal_data) for design in values(library)]
    components = [
        join([comp.id for comp in design.components], ", ") for
        design in values(library)
    ]
    df = DataFrame(
        cable_id=collect(ids),
        nominal_data=nominal_data,
        components=components,
    )
    return (df)
end