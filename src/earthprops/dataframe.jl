"""
$(TYPEDSIGNATURES)

Generates a `DataFrame` summarizing basic properties of earth layers from an [`EarthModel`](@ref).

# Arguments

- `earth_model`: Instance of [`EarthModel`](@ref) containing earth layers.

# Returns

- A `DataFrame` with columns:
  - `rho_g`: Base (DC) resistivity of each layer \\[Ω·m\\].
  - `epsr_g`: Base (DC) relative permittivity of each layer \\[dimensionless\\].
  - `mur_g`: Base (DC) relative permeability of each layer \\[dimensionless\\].
  - `thickness`: Thickness of each layer \\[m\\].

# Examples

```julia
df = $(FUNCTIONNAME)(earth_model)
println(df)
```
"""
function DataFrames.DataFrame(earth_model::EarthModel)
  layers = earth_model.layers

  base_rho_g = [layer.base_rho_g for layer in layers]
  base_epsr_g = [layer.base_epsr_g for layer in layers]
  base_mur_g = [layer.base_mur_g for layer in layers]
  thickness = [layer.t for layer in layers]

  return DataFrame(
    rho_g=base_rho_g,
    epsr_g=base_epsr_g,
    mur_g=base_mur_g,
    thickness=thickness,
  )
end