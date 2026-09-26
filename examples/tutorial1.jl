#=
# Tutorial 1 - Using the materials library

Create, inspect, modify, save, and reload a materials library for cable
calculations. The examples use electromagnetic property values from CIGRE
TB-531 [cigre531](@cite) and IEC 60287 [IEC60287](@cite).
=#

#=
**Tutorial outline**
```@contents
Pages = [
    "tutorial1.md",
]
Depth = 2:3
```
=#

#=
## Getting started

[`MaterialsLibrary`](@ref) stores [`Material`](@ref) values under string names.
Its default constructor includes the package's built-in material records.
The library supports ordinary dictionary operations for lookup, replacement,
and removal.
=#

using LineCableModels: Material, MaterialsLibrary, add!, load!, save

# Create and inspect a library containing the built-in materials:
materials = MaterialsLibrary()

# The library display summarizes the stored names and material classes.
# Retrieve an individual material to inspect its properties:
materials["copper"]

#=
## Material properties

The material declarations below specify the following properties:

| Keyword | Meaning | Unit |
|:--|:--|:--|
| `kind` | Material class: `:conductor`, `:insulator`, or `:semicon` in this tutorial | None |
| `rho` | Electrical resistivity at `T0` | Ω·m |
| `eps_r` | Relative permittivity | Dimensionless |
| `mu_r` | Relative permeability | Dimensionless |
| `T0` | Reference temperature | °C |
| `alpha` | Temperature coefficient of resistivity | 1/°C |

A material's library name identifies the entry; `kind` describes its physical
class. For example, several entries can have `kind=:conductor` while retaining
different resistivities and names.

Construct a material with named properties, then use [`add!`](@ref) to insert
it under a new name. `add!` rejects an existing name rather than silently
replacing its value.
=#

#=
## Adding conductor materials

The following records retain the conductor values used in this example from
[cigre531](@cite) and [IEC60287](@cite). The corrected copper and aluminum
resistivities are attributed to IEC 60287-3-2.

The built-in copper, aluminum, lead, and steel entries remain unchanged.
Distinct names keep these additional records available alongside them.
=#

copper_corrected = Material(
    kind = :conductor,
    rho = 1.835e-8,
    eps_r = 1.0,
    mu_r = 0.999994,
    T0 = 20.0,
    alpha = 0.00393,
)
add!(materials, "copper_corrected", copper_corrected);

aluminum_corrected = Material(
    kind = :conductor,
    rho = 3.03e-8,
    eps_r = 1.0,
    mu_r = 0.999994,
    T0 = 20.0,
    alpha = 0.00403,
)
add!(materials, "aluminum_corrected", aluminum_corrected);

# Lead or lead alloy:
lead = Material(
    kind = :conductor,
    rho = 21.4e-8,
    eps_r = 1.0,
    mu_r = 0.999983,
    T0 = 20.0,
    alpha = 0.00400,
)
add!(materials, "lead_reference", lead);

steel = Material(
    kind = :conductor,
    rho = 13.8e-8,
    eps_r = 1.0,
    mu_r = 300.0,
    T0 = 20.0,
    alpha = 0.00450,
)
add!(materials, "steel_reference", steel);

bronze = Material(
    kind = :conductor,
    rho = 3.5e-8,
    eps_r = 1.0,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.00300,
)
add!(materials, "bronze", bronze);

stainless_steel = Material(
    kind = :conductor,
    rho = 70.0e-8,
    eps_r = 1.0,
    mu_r = 500.0,
    T0 = 20.0,
    alpha = 0.0,
)
add!(materials, "stainless_steel", stainless_steel);

#=
## Adding insulation and semiconducting materials

The following declarations retain the insulation and semiconducting material
values used in this example, with dielectric properties attributed to Table 6
of [cigre531](@cite).
=#

# EPR: ethylene propylene rubber.
epr = Material(
    kind = :insulator,
    rho = 1e15,
    eps_r = 3.0,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.005,
)
add!(materials, "epr", epr);

# PVC: polyvinyl chloride.
pvc = Material(
    kind = :insulator,
    rho = 1e15,
    eps_r = 8.0,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.1,
)
add!(materials, "pvc", pvc);

# Laminated paper propylene.
laminated_paper = Material(
    kind = :insulator,
    rho = 1e15,
    eps_r = 2.8,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.0,
)
add!(materials, "laminated_paper", laminated_paper);

# Carbon-polyethylene semiconducting compound.
carbon_pe = Material(
    kind = :semicon,
    rho = 0.06,
    eps_r = 1e3,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.0,
)
add!(materials, "carbon_pe", carbon_pe);

# Conductive paper used as a semiconducting layer.
conductive_paper = Material(
    kind = :semicon,
    rho = 18.5,
    eps_r = 8.6,
    mu_r = 1.0,
    T0 = 20.0,
    alpha = 0.0,
)
add!(materials, "conductive_paper", conductive_paper);

# Inspect the populated library:
materials

#=
## Modifying a material entry

`Material` values are immutable. To change properties, construct a new value
from an existing material and override the desired keywords. Unspecified
properties are retained.

Indexed assignment explicitly replaces a library entry. This differs from
`add!`, which requires a new name. Use a temporary copper entry to demonstrate
the operation without altering any of the reference records above.
=#

add!(materials, "copper_trial", materials["copper"]);

materials["copper_trial"] = Material(
    materials["copper_trial"];
    rho = 1.835e-8,
);

# The trial has the new resistivity; the built-in copper is unchanged:
materials["copper_trial"]
materials["copper"]

# Remove the temporary entry before saving the reference library:
delete!(materials, "copper_trial");

#=
## Removing materials

Use `delete!` with the stored name. Removing an entry does not delete a material
kept under another name or change a value already retrieved from the library.
The following example adds and removes an extra EPR entry while retaining
`"epr"`.
=#

add!(materials, "epr_dupe", epr);
delete!(materials, "epr_dupe");

# Inspect the library after removing the duplicate:
materials

#=
## Saving the materials library to JSON

Save the current library to `materials_library.json` beside this tutorial.
The file includes both the built-in records and the additional entries.
Saving to an existing filename replaces that file's contents.
=#

output_file = joinpath(@__DIR__, "materials_library.json")
save(materials; file_name = output_file);

#=
## Loading and retrieving materials

Start with an empty library when the saved file supplies all its entries.
`load!` replaces the destination library's contents with the file's records;
it does not merge them with built-in defaults. A parsing or validation failure
leaves the destination unchanged.
=#

materials_from_json = MaterialsLibrary(add_defaults = false)
load!(materials_from_json; file_name = output_file);
materials_from_json

# Retrieve the corrected copper by its stored name. The material's native
# display shows its electromagnetic properties and reference temperature:
copper = materials_from_json["copper_corrected"]

#=
### Reading individual properties

Named material properties are available directly. Their units are those listed
in the material-property table above; reading a property does not perform a
calculation or modify the library.
=#

# Electrical resistivity at the reference temperature, in Ω·m:
copper.rho

# Relative permittivity and permeability:
copper.eps_r
copper.mu_r

# Reference temperature in °C and resistivity coefficient in 1/°C:
copper.T0
copper.alpha

#=
### Optional lookup

Indexed lookup raises `KeyError` when a required entry is absent. Use `get`
with an explicit fallback for an optional lookup; it does not insert a record.
=#

get(materials_from_json, "unlisted_material", nothing)
