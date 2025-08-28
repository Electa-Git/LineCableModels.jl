#=
# Tutorial 1 - Using the  materials library

This tutorial demonstrates how to manage material properties for power cable modeling using the package [`LineCableModels.jl`](@ref). Accurate knowledge of electromagnetic properties is essential for reliable cable design and analysis.

Beyond showcasing the API, this guide serves as a practical reference by providing standard property values from recognized industry sources like CIGRE TB-531 [cigre531](@cite) and IEC 60287 [IEC60287](@cite) that can be stored and consistently applied across multiple design iterations and simulation studies.
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
##   Getting started
=#

# Load the package:
using DataFrames
using LineCableModels
fullfile(filename) = joinpath(@__DIR__, filename); #hide
setup_logging!(0); #hide

#=
The [`MaterialsLibrary`](@ref) is a container for storing electromagnetic properties of 
different materials used in power cables. By default, it initializes with several common 
materials with their standard properties.
=#

# Initialize a [`MaterialsLibrary`](@ref) with default values:
materials = MaterialsLibrary()

# Inspect the contents of the materials library:
materials_df = DataFrame(materials)

#=
The function [`DataFrame`](@ref) returns a `DataFrame` with all materials and their properties, namely: electrical resistivity, relative permittivity, relative permeability, reference temperature, and temperature coefficient.
=#

# ##   Adding new materials
#=
!!! note "Note"
	New materials can be added to the library using the [`Material`](@ref) constructor followed by [`add!`](@ref).

It might be useful to add other conductor materials with corrected properties based on recognized standards [cigre531](@cite) [IEC60287](@cite).
=#

copper_corrected = Material(1.835e-8, 1.0, 0.999994, 20.0, 0.00393) # Copper with corrected resistivity from IEC 60287-3-2
add!(materials, "copper_corrected", copper_corrected)
aluminum_corrected = Material(3.03e-8, 1.0, 0.999994, 20.0, 0.00403) # Aluminum with corrected resistivity from IEC 60287-3-2
add!(materials, "aluminum_corrected", aluminum_corrected)
lead = Material(21.4e-8, 1.0, 0.999983, 20.0, 0.00400) # Lead or lead alloy
add!(materials, "lead", lead)
steel = Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450) # Steel
add!(materials, "steel", steel)
bronze = Material(3.5e-8, 1.0, 1.0, 20.0, 0.00300) # Bronze
add!(materials, "bronze", bronze)
stainless_steel = Material(70.0e-8, 1.0, 500.0, 20.0, 0.0) # Stainless steel
add!(materials, "stainless_steel", stainless_steel)

#=
When modeling cables for EMT analysis, one might be concerned with the impact of insulators and semiconductive layers on cable constants. Common insulation materials and semicons with different dielectric properties are reported in Table 6 of [cigre531](@cite). Let us include some of these materials in the [`MaterialsLibrary`](@ref) to help our future selves.
=#

epr = Material(1e15, 3.0, 1.0, 20.0, 0.005) # EPR (ethylene propylene rubber)
add!(materials, "epr", epr)
pvc = Material(1e15, 8.0, 1.0, 20.0, 0.1) # PVC (polyvinyl chloride)
add!(materials, "pvc", pvc)
laminated_paper = Material(1e15, 2.8, 1.0, 20.0, 0.0) # Laminated paper propylene
add!(materials, "laminated_paper", laminated_paper)
carbon_pe = Material(0.06, 1e3, 1.0, 20.0, 0.0) # Carbon-polyethylene compound (semicon)
add!(materials, "carbon_pe", carbon_pe)
conductive_paper = Material(18.5, 8.6, 1.0, 20.0, 0.0) # Conductive paper layer (semicon)
add!(materials, "conductive_paper", conductive_paper)

# ##  Removing materials
#=
!!! note "Note"
	Materials can be removed from the library with the [`delete!`](@ref) function.
=#

# Add a duplicate material by accident:
add!(materials, "epr_dupe", epr)

# And now remove it using the [`delete!`](@ref) function:
delete!(materials, "epr_dupe")

# Examine the updated library after removing the duplicate:
println("Material properties compiled from CIGRE TB-531 and IEC 60287:")
materials_df = DataFrame(materials)

# ##  Saving the materials library to JSON
output_file = fullfile("materials_library.json")
save(materials, file_name=output_file);


# ##  Retrieving materials for use
#=
!!! note "Note"
	To load from an existing JSON file, instantiate a new [`MaterialsLibrary`](@ref) followed by a call to the [`load!`](@ref) method. Materials can be retrieved from the library using the [`get`](@ref) function.
=#

# Initialize a new [`MaterialsLibrary`](@ref) and load from the JSON file:
materials_from_json = MaterialsLibrary()
load!(materials_from_json, file_name=output_file)

# Retrieve a material and display the object:
copper = get(materials_from_json, "copper_corrected")

# Access the material properties:
println("Retrieved copper_corrected material properties:")
println("Resistivity: $(copper.rho) Ω·m")
println("Relative permittivity: $(copper.eps_r)")
println("Relative permeability: $(copper.mu_r)")
println("Reference temperature: $(copper.T0) °C")
println("Temperature coefficient: $(copper.alpha) 1/°C")

# ##  Conclusion
#=
This tutorial has demonstrated how to:

1. Initialize a [`MaterialsLibrary`](@ref) with default [`Material`](@ref) objects.
2. Add new materials with specific properties.
3. Remove duplicate materials.
4. Save the library to a file for future use.
5. Retrieve materials for use in cable modeling.

The [`MaterialsLibrary`](@ref) provides a flexible and traceable framework to manage material properties for accurate power cable modeling. Custom [`Material`](@ref) objects can be defined and used to match specific manufacturer data or standards requirements.
=#
