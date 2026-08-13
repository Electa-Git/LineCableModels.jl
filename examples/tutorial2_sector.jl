#=
# Tutorial 2a - Building a sector-shaped cable design

This tutorial demonstrates how to model a typical low-voltage three-core power cable with sector-shaped conductors
using the [`LineCableModels.jl`](@ref) package. The objective is to build a complete representation of a three-core 1 kV cable with 95 mm² aluminum sector-shaped conductors and a concentric copper neutral.
=#

#=
**Tutorial outline**
```@contents
Pages = [
	"tutorial2_sector.md",
]
Depth = 2:3
```
=#

#=
## Introduction

Three-core power cables with sector-shaped conductors are common in low-voltage distribution networks. Their compact design allows for efficient use of space. This tutorial will guide you through creating a detailed [`CableDesign`](@ref) for such a cable.

This tutorial covers:

1. Defining materials with corrected resistivity.
2. Creating sector-shaped conductors using [`SectorParams`](@ref) and [`Sector`](@ref).
3. Assembling a multi-core [`CableDesign`](@ref).
4. Modeling a concentric neutral wire array.
5. Previewing the final cable design.
=#

#=
## Getting started
=#

# Load the package and set up the environment:
using LineCableModels
using DataFrames
import LineCableModels.BackendHandler: renderfig #hide
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_verbosity!(0); #hide

#=
## Cable and Material Data

We start by defining the materials. We will create a custom aluminum material with a resistivity corrected based on its nominal DC resistance, and a PVC material for insulation.
=#

# Initialize materials library and add a PVC material
materials = MaterialsLibrary(add_defaults=true)
pvc = Material(Inf, 8.0, 1.0, 20.0, 0.1) # simple PVC
add!(materials, "pvc", pvc)
copper = get(materials, "copper")
aluminum = get(materials, "aluminum")


#=
## Sector-Shaped Core Conductors

The core of the cable consists of three sector-shaped aluminum conductors. We define the geometry using `SectorParams` based on datasheet or standard values.
=#

# === Sector (core) geometry (table data) ===
# Based on Urquhart's paper for a 3-core 95mm^2 cable
n_sectors = 3
r_back_mm = 10.24       # sector radius b
d_sector_mm = 9.14        # sector depth s
r_corner_mm = 1.02      # corner radius c
theta_cond_deg = 119.0    # sector angle φ
ins_thick = 1.1e-3        # core insulation thickness

sector_params = SectorParams(
    n_sectors,
    r_back_mm / 1000,
    d_sector_mm / 1000,
    r_corner_mm / 1000,
    theta_cond_deg,
    ins_thick
)

#=
With the sector parameters defined, we can create the individual `Sector` conductors and their insulation. Each sector is rotated to form the 3-core bundle.
=#

rot_angles = (0.0, 120.0, 240.0)
sectors = [Sector(sector_params, ang, aluminum) for ang in rot_angles]
insulators = [SectorInsulator(sectors[i], ins_thick, pvc) for i in 1:3]

components = [
    CableComponent("core1", ConductorGroup(sectors[1]), InsulatorGroup(insulators[1])),
    CableComponent("core2", ConductorGroup(sectors[2]), InsulatorGroup(insulators[2])),
    CableComponent("core3", ConductorGroup(sectors[3]), InsulatorGroup(insulators[3]))
]

#=
## Concentric Neutral and Outer Jacket

The cable includes a concentric neutral conductor made of copper wires and an outer PVC jacket.
=#

# === Concentric neutral (30 wires) ===
n_neutral = 30
r_strand = 0.79e-3
R_N = 14.36e-3 # radius to center of neutral wires
R_O = 17.25e-3 # outer radius of the cable

inner_radius_neutral = R_N - r_strand
outer_jacket_thickness = R_O - (R_N + r_strand)

neutral_wires = WireArray(
    inner_radius_neutral,
    Diameter(2*r_strand),
    n_neutral,
    0.0, # lay ratio
    copper
)

neutral_jacket = Insulator(neutral_wires, Thickness(outer_jacket_thickness), pvc)
neutral_component = CableComponent("neutral", ConductorGroup(neutral_wires), InsulatorGroup(neutral_jacket))

#=
## Assembling the Cable Design

Now we assemble the complete `CableDesign` by adding all the components.
=#

design = CableDesign("NAYCWY_O_3x95_30x2_5", components[1])
add!(design, components[2])
add!(design, components[3])
add!(design, neutral_component)

#=
## Examining the Cable Design

We can now display a summary of the cable design and preview it graphically.
=#

println("Cable design summary:")
detailed_df = DataFrame(design, :detailed)
display(detailed_df)

println("Previewing cable design...")
plt, _ = preview(design)
plt #hide

#= 
## Storing in a Library

Finally, we can store the cable design in a `CablesLibrary` for future reference.
=#

library = CablesLibrary()
add!(library, design)
library_df = DataFrame(library)

# Save to file for later use:
output_file = fullfile("cables_library.json")
save(library, file_name = output_file);

#=
## Conclusion

This tutorial has demonstrated how to model a three-core cable with sector-shaped conductors. Key takeaways include:

1.  Creating custom materials with corrected properties.
2.  Defining complex conductor shapes like sectors.
3.  Assembling a multi-core cable design component by component.
4.  Visualizing the final design for verification.

This detailed modeling capability allows for accurate analysis of various cable configurations.
=#

