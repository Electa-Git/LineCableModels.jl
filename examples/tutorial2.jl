#=
# Tutorial 2 - Building a cable design

This tutorial demonstrates how to model a typical medium-voltage single-core power cable
using the [`LineCableModels.jl`](@ref) package. The objective is to build a complete representation of a single-core 18/30 kV cable with a 1000 mm² aluminum conductor and 35 mm² copper screen.
=#

#=
**Tutorial outline**
```@contents
Pages = [
    "tutorial2.md",
]
Depth = 2:3
```
=#

#=
## Introduction

Single-core power cables have a complex structure consisting of multiple concentric layers, each with specific geometric and material properties -- for example, a cable of type NA2XS(FL)2Y 18/30 [is shown here](https://www.google.com/search?udm=2&q=%22NA2XS(FL)2Y%2018/30%20kV%20cable%22). Prior to building actual transmission line models that incorporate cables as part of the transmission system, e.g. for EMT simulations, power flow, harmonics, protection studies etc., it is necessary to determine the base (or DC) electrical parameters of the cable itself.

This tutorial covers:

1. Creating a detailed [`CableDesign`](@ref LineCableModels.DataModel.CableDesign) with all its components.
2. Examining the main electrical parameters (R, L, C) of the cable core [`ConductorGroup`](@ref LineCableModels.DataModel.ConductorGroup) and main [`InsulatorGroup`](@ref LineCableModels.DataModel.InsulatorGroup).
3. Examining the equivalent electromagnetic properties of every [`CableComponent`](@ref LineCableModels.DataModel.CableComponent) (core, sheath, jacket).
4. Saving the cable design to a [`CablesLibrary`](@ref) for future use.
4. Assigning [`CableDesign`](@ref LineCableModels.DataModel.CableDesign) objects to a [`LineCableSystem`](@ref LineCableModels.DataModel.LineCableSystem) and exporting the model to PSCAD for EMT analysis.
=#

#=
## Getting started
=#

# Load the public modeling API and the packages used for presentation:
using LineCableModels
import CairoMakie
using DataFrames
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_verbosity!(0); #hide
set_backend!(:cairo); #hide

# Initialize materials library with default values:
materials = MaterialsLibrary(add_defaults = true)
materials_df = DataFrame(materials)

#=
```julia
# Alternatively, it can be loaded from the example file built in the previous tutorial:
load!(materials, file_name = "materials_library.json")
```
=#

#=
## Cable dimensions

The cable under consideration is a medium-voltage, stranded aluminum conductor cable with XLPE insulation, copper wire concentric screens, water-blocking tape, and PE jacket that is rated for 18/30 kV systems. This information is typically found in the cable datasheet and is fully described in the code type under standards HD 620 10C [CENELEC_HD620_S3_2023](@cite) or DIN VDE 0276-620 [VDE_DIN_VDE_0276_620_2024](@cite):

```
NA2XS(FL)2Y
-----------
│ │   │  │
│ │   │  └── 2Y: Outer sheath of polyethylene (PE)
│ │   └── (FL): Longitudinal watertight protection
│ │
│ └── 2XS: XLPE insulation with screen of copper wires
└── NA: Aluminum conductor
```
=#

# After some research, it is found that a typical cable of this type has the following configuration:
num_co_wires = 61  # number of core wires
num_sc_wires = 49  # number of screen wires
d_core = 38.1e-3   # nominal core overall diameter
d_w = 4.7e-3       # nominal strand diameter of the core
t_sc_in = 0.6e-3   # nominal internal semicon thickness
t_ins = 8e-3       # nominal main insulation thickness
t_sc_out = 0.3e-3  # nominal external semicon thickness
d_ws = .95e-3      # nominal wire screen diameter
t_cut = 0.1e-3     # nominal thickness of the copper tape (around wire screens)
w_cut = 10e-3      # nominal width of copper tape
t_wbt = .3e-3      # nominal thickness of the water blocking tape
t_sct = .3e-3      # nominal thickness of the semiconductive tape
t_alt = .15e-3     # nominal thickness of the aluminum tape
t_pet = .05e-3     # nominal thickness of the pe face in the aluminum tape
t_jac = 2.4e-3     # nominal PE jacket thickness

d_overall = d_core #hide
layers = [] #hide
push!(layers, ("Conductor", missing, d_overall * 1000)) #hide
d_overall += 2 * t_sct #hide
push!(layers, ("Inner semiconductive tape", t_sct * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc_in #hide
push!(layers, ("Inner semiconductor", t_sc_in * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_ins #hide
push!(layers, ("Main insulation", t_ins * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc_out #hide
push!(layers, ("Outer semiconductor", t_sc_out * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sct #hide
push!(layers, ("Outer semiconductive tape", t_sct * 1000, d_overall * 1000)) #hide
d_overall += 2 * d_ws #hide
push!(layers, ("Wire screen", d_ws * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_cut #hide
push!(layers, ("Copper tape", t_cut * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_wbt #hide
push!(layers, ("Water-blocking tape", t_wbt * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_alt #hide
push!(layers, ("Aluminum tape", t_alt * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_pet #hide
push!(layers, ("PE with aluminum face", t_pet * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_jac #hide
push!(layers, ("PE jacket", t_jac * 1000, d_overall * 1000)); #hide

# The cable structure is summarized in a table for better visualization, with dimensions in milimiters:
df = DataFrame( #hide
    layer = first.(layers), #hide
    thickness = [ #hide
                 ismissing(t) ? "-" : round(t, sigdigits = 2) for t in getindex.(layers, 2) #hide
                 ], #hide
    diameter = [round(d, digits = 2) for d in getindex.(layers, 3)] #hide
) #hide

#=
## Describing the cable

!!! note "Materialized object hierarchy"
    [`CableBuilder`](@ref) is the user-facing description of the cable. When the
    description is resolved, the package builds the existing physical hierarchy
    shown below. That hierarchy validates every radius and incrementally computes
    the equivalent electrical properties used by the numerical methods.

```
CableDesign
├── CableComponent
│   ├── conductor_group::ConductorGroup <: AbstractConductorPart
│   │   ├── conductor_props::Material
│   │   └── layers::Vector{AbstractConductorPart}
│   │       ├── CircStrands
│   │       ├── Tubular
│   │       ├── Strip
│   │       └── …
│   └── insulator_group::InsulatorGroup <: AbstractInsulatorPart
│       ├── insulator_props::Material
│       └── layers::Vector{AbstractInsulatorPart}
│           ├── Insulator
│           ├── Semicon
│           └── …
⋮
├── CableComponent
│   ├── …
⋮   ⋮
```

### Cable designs

The materialized [`CableDesign`](@ref LineCableModels.DataModel.CableDesign) is the main container for all cable
components. Users normally obtain it by resolving a [`CableBuilder`](@ref)
description rather than constructing each intermediate object manually.

### Cable components

Each [`CableComponent`](@ref LineCableModels.DataModel.CableComponent) represents a functional group of the cable (core,
sheath, armor, or jacket), organized into conductive and insulating layers with
their respective effective material properties.

### Conductor groups

The conductive layers calculate equivalent resistance (R) and inductance (L)
while retaining the geometry of stranded wires, tubular conductors, and tapes.

#### AbstractConductorPart implementations

- `Conductor.Wires` describes stranded cores and screens with helical patterns
  and circular cross-sections.
- `Conductor.Tubular` describes a solid core or an annular conductor.
- `Conductor.Strip` describes a helically wound rectangular tape.

### Insulator groups

The insulating layers calculate the equivalent capacitance (C) and conductance
(G) of a concentric stack.

#### AbstractInsulatorPart implementations

- `Insulator.Tubular` describes a dielectric layer.
- `Insulator.Semicon` describes a semiconducting layer with intermediate
  resistivity and high permittivity.

!!! note "Equivalent circuit parameters"
    The hierarchy calculates equivalent circuit parameters by:

    1. Computing geometry-specific parameters for every conductive and insulating layer.
    2. Aggregating these into equivalent parameters for each cable component.
    3. Converting the composite structure into an equivalent coaxial model by matching lumped circuit quantities (R, L, C, G) to effective electromagnetic properties (ρ, ε, µ) at the [`CableComponent`](@ref LineCableModels.DataModel.CableComponent) level. The effective properties are stored in dedicated [`Material`](@ref) objects.
=#

#=
## Core and main insulation

The core consists of a central wire and four concentric AAAC layers with 61
wires arranged in a (1/6/12/18/24) pattern. The respective lay ratios are
(15/13.5/12.5/11) [CENELEC50182](@cite). `Conductor.Wires` retains the helical
pattern required for the resistance and GMR corrections.
=#

# Select reusable materials from the library:
aluminum = Material(materials, :aluminum)
copper = Material(materials, :copper)
polyacrylate = Material(materials, :polyacrylate)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe)

#=
!!! tip "Convenience methods"
    Cable datasheets usually report wire diameters, while the builder accepts
    `wire_radius`. The conversion `d_w / 2` therefore remains visible at the
    point of use. Subsequent layers are stacked over the current outer radius.
=#

# Describe the central wire and the four stranded layers:
core_wire_counts = (1, 6, 12, 18, 24)
@assert sum(core_wire_counts) == num_co_wires
core_conductors = (
    Conductor.Wires(
        :core; wire_radius=d_w / 2, num_wires=1, lay_ratio=0.0, material=aluminum),
    Conductor.Wires(
        :core; wire_radius=d_w / 2, num_wires=6, lay_ratio=15.0, material=aluminum),
    Conductor.Wires(
        :core; wire_radius=d_w / 2, num_wires=12, lay_ratio=13.5, material=aluminum),
    Conductor.Wires(
        :core; wire_radius=d_w / 2, num_wires=18, lay_ratio=12.5, material=aluminum),
    Conductor.Wires(
        :core; wire_radius=d_w / 2, num_wires=24, lay_ratio=11.0, material=aluminum),
)

#=
### Inner semiconductor

The inner semiconductor layer ensures uniform electric field distribution between
the conductor and insulation, eliminating air gaps and reducing field concentrations. An optional semiconductive tape is often used to ensure core uniformity and enhanced adherence.
=#

#=
!!! tip "Convenience methods"
    Annular layers accept either `radius=...` for an absolute outer radius or
    `thickness=...` for a radial increment. Exactly one is required. A thickness
    declaration automatically starts at the current outer radius.
=#

# Describe the inner semiconductive tape and semiconductor (1000 Ω·m as per IEC 840):
inner_layers = (
    Insulator.Semicon(:core; thickness=t_sct, material=polyacrylate),
    Insulator.Semicon(:core; thickness=t_sc_in, material=semicon1),
)

#=
### Main insulation

XLPE (cross-linked polyethylene) is the standard insulation material for modern
medium and high voltage cables due to its excellent dielectric properties.
=#

# Describe the main insulation layer:
main_insulation = Insulator.Tubular(:core; thickness=t_ins, material=pe)

#=
### Outer semiconductor

Similar to the inner semiconductor, the outer semiconductor provides a uniform
transition from insulation to the metallic screen.
=#

# Describe the outer semiconductor (500 Ω·m as per IEC 840) and tape:
outer_layers = (
    Insulator.Semicon(:core; thickness=t_sc_out, material=semicon2),
    Insulator.Semicon(:core; thickness=t_sct, material=polyacrylate),
)

# Assemble all declarations associated with the core component:
core_parts = (core_conductors, inner_layers, main_insulation, outer_layers)

#=
With the core parts properly defined, the [`CableDesign`](@ref LineCableModels.DataModel.CableDesign) object is initialized with nominal data from the datasheet. This includes voltage ratings and reference electrical parameters that will be used to benchmark the design.
=#

# Record the nominal values reported by the datasheet:
cable_id = "18kV_1000mm2"
datasheet_info = (
    designation_code = "NA2XS(FL)2Y",
    U0 = 18.0,                        # Phase-to-ground voltage [kV]
    U = 30.0,                         # Phase-to-phase voltage [kV]
    conductor_cross_section = 1000.0, # [mm²]
    screen_cross_section = 35.0,      # [mm²]
    resistance = 0.0291,              # DC resistance [Ω/km]
    capacitance = 0.39,               # Capacitance [μF/km]
    inductance = 0.3,                 # Inductance in trifoil [mH/km]
)

# Resolve the deterministic core description into a materialized cable design:
core_design = only(CableBuilder(cable_id, core_parts; nominal=datasheet_info))

# At this point, it becomes possible to preview the cable design:
plt1 = preview(
    core_design,
    display_plot=false, #hide
    controls=false, #hide
)
plt1.figure #hide

#=
### Wire screens

The metallic screen (typically copper) serves multiple purposes:
- Provides a return path for fault currents.
- Ensures radial symmetry of the electric field.
- Acts as electrical shielding.
- Provides mechanical protection.
=#

# Describe the wire screen on top of the previous component:
lay_ratio = 10.0 # typical value for wire screens
sheath_parts = (
    Conductor.Wires(
        :sheath;
        wire_radius=d_ws / 2,
        num_wires=num_sc_wires,
        lay_ratio,
        material=copper,
    ),
    Conductor.Strip(
        :sheath;
        thickness=t_cut,
        width=w_cut,
        lay_ratio,
        material=copper,
    ),
    Insulator.Semicon(:sheath; thickness=t_wbt, material=polyacrylate),
)

# Resolve and examine the core plus metallic screen:
screened_design = only(CableBuilder(
    cable_id, core_parts, sheath_parts; nominal=datasheet_info))

# Examine the newly added components:
plt2 = preview(
    screened_design,
    display_plot=false, #hide
    controls=false, #hide
)
plt2.figure #hide

#=
### Outer jacket components

Modern cables often include an aluminum tape as moisture barrier
and PE (polyethylene) outer jacket for mechanical protection.
=#

# Describe the aluminum moisture barrier, bonded PE face, and outer PE jacket:
jacket_parts = (
    Conductor.Tubular(:jacket; thickness=t_alt, material=aluminum),
    Insulator.Tubular(:jacket; thickness=t_pet, material=pe),
    Insulator.Tubular(:jacket; thickness=t_jac, material=pe),
)

#=
!!! tip "Convenience methods"
    A component name groups its conductive and insulating declarations. The
    builder resolves components in radial order and uses the previous
    component's outer radius as the next component's inner boundary.
=#

# Resolve the complete cable description:
cable_design = only(CableBuilder(
    cable_id,
    core_parts,
    sheath_parts,
    jacket_parts;
    nominal=datasheet_info,
))

# Inspect the finished cable design:
plt3 = preview(
    cable_design,
    display_plot=false, #hide
    controls=false, #hide
)
plt3.figure #hide

#=
## Examining the cable parameters (RLC)

In this section, the cable design is examined and the calculated parameters are compared with datasheet values. [`LineCableModels.jl`](@ref) provides methods to analyze the design in different levels of detail.
=#

# Calculate the cable constants explicitly. DataFrame presentation is applied
# only after the numerical result exists:
constants = compute!(CableConstantsProblem(cable_design), Formulation())

# Compare the calculated values with the datasheet information in the units
# conventionally used by cable manufacturers:
core_df = DataFrame(
    parameter = ["R", "L", "C"],
    calculated = [constants.R * 1e3, constants.L * 1e6, constants.C * 1e9],
    datasheet = [
        datasheet_info.resistance,
        datasheet_info.inductance,
        datasheet_info.capacitance,
    ],
    unit = ["Ω/km", "mH/km", "μF/km"],
)

# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)

# Get detailed description of all cable parts:
detailed_df = DataFrame(cable_design, :detailed)

#=
## Saving the cable design

!!! note "Cables library"
    Designs can be saved to a library for future use. The [`CablesLibrary`](@ref) is a container for storing multiple cable designs, allowing for easy access and reuse in different projects. Library management is performed using `DataFrame`, [`add!`](@ref), and [`save`](@ref).
=#

# Store the cable design and inspect the library contents:
library = CablesLibrary()
add!(library, cable_design)
library_df = DataFrame(library)

# Save to file for later use:
output_file = fullfile("cables_library.json")
save(library, file_name = output_file);

# Load the saved design into a fresh library and retrieve it by identifier:
loaded_library = CablesLibrary()
load!(loaded_library, file_name=output_file)
loaded_design = get(loaded_library, cable_id)
loaded_library_df = DataFrame(loaded_library)

#=
### Defining a cable system

!!! note "Cable systems"
    [`SystemBuilder`](@ref) combines a cable design, its positions, length,
    operating temperature, earth properties, and analysis frequencies. Resolving
    this description produces the complete line-parameter problem used by
    [`compute!`](@ref), preview, and export routines.
=#

#=
### Earth model

The earth return path significantly affects cable impedance calculations and needs to be properly modeled. In this tutorial, only a basic model with typical soil properties is defined. This will be further elaborated in the subsequent tutorials.
=#

# Define a frequency scan and typical homogeneous-soil properties:
f = collect(10.0 .^ range(0, stop=6, length=10)) # 1 Hz to 1 MHz
earth = Earth(rho=100.0, eps_r=10.0, mu_r=1.0)

#=
### Three-phase system in trifoil configuration

This section ilustrates the construction of a cable system with three identical cables arranged in a trifoil formation.
=#

# Describe three cables touching in trifoil at 1 m burial depth. The spacing is
# the center-to-center distance:
formation = trifoil(
    x=0.0,
    y=-1.0,
    spacing=70e-3,
    phases=(
        :core => (1, 2, 3),
        :sheath => 0,
        :jacket => 0,
    ),
)

# Combine the loaded design, formation, earth, and frequency scan:
problem = only(SystemBuilder(
    "18kV_1000mm2_trifoil",
    loaded_design,
    formation;
    length=1000.0,
    temperature=20.0,
    earth,
    frequencies=f,
))
cable_system = problem.system
earth_params = problem.earth_props

# Earth model base properties:
earthmodel_df = DataFrame(earth_params)

#=
!!! note "Phase mapping"
    The `phases` declaration maps each cable component to its electrical phase.
    The core tuple `(1, 2, 3)` assigns phases A, B, and C to the three cables.
    Components mapped to phase 0 are grounded and subsequently Kron-eliminated.
    Components assigned to the same positive phase are bundled.
=#

#=
### Cable system preview

In this section the complete three-phase cable system is examined.
=#

# Display system details:
system_df = DataFrame(cable_system)

# Visualize the cross-section of the three-phase system:
plt4 = preview(
    cable_system,
    earth_model=earth_params,
    zoom_factor=2.0,
    display_plot=false, #hide
    controls=false, #hide
)
plt4.figure #hide

#=
## PSCAD & ATPDraw export

The final step showcases how to export the model for electromagnetic transient simulations in EMT-type software.
=#

# Export to PSCAD input file:
output_file = fullfile("pscad_export.pscx")
export_file = export_data(:pscad, cable_system, earth_params, file_name = output_file);

# Export to ATPDraw project file (XML):
output_file = fullfile("atp_export.xml")
export_file = export_data(:atp, cable_system, earth_params, file_name = output_file);

#=
## Conclusion

This tutorial has demonstrated how to:

1. Describe and preview a complex power cable with multiple concentric layers.
2. Calculate and compare its base parameters (R, L, C) with datasheet values.
3. Save the design, load it in a fresh library, and reuse it.
4. Build and preview a three-phase cable system in trifoil arrangement.
5. Export the physical model for PSCAD and ATPDraw.

[`LineCableModels.jl`](@ref) provides detailed routines for power cable modeling
with a physically meaningful representation of all cable components. This approach
ensures that electromagnetic parameters are calculated with high precision. Now you can go ahead and run these cable simulations like a boss!
=#
