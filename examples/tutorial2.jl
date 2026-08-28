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

1. Building a [`CableDesign`](@ref) from physical regions and placement rules.
2. Examining the resolved geometry and the retained electrical terminals.
3. Calculating the cable's base resistance, inductance, and capacitance.
4. Saving the design to a [`CablesLibrary`](@ref) for future use.
5. Placing designs in a [`LineCableSystem`](@ref) and exporting the system for EMT analysis.
=#

#=
## Getting started
=#

# Load the public modeling API and the packages used for presentation:
using LineCableModels
import CairoMakie
using DataFrames
fullfile(filename) = joinpath(@__DIR__, filename); #hide
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

layer_names = ( #hide
    "Conductor", "Inner semiconductive tape", "Inner semiconductor", #hide
    "Main insulation", "Outer semiconductor", "Outer semiconductive tape", #hide
    "Wire screen", "Copper tape", "Water-blocking tape", "Aluminum tape", #hide
    "PE with aluminum face", "PE jacket" #hide
) #hide
layer_thicknesses = ( #hide
    missing, t_sct, t_sc_in, t_ins, t_sc_out, t_sct, d_ws, t_cut, t_wbt, #hide
    t_alt, t_pet, t_jac #hide
) #hide
radial_increments = ( #hide
    0.0, t_sct, t_sc_in, t_ins, t_sc_out, t_sct, d_ws, t_cut, t_wbt, #hide
    t_alt, t_pet, t_jac #hide
) #hide
layer_diameters = d_core .+ 2 .* cumsum(radial_increments) #hide

# The cable structure is summarized in a table for better visualization, with dimensions in milimiters:
df = DataFrame( #hide
    layer = layer_names, #hide
    thickness = [ #hide
                 ismissing(t) ? "-" : round(1000t, sigdigits = 2) for t in layer_thicknesses #hide
                 ], #hide
    diameter = round.(1000 .* layer_diameters, digits = 2) #hide
) #hide

#=
## Describing the cable

The design is one outward physical declaration. `Conductor`, `Insulator`, and
`Semiconductor` provide the recurring cable constructions; [`Stack`](@ref)
states their order; and [`build`](@ref) completes the design once. No partial
design is needed between layers.
=#

#=
## Core and main insulation

The core consists of a central wire and four concentric AAAC layers with 61
wires arranged in a (1/6/12/18/24) pattern. The respective lay ratios are
(15/13.5/12.5/11) [CENELEC50182](@cite). `Conductor.Stranded` retains every
physical strand and the layer-specific longitudinal paths required for
resistance and GMR calculations.
=#

# Select reusable materials from the library:
aluminum = Material(materials, :aluminum)
copper = Material(materials, :copper)
polyacrylate = Material(materials, :polyacrylate)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe)

# State the actual strand count and lay of every noncentral layer. These are the
# physical declarations reported for this conductor, not inputs from which the
# API must infer another layout:
stranded_core = Conductor.Stranded(
    :core,
    DiskDefinition(d_w / 2),
    aluminum;
    counts = (1, 6, 12, 18, 24),
    lay = (15.0, 13.5, 12.5, 11.0)
)

#=
### Inner semiconductor

The inner semiconductor layer ensures uniform electric field distribution between
the conductor and insulation, eliminating air gaps and reducing field concentrations. An optional semiconductive tape is often used to ensure core uniformity and enhanced adherence.
=#

#=
!!! tip "Contextual shells"
    [`ShellDefinition`](@ref) stores a radial thickness rather than a duplicated inner
    radius. [`Stack`](@ref) resolves it against the preceding outward boundary.
=#

#=
### Main insulation

XLPE (cross-linked polyethylene) is the standard insulation material for modern
medium and high voltage cables due to its excellent dielectric properties.
=#

#=
### Outer semiconductor

Similar to the inner semiconductor, the outer semiconductor provides a uniform
transition from insulation to the metallic screen.
=#

#=
### Wire screens

The metallic screen (typically copper) serves multiple purposes:
- Provides a return path for fault currents.
- Ensures radial symmetry of the electric field.
- Acts as electrical shielding.
- Provides mechanical protection.
=#

#=
### Outer jacket regions

Modern cables often include an aluminum tape as moisture barrier
and PE (polyethylene) outer jacket for mechanical protection.
=#

# The wire screen, copper tape, and aluminum barrier are absolute geometries.
# Their coordinates are stated once; contextual dielectric shells still need
# only their physical thickness.
conductor_outer = 9d_w / 2
screen_wire_locus = conductor_outer + t_sct + t_sc_in + t_ins +
                    t_sc_out + t_sct + d_ws / 2
screen_tape_inner = screen_wire_locus + d_ws / 2
screen_tape_outer = screen_tape_inner + t_cut
screen_tape_span = w_cut / ((screen_tape_inner + screen_tape_outer) / 2)
barrier_inner = screen_tape_outer + t_wbt
barrier_outer = barrier_inner + t_alt

# Declare the complete cable in physical order. Every item is an ordinary v1
# physical object and the declaration is completed by one `build` call.
cable_root = Stack(
    stranded_core,
    Semiconductor.Shell(:core_semicon_tape_inner, polyacrylate; t = t_sct),
    Semiconductor.Shell(:core_semicon_inner, semicon1; t = t_sc_in),
    Insulator.Shell(:core_insulation, pe; t = t_ins),
    Semiconductor.Shell(:core_semicon_outer, semicon2; t = t_sc_out),
    Semiconductor.Shell(:core_semicon_tape_outer, polyacrylate; t = t_sct),
    Conductor.Wires(
        :sheath,
        DiskDefinition(d_ws / 2),
        copper;
        n = num_sc_wires,
        r = screen_wire_locus,
        lay = 10.0
    ),
    Conductor.Strip(
        :sheath,
        SectorDefinition(
            screen_tape_inner,
            screen_tape_outer,
            -screen_tape_span / 2,
            screen_tape_span
        ),
        copper;
        lay = 10.0
    ),
    Semiconductor.Shell(:sheath_water_blocking, polyacrylate; t = t_wbt),
    Conductor.Tubular(:jacket, AnnulusDefinition(barrier_inner, barrier_outer), aluminum),
    Insulator.Shell(:jacket_pe_face, pe; t = t_pet),
    Insulator.Shell(:jacket_insulation, pe; t = t_jac)
)

# Keep catalogue data beside the physical model rather than inside it:
cable_id = "18kV_1000mm2"
datasheet_info = (
    designation_code = "NA2XS(FL)2Y",
    U0 = 18.0,                        # Phase-to-ground voltage [kV]
    U = 30.0,                         # Phase-to-phase voltage [kV]
    conductor_cross_section = 1000.0, # [mm²]
    screen_cross_section = 35.0,      # [mm²]
    resistance = 0.0291,              # DC resistance [Ω/km]
    capacitance = 0.39,               # Capacitance [μF/km]
    inductance = 0.3                  # Inductance in trefoil [mH/km]
)
cable_design = build(
    CableDesign,
    cable_id,
    cable_root
)
cable_library = CablesLibrary()
add!(cable_library, cable_design; catalogue = datasheet_info)

# Inspect the one completed physical design:
cable_plot = preview(
    cable_design,
    display_plot = false, #hide
    controls = false #hide
)
cable_plot.figure #hide

#=
## Examining the cable parameters (RLC)

In this section, the cable design is examined and the calculated parameters are compared with datasheet values. [`LineCableModels.jl`](@ref) provides methods to analyze the design in different levels of detail.
=#

# Calculate the cable constants explicitly. Scientific extraction and tabular
# presentation are separate consumers of the completed result:
constants = CableConstants(cable_design)

# The native DataFrames adapter returns the complete R/L/C result in its stored
# per-length units:
constants_table = DataFrame(constants)

# Materialize the homogeneous equivalent only when that design is
# itself the requested product:
equivalent_design = equivalent(cable_design; new_id = cable_id * "_equivalent")
equivalent_regions = DataFrame(equivalent_design)

# `observables` publishes detached values in the units conventionally used by
# cable manufacturers. The two rows identify real comparison sources; the
# physical quantities remain separate columns:
published_constants = observables(constants, (R, L, C))
datasheet_comparison = DataFrame(
    source = ("calculated", "datasheet"),
    R = (published_constants[1].values, datasheet_info.resistance),
    L = (published_constants[2].values, datasheet_info.inductance),
    C = (published_constants[3].values, datasheet_info.capacitance)
)
comparison_units = map(payload -> payload.unit, published_constants)

# Inspect every resolved physical region:
regions_df = DataFrame(cable_design)

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
load!(loaded_library, file_name = output_file)
loaded_design = get(loaded_library, cable_id)
loaded_library_df = DataFrame(loaded_library)

#=
### Defining a cable system

!!! note "Cable systems"
    [`LineParametersProblem`](@ref LineCableModels.Engine.LineParametersProblem)
    accepts completed designs and their placements directly. It builds the
    corresponding `LineCableSystem` once, then adds operating temperature,
    earth properties, and analysis frequencies.
=#

#=
### Earth model

The earth return path significantly affects cable impedance calculations and needs to be properly modeled. In this tutorial, only a basic model with typical soil properties is defined. This will be further elaborated in the subsequent tutorials.
=#

# Define a frequency scan and typical homogeneous-soil properties:
f = collect(10.0 .^ range(0, stop = 6, length = 10)) # 1 Hz to 1 MHz
earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0)

#=
### Three-phase system in trefoil configuration

This section ilustrates the construction of a cable system with three identical cables arranged in a trefoil formation.
=#

# Describe three cables touching in trefoil at 1 m burial depth. The spacing is
# the center-to-center distance:
formation = trefoil(
    x = 0.0,
    y = -1.0,
    spacing = 70e-3
)

# Retain one core phase per cable and ground the sheath and jacket terminals:
connections = [Dict(:core => phase, :sheath => 0, :jacket => 0) for phase in 1:3]
problem = LineParametersProblem(
    fill(loaded_design, 3),
    formation;
    connections,
    system_id = "18kV_1000mm2_trefoil",
    line_length = 1000.0,
    temperature = 20.0,
    earth_props = earth,
    frequencies = f
)
cable_system = problem.system
earth_params = problem.earth_props

# Earth model base properties:
earthmodel_df = DataFrame(earth_params)

#=
!!! note "Phase mapping"
    Each connection mapping uses retained terminal names from `terminal_order`.
    Core assignments `1`, `2`, and `3` select phases A, B, and C. A zero
    assignment marks a terminal for Kron elimination. Reusing one positive
    assignment bundles terminals.
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
    earth_model = earth_params,
    zoom_factor = 2.0,
    display_plot = false, #hide
    controls = false #hide
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
4. Build and preview a three-phase cable system in trefoil arrangement.
5. Export the physical model for PSCAD and ATPDraw.

[`LineCableModels.jl`](@ref) provides detailed routines for power cable modeling
with a physically meaningful representation of its regions and terminal groups. This approach
ensures that electromagnetic parameters are calculated with high precision. Now you can go ahead and run these cable simulations like a boss!
=#
