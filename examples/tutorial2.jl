#=
# Tutorial 2 - Building a cable design

Build an 18/30 kV single-core cable with a 1000 mm² aluminum conductor and
a 35 mm² copper screen, calculate its cable constants, and place three cables
in a line system.
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

The cable consists of concentric conductive, semiconducting, and insulating
regions. Their dimensions and material properties determine the resistance,
inductance, and capacitance calculated below. `CableConstants` evaluates
these quantities at a specified frequency and temperature; its resistance
is not necessarily the DC resistance.

This tutorial covers:

1. Building a [`CableDesign`](@ref) from physical regions and placement rules.
2. Examining the resolved geometry and the retained electrical terminals.
3. Calculating the cable's resistance, inductance, and capacitance.
4. Saving the design to a [`CablesLibrary`](@ref) for future use.
5. Placing designs in a [`LineCableSystem`](@ref) and exporting the system for EMT analysis.
=#

#=
## Getting started
=#

# Load the public modeling API and the packages used for presentation:
using LineCableModels
import LineCableModels: homogenize
import CairoMakie
using DataFrames
fullfile(filename) = joinpath(@__DIR__, filename); #hide

# Initialize materials library with default values:
materials = MaterialsLibrary(add_defaults = true)
materials

#=
```julia
# Alternatively, it can be loaded from the example file built in the previous tutorial:
load!(materials, file_name = "materials_library.json")
```
=#

#=
## Cable dimensions

The 18/30 kV cable has a stranded aluminum conductor, XLPE insulation,
concentric copper wire screen, water-blocking tape, and PE jacket.
Its designation describes these construction features under HD 620 10C
[CENELEC_HD620_S3_2023](@cite) and DIN VDE 0276-620
[VDE_DIN_VDE_0276_620_2024](@cite):

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

# The example uses the following dimensions:
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
t_jac = 2.4e-3;    # nominal PE jacket thickness

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

# The cable structure is summarized in a row-wise table with dimensions in millimeters:
cable_dimensions = DataFrame(
    "layer" => collect(layer_names),
    "thickness [mm]" => [ismissing(t) ? missing : round(1000t, sigdigits = 2)
     for t in layer_thicknesses],
    "diameter [mm]" => collect(round.(1000 .* layer_diameters, digits = 2))
)

#=
## Describing the cable

Declare the regions in order from the cable center outward.
[`@terminal`](@ref) groups conductive regions into an electrical terminal;
[`@cable`](@ref) constructs the cable design.
=#

#=
## Core and main insulation

The wire diameter and 38.1 mm finished core diameter determine how many
complete strand layers fit. Compaction preserves each wire's area. Four outer
layers fit, with `6k` wires in layer `k` and the lay ratios specified below.
=#

# Select reusable materials from the library:
aluminum = Material(materials, :aluminum)
copper = Material(materials, :copper)
polyacrylate = Material(materials, :polyacrylate)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe);

# Specify the wire shape, lay ratios, and finished core diameter:
stranded_core = stranded(
    aluminum;
    shape = Disk(d_w / 2),
    lay = LayRatio(15.0, 13.5, 12.5, 11.0),
    compact = true,
    boundary = Disk(d_core / 2)
);

#=
### Inner semiconductor

The inner semiconducting layer smooths the conductor–insulation interface
and reduces electric-field concentrations around individual strands.
Semiconducting tape covers the stranded core beneath this layer.
=#

#=
!!! tip "Physical order"
    [`insulation`](@ref), [`screen`](@ref), [`sheath`](@ref), and
    [`jacket`](@ref) state their physical thickness. Construction resolves each
    one against the preceding outward boundary.
=#

#=
### Main insulation

The main insulation separates the conductor from the outer semiconducting
layer. This example assigns the library's PE material to that region.
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

# The wire screen needs its physical center locus. The copper tape retains its
# measured rectangular width and thickness; placement bends it around the
# preceding cable boundary without changing its cross-sectional area.
conductor_outer = d_core / 2
screen_wire_locus = conductor_outer + t_sct + t_sc_in + t_ins +
                    t_sc_out + t_sct + d_ws / 2

# Keep catalog data beside the physical model rather than inside it:
cable_id = "18kV_1000mm2"
datasheet_info = DatasheetInfo(
    designation_code = "NA2XS(FL)2Y",
    U0 = 18.0,                        # Phase-to-ground voltage [kV]
    U = 30.0,                         # Phase-to-phase voltage [kV]
    conductor_cross_section = 1000.0, # [mm²]
    screen_cross_section = 35.0,      # [mm²]
    resistance = 0.0291,              # DC resistance [Ω/km]
    capacitance = 0.39,               # Capacitance [μF/km]
    inductance = 0.3                  # Inductance in trefoil [mH/km]
)

# Declare and complete the cable in one outward block. `@terminal` is used only
# where the enclosed conductive regions form one electrical terminal.
cable_design = @cable cable_id begin
    @terminal :core begin
        stranded_core
        screen(polyacrylate; t = t_sct)
        screen(semicon1; t = t_sc_in)
        insulation(pe; t = t_ins)
        screen(semicon2; t = t_sc_out)
        screen(polyacrylate; t = t_sct)
    end
    @terminal :sheath begin
        wires(
            copper;
            shape = Disk(d_ws / 2),
            n = num_sc_wires,
            r = screen_wire_locus,
            lay = LayRatio(10),
            tag = :screen_wire
        )
        tape(
            copper;
            section = Rectangle(w_cut, t_cut),
            n = 1,
            lay = LayRatio(10),
            tag = :copper_tape
        )
    end
    screen(polyacrylate; t = t_wbt)
    @terminal :jacket begin
        sheath(aluminum; t = t_alt)
    end
    jacket(pe; t = t_pet)
    jacket(pe; t = t_jac)
end;
cable_library = CablesLibrary()
add!(cable_library, cable_design; catalogue = datasheet_info);

# Inspect the one completed physical design:
cable_plot = preview(
    cable_design,
    display_plot = false, #hide
    controls = false #hide
)
cable_plot.figure #hide

#=
## Examining the cable parameters (RLC)

Calculate resistance, inductance, and capacitance, then compare them with the
datasheet values.
=#

# Calculate the cable constants at the defaults of 50 Hz and 20 °C:
constants = CableConstants(cable_design);
constants

# Select R, L, and C and convert them to a table:
constants_table = DataFrame(observables(constants, (R, L, C)))

# Construct the homogeneous equivalent cable design:
equivalent_design = homogenize(cable_design; new_id = cable_id * "_equivalent")
equivalent_summary = equivalent_design

# `observables` publishes detached values in the units conventionally used by
# cable manufacturers. Compare the calculated and datasheet values:
published_constants = observables(constants, (R, L, C));
datasheet_comparison = DataFrame(
    source = ("calculated", "datasheet"),
    R = (published_constants[1].values, datasheet_info.resistance),
    L = (published_constants[2].values, datasheet_info.inductance),
    C = (published_constants[3].values, datasheet_info.capacitance)
)
comparison_units = map(payload -> payload.unit, published_constants);

# Inspect the completed physical design through its bounded Base display:
cable_design

#=
## Saving the cable design

!!! note "Cables library"
    Designs can be saved to a library for future use. The [`CablesLibrary`](@ref)
    stores multiple cable designs and is managed through [`add!`](@ref),
    ordinary collection operations, and
    [`save`](@ref LineCableModels.ImportExport.save).
=#

# Store the cable design and inspect the library contents:
library = CablesLibrary()
add!(library, cable_design);
library

# Save to file for later use:
output_file = fullfile("cables_library.json")
save(library, file_name = output_file);

# Load the saved design into a fresh library and retrieve it by identifier:
loaded_library = CablesLibrary()
load!(loaded_library, file_name = output_file)
loaded_design = loaded_library[cable_id];
loaded_library

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
earth = homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 1.0);

#=
### Three-phase system in trefoil configuration

This section ilustrates the construction of a cable system with three identical cables arranged in a trefoil formation.
=#

# Describe three cables touching in trefoil at 1 m burial depth. The connection
# schedules assign one core phase per cable and ground the metallic screens:
placements = @trefoil loaded_design spacing=70e-3 center=(0.0, -1.0) core=(1, 2, 3) sheath=0 jacket=0
cable_system = build(
    LineCableSystem,
    placements;
    environment = earth,
    system_id = "18kV_1000mm2_trefoil",
    line_length = 1000.0
)
problem = LineParametersProblem(
    cable_system;
    temperature = 20.0,
    earth_props = earth,
    frequencies = f
)
earth_params = problem.earth_props;

# Inspect the static earth declaration:
earth_params

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
cable_system

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

Export the cable system for electromagnetic transient simulations in PSCAD
and ATPDraw.
=#

# Export to PSCAD input file:
output_file = fullfile("pscad_export.pscx")
export_file = export_data(:pscad, cable_system, earth_params, file_name = output_file);

# Export to ATPDraw project file (XML):
output_file = fullfile("atp_export.xml")
export_file = export_data(:atp, cable_system, earth_params, file_name = output_file);
