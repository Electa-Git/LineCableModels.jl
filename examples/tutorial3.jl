#=
# Tutorial 3 - Computing line parameters

This case file demonstrates how to model an armored high-voltage single-core power cable
using the [`LineCableModels.jl`](@ref) package. The objective is to build a complete representation of a single-core 525 kV cable with a 1600 mm² copper conductor, 1.2 mm tubular lead sheath and 68 x 6 mm galvanized steel armor, based on the design described in [Karmokar2025](@cite).
=#

#=
**Tutorial outline**
```@contents
Pages = [
    "tutorial3.md",
]
Depth = 2:3
```
=#

#=
## Introduction
HVDC cables are constructed around a central conductor enclosed by a triple-extruded insulation system (inner/outer semi-conductive layers and main insulation). A metallic screen and protective outer sheath are then applied for land cables. Subsea designs add galvanized steel wire armor over this structure to provide mechanical strength against water pressure. A reference design for a 525 kV HVDC cable [is shown here](https://nkt.widen.net/content/pnwgwjfudf/pdf/Extruded_DC_525kV_DS_EN_DEHV_HV_DS_DE-EN.pdf).
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

# Initialize library and the required materials for this design:
materials = MaterialsLibrary(add_defaults = true)

# Inspect the contents of the materials library:
materials_df = DataFrame(materials)

#=
## Cable dimensions

The cable under consideration is a high-voltage, stranded copper conductor cable with XLPE insulation, water-blocking tape, lead tubular screens, PE inner sheath, PP bedding, steel armor and PP jacket, rated for 525 kV HVDC systems. This information is typically found in the cable datasheet and is based on the design studied in [Karmokar2025](@cite).

The cable is found to have the following configuration:
=#

num_co_wires = 127 # number of core wires
num_ar_wires = 68  # number of armor wires
d_core = 0.0463    # nominal core overall diameter
d_w = 3.6649e-3    # nominal strand diameter of the core (minimum value to match datasheet)
t_sc_in = 2e-3     # nominal internal semicon thickness
t_ins = 26e-3      # nominal main insulation thickness
t_sc_out = 1.8e-3  # nominal external semicon thickness
t_wbt = .3e-3      # nominal thickness of the water blocking tape
t_sc = 3.3e-3      # nominal lead screen thickness
t_pe = 3e-3        # nominal PE inner sheath thickness
t_bed = 3e-3       # nominal thickness of the PP bedding
d_wa = 5.827e-3    # nominal armor wire diameter
t_jac = 10e-3      # nominal PP jacket thickness

d_overall = d_core #hide
layers = [] #hide
push!(layers, ("Conductor", missing, d_overall * 1000)) #hide
d_overall += 2 * t_sc_in #hide
push!(layers, ("Inner semiconductor", t_sc_in * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_ins #hide
push!(layers, ("Main insulation", t_ins * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc_out #hide
push!(layers, ("Outer semiconductor", t_sc_out * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_wbt #hide
push!(layers, ("Swellable tape", t_wbt * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc #hide
push!(layers, ("Lead screen", t_sc * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_pe #hide
push!(layers, ("PE inner sheath", t_pe * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_bed #hide
push!(layers, ("PP bedding", t_bed * 1000, d_overall * 1000)) #hide
d_overall += 2 * d_wa #hide
push!(layers, ("Stranded wire armor", d_wa * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_jac #hide
push!(layers, ("PP jacket", t_jac * 1000, d_overall * 1000)); #hide

# The cable structure is summarized in a table for better visualization, with dimensions in milimiters:
df = DataFrame( #hide
    layer = first.(layers), #hide
    thickness = [ #hide
                 ismissing(t) ? "-" : round(t, sigdigits = 2) for t in getindex.(layers, 2) #hide
                 ], #hide
    diameter = [round(d, digits = 2) for d in getindex.(layers, 3)] #hide
) #hide

#=
## Core and main insulation

The conductor has a central wire and six concentric layers. `Conductor.Stranded`
generates the (1/6/12/18/24/30/36) wire pattern while retaining the specified
lay ratio for the helical corrections.
=#

# Select reusable materials from the library:
copper = Material(materials, :copper)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe)
polyacrylate = Material(materials, :polyacrylate)
lead = Material(materials, :lead)
pp = Material(materials, :pp)
steel = Material(materials, :steel)

# Check the reported strand count and describe the complete conductor:
n_strands = 6 # Strands per layer
n_layers = 6 # Layers of strands
@assert 1 + n_strands * sum(1:n_layers) == num_co_wires
core_conductor = Conductor.Stranded(
    :core;
    layers = n_layers + 1,
    wire_radius = d_w / 2,
    num_wires = n_strands,
    lay_ratio = 11.0,
    material = copper
)

#=
### Inner semiconductor

Inner semiconductor (1000 Ω.m as per IEC 840):
=#

inner_semiconductor = Insulator.Semicon(:core; thickness = t_sc_in, material = semicon1)

#=
### Main insulation

Add the insulation layer:
=#

main_insulation = Insulator.Tubular(:core; thickness = t_ins, material = pe)

#=
### Outer semiconductor

Outer semiconductor (500 Ω.m as per IEC 840):
=#

outer_semiconductor = Insulator.Semicon(:core; thickness = t_sc_out, material = semicon2)

# Water blocking (swellable) tape:
swellable_tape = Insulator.Semicon(:core; thickness = t_wbt, material = polyacrylate)

# Group the declarations associated with the core component:
core_parts = (
    core_conductor,
    inner_semiconductor,
    main_insulation,
    outer_semiconductor,
    swellable_tape
)

cable_id = "525kV_1600mm2"
datasheet_info = (
    designation_code = "(N)2XH(F)RK2Y",
    U0 = 500.0,                        # Phase (pole)-to-ground voltage [kV]
    U = 525.0,                         # Phase (pole)-to-phase (pole) voltage [kV]
    conductor_cross_section = 1600.0,  # [mm²]
    screen_cross_section = 1000.0,     # [mm²]
    resistance = nothing,              # DC resistance [Ω/km]
    capacitance = nothing,             # Capacitance [μF/km]
    inductance = nothing               # Inductance in trifoil [mH/km]
)

#=
### Lead screen/sheath

The lead sheath is described by its radial thickness, followed by the PE inner
sheath and PP bedding. These declarations are stacked over the resolved core.
=#

sheath_parts = (
    Conductor.Tubular(:sheath; thickness = t_sc, material = lead),
    Insulator.Tubular(:sheath; thickness = t_pe, material = pe),
    Insulator.Tubular(:sheath; thickness = t_bed, material = pp)
)

#=
### Armor and outer jacket components

=#

# Describe the armor wires and the PP outer jacket:
lay_ratio = 10.0 # typical value for wire screens
armor_parts = (
    Conductor.Wires(
        :armor;
        wire_radius = d_wa / 2,
        num_wires = num_ar_wires,
        lay_ratio,
        material = steel
    ),
    Insulator.Tubular(:armor; thickness = t_jac, material = pp)
)

# Resolve the complete deterministic cable description:
cable_design = only(CableBuilder(
    cable_id,
    core_parts,
    sheath_parts,
    armor_parts;
    nominal = datasheet_info
))

# Inspect the finished cable design:
plt1 = preview(
    cable_design,
    display_plot = false, #hide
    controls = false #hide
)
plt1.figure #hide

#=
The preview keeps its material scales visible when the side dock is short. If
the complete layer legend does not fit, it ends with `(...)`; resizing an
interactive GLMakie or WGLMakie window restores entries as space becomes
available. SVG export always includes the complete legend.
=#

#=
## Examining the cable parameters (RLC)

=#

# Calculate and summarize the cable constants explicitly:
constants = compute(
    CableConstantsProblem(cable_design),
    Formulation()
)
core_df = DataFrame(constants)

# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)

#=
## Saving the cable design

Load an existing [`CablesLibrary`](@ref) file or create a new one:
=#

library = CablesLibrary()
library_file = fullfile("cables_library.json")
isfile(library_file) && load!(library, file_name = library_file)
add!(library, cable_design)
library_df = DataFrame(library)

# Save to file for later use:
save(library, file_name = library_file);

# Verify that the saved cable can be recovered in a fresh session:
loaded_library = CablesLibrary()
load!(loaded_library, file_name = library_file)
loaded_design = get(loaded_library, cable_id)

#=
## Defining a cable system

=#

#=
### Earth model

Define a static earth model and a logarithmic frequency scan. Earth properties are
declared independently of frequency; they are evaluated when the complete
problem is resolved.
=#

f = collect(10.0 .^ range(0, stop = 6, length = 61)) # 1 Hz to 1 MHz
earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0)

#=
### Underground bipole configuration

=#

# Define the coordinates and phase mapping for both poles:
xp, xn, y0 = -0.5, 0.5, -1.0;
positions = (
    at(x = xp, y = y0, phases = (:core => 1, :sheath => 0, :armor => 0)),
    at(x = xn, y = y0, phases = (:core => 2, :sheath => 0, :armor => 0))
)

# Build the complete bipole line-parameter problem from the loaded design:
problem = only(SystemBuilder(
    "525kV_1600mm2_bipole",
    loaded_design,
    positions;
    length = 1000.0,
    temperature = 20.0,
    earth,
    frequencies = f
))
cable_system = problem.system
earth_params = problem.earth_props

# Inspect the frequency-dependent earth model produced for this problem:
earthmodel_df = DataFrame(earth_params)

#=
### Cable system preview

In this section the complete bipole cable system is examined.
=#

# Display system details:
system_df = DataFrame(cable_system)

# Visualize the cross-section of the three-phase system:
plt2 = preview(
    cable_system,
    earth_model = earth_params,
    zoom_factor = 2.0,
    display_plot = false, #hide
    controls = false #hide
)
plt2.figure #hide

#=
## PSCAD & ATPDraw export
Export to PSCAD input file:
=#

output_file = fullfile("pscad_export.pscx")
export_file = export_data(:pscad, cable_system, earth_params, file_name = output_file);

# Export to ATPDraw project file (XML):
output_file = fullfile("atp_export.xml")
export_file = export_data(:atp, cable_system, earth_params, file_name = output_file);

#=
## Frequency-dependent line parameters

[`Formulation`](@ref) selects the physical and numerical methods. The default
analytical formulation uses the scaled-Bessel internal-impedance method, lossless
insulation impedance/admittance, and the Papadopoulos earth-return methods. The
same formulation value can be reused for ordinary, parametric, or Monte Carlo
execution.
=#

# Define the formulation and run the 1 Hz–1 MHz frequency scan:
formulation = Formulation()
@time line_parameters = compute(
    problem,
    formulation;
    options = (verbosity = (default = 0,),)
);

# Obtain the series and shunt results in per-kilometre units:
series_rl, shunt_gc = DataFrame(
    line_parameters, (R, L, G, C); length_unit = :kilo, tol = 1e-9);

# Display the series resistance and inductance table:
series_rl[1, 1]

# Display the corresponding shunt conductance and capacitance table:
shunt_gc[1, 1]

# Plot the R/L and G/C frequency responses on logarithmic frequency axes. The
# accessors select the displayed quantities; each matrix family occupies one
# page with its two quantities side by side:
rlcg_plots = CairoMakie.plot(
    line_parameters,
    (R, L, G, C);
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 450),
    display_plot = false, #hide
    controls = false #hide
)
rlcg_plots[1].figure #hide

# The shunt conductance and capacitance responses:
rlcg_plots[2].figure #hide

# Plot the real and imaginary parts of the complex Z/Y matrices. Each generated
# page places the two components side by side so their frequency dependence can
# be compared directly:
zy_plots = CairoMakie.plot(
    line_parameters;
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 450),
    display_plot = false, #hide
    controls = false #hide
)
zy_plots[1].figure #hide

# The shunt-admittance page uses the same real/imaginary arrangement:
zy_plots[2].figure #hide

# Export ZY matrices to ATPDraw
output_file = fullfile("ZY_export.xml")
export_file = export_data(
    :atp, line_parameters; file_name = output_file, cable_system);

# Obtain the symmetrical components via Fortescue transformation
Tv, sequence_parameters = Fortescue(tol = 1e-5)(line_parameters);

# Obtain the transformed series and shunt matrices:
series_zy, shunt_zy = DataFrame(
    sequence_parameters; length_unit = :kilo, tol = 1e-9);

# Display the transformed series matrix:
series_zy[1, 1]

# Display the transformed shunt matrix:
shunt_zy[1, 1]

# Obtain the corresponding lumped circuit quantities:
series_rl012, shunt_gc012 = DataFrame(
    sequence_parameters,
    (R, L, G, C);
    length_unit = :kilo,
    tol = 1e-9
);

# Display the sequence-domain series table:
series_rl012[1, 1]

# Display the sequence-domain shunt table:
shunt_gc012[1, 1]

# Plot the sequence-domain R/L and G/C responses:
sequence_plots = CairoMakie.plot(
    sequence_parameters,
    (R, L, G, C);
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 450),
    display_plot = false, #hide
    controls = false #hide
)
sequence_plots[1].figure #hide

# The sequence-domain shunt response:
sequence_plots[2].figure #hide

#=
## Conclusion

This tutorial has demonstrated how to:

1. Describe and preview a detailed armored HVDC cable.
2. Save and reload the cable design.
3. Place two cables in an underground bipole system.
4. Preview and export the physical system for PSCAD and ATPDraw.
5. Compute, tabulate, plot, transform, and export frequency-dependent line parameters.
=#
