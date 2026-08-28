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

The conductor has a central wire and six concentric layers. Repeated `Group`
values with `Ring` patterns retain the (1/6/12/18/24/30/36) individual wires,
and `Helix` records the path needed for the helical corrections.
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
strand_radius = d_w / 2
parts, radius = let
    current_radius = zero(strand_radius)
    items = AbstractCablePart[]
    for layer in 0:n_layers
        count = layer == 0 ? 1 : layer * n_strands
        centre_radius = count == 1 ? zero(current_radius) :
                        current_radius + strand_radius
        push!(items, Group(
            :core,
            Region(Symbol(:core_strands_, layer + 1), Disk(strand_radius), copper);
            pattern = Ring(count; r = centre_radius),
            path = count == 1 ? nothing : Helix(LayRatio(11.0))
        ))
        current_radius = count == 1 ? strand_radius :
                         current_radius + 2strand_radius
    end
    items, current_radius
end

#=
### Inner semiconductor

Inner semiconductor (1000 Ω.m as per IEC 840):
=#

inner_semiconductor = Region(:core_semicon_inner, Shell(t_sc_in), semicon1)

#=
### Main insulation

Add the insulation layer:
=#

main_insulation = Region(:core_insulation, Shell(t_ins), pe)

#=
### Outer semiconductor

Outer semiconductor (500 Ω.m as per IEC 840):
=#

outer_semiconductor = Region(:core_semicon_outer, Shell(t_sc_out), semicon2)

# Water blocking (swellable) tape:
swellable_tape = Region(:core_water_blocking, Shell(t_wbt), polyacrylate)

# Append the dielectric regions in their physical order:
append!(parts, (
    inner_semiconductor, main_insulation, outer_semiconductor, swellable_tape
))
radius += t_sc_in + t_ins + t_sc_out + t_wbt

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

lead_outer = radius + t_sc
push!(parts, Group(
    :sheath,
    Region(:sheath_lead_screen, Annulus(radius, lead_outer), lead)
))
push!(parts, Region(:sheath_inner, Shell(t_pe), pe))
push!(parts, Region(:sheath_bedding, Shell(t_bed), pp))
radius = lead_outer + t_pe + t_bed

#=
### Armor and outer jacket regions

=#

# Describe the armor wires and the PP outer jacket:
lay_ratio = 10.0 # typical value for wire screens
armor_radius = d_wa / 2
armor_centre = radius + armor_radius
push!(parts, Group(
    :armor,
    Region(:armor_wires, Disk(armor_radius), steel);
    pattern = Ring(num_ar_wires; r = armor_centre),
    path = Helix(LayRatio(lay_ratio))
))
push!(parts, Region(:armor_jacket, Shell(t_jac), pp))

# Resolve the complete deterministic cable description:
cable_design = CableDesign(
    Stack(parts);
    cable_id,
    nominal_data = NominalData(; datasheet_info...)
)

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

# Calculate the cable constants explicitly:
constants = compute(
    CableConstantsProblem(cable_design),
    Formulation()
)

# Publish detached scientific observations without reaching into result fields:
published_constants = observables(constants, (R, L, C))

# Use the native DataFrames adapter when a complete tabular result is wanted:
constants_table = DataFrame(constants)

# Obtain the radial analytical equivalence of the retained terminals:
terminals_df = DataFrame(cable_design, :terminals)

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

# Define the coordinates and terminal assignments for both poles:
xp, xn, y0 = -0.5, 0.5, -1.0;
positions = [Pose2(xp, y0), Pose2(xn, y0)]
connections = [
    Dict(:core => 1, :sheath => 0, :armor => 0),
    Dict(:core => 2, :sheath => 0, :armor => 0)
]

# Materialize the physical system before adding calculation-only state:
cable_system = LineCableSystem(
    fill(loaded_design, 2);
    positions,
    connections,
    system_id = "525kV_1600mm2_bipole",
    line_length = 1000.0
)
problem = LineCableModels.Engine.LineParametersProblem(
    cable_system;
    temperature = 20.0,
    earth_props = earth,
    frequencies = f
)
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

# The lossless insulation model makes shunt conductance mathematically zero.
# Inspect the native floating-point residual without display clipping:
conductance_residual = extrema(@observe line_parameters G[1, 1, :])

# Obtain one wide table. Frequency and matrix coordinates come from the result;
# each observable request adds its own quantity column:
rlgc_table = DataFrame(
    line_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
    length_unit = :kilo
);

# Select the first matrix term with ordinary DataFrames transformations:
first_term_table = subset(
    rlgc_table,
    :row => ByRow(==(1)),
    :column => ByRow(==(1))
)
first(first_term_table, 12)

# Plot the R/L and G/C frequency responses on logarithmic frequency axes. The
# accessors select the displayed quantities; each matrix family occupies one
# page with its two quantities side by side:
rlcg_plots = CairoMakie.plot(
    line_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
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

# Read one transformed Z/Y term through the same observation boundary:
sequence_impedance = @observe sequence_parameters Z[1, 1, :]
sequence_admittance = @observe sequence_parameters Y[1, 1, :]

# Obtain the complete transformed quantities through the native table adapter:
sequence_table = DataFrame(
    sequence_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
    length_unit = :kilo
);

# Display a compact slice with ordinary DataFrames transformations:
first_sequence_term = subset(
    sequence_table,
    :row => ByRow(==(1)),
    :column => ByRow(==(1))
)
first(first_sequence_term, 12)

# Plot the sequence-domain R/L and G/C responses:
sequence_plots = CairoMakie.plot(
    sequence_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
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
