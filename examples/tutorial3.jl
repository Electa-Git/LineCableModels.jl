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
import LineCableModels: homogenize
import CairoMakie
using DataFrames
using LinearAlgebra: diag
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_backend!(:cairo); #hide

# Initialize library and the required materials for this design:
materials = MaterialsLibrary(add_defaults = true)

# Inspect the contents of the materials library:
materials

#=
## Cable dimensions

The cable under consideration is a high-voltage, stranded copper conductor cable with XLPE insulation, water-blocking tape, lead tubular screens, PE inner sheath, PP bedding, steel armor and PP jacket, rated for 525 kV HVDC systems. This information is typically found in the cable datasheet and is based on the design studied in [Karmokar2025](@cite).

The cable is found to have the following configuration:
=#

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
t_jac = 10e-3;     # nominal PP jacket thickness

layer_names = ( #hide
    "Conductor", "Inner semiconductor", "Main insulation", #hide
    "Outer semiconductor", "Swellable tape", "Lead screen", #hide
    "PE inner sheath", "PP bedding", "Stranded wire armor", "PP jacket" #hide
) #hide
layer_thicknesses = ( #hide
    missing, t_sc_in, t_ins, t_sc_out, t_wbt, t_sc, t_pe, t_bed, d_wa, t_jac #hide
) #hide
radial_increments = ( #hide
    0.0, t_sc_in, t_ins, t_sc_out, t_wbt, t_sc, t_pe, t_bed, d_wa, t_jac #hide
) #hide
layer_diameters = d_core .+ 2 .* cumsum(radial_increments) #hide

# The cable structure is summarized in a row-wise table with dimensions in millimeters:
cable_dimensions = DataFrame(
    "layer" => collect(layer_names),
    "thickness [mm]" => [
        ismissing(t) ? missing : round(1000t, sigdigits = 2)
        for t in layer_thicknesses
    ],
    "diameter [mm]" => collect(round.(1000 .* layer_diameters, digits = 2))
)

#=
## Core and main insulation

The conductor has a central wire and six concentric courses. [`stranded`](@ref)
retains the (1/6/12/18/24/30/36) individual wires and the longitudinal paths
needed for the helical corrections.
=#

# Select reusable materials from the library:
copper = Material(materials, :copper)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe)
polyacrylate = Material(materials, :polyacrylate)
lead = Material(materials, :lead)
pp = Material(materials, :pp)
steel = Material(materials, :steel);

# State the actual population and lay of every noncentral conductor course:
stranded_core = stranded(
    copper;
    shape = Disk(d_w / 2),
    layers = 6,
    n = (6, 12, 18, 24, 30, 36),
    lay = LayRatio.((11.0, 11.0, 11.0, 11.0, 11.0, 11.0))
);

#=
### Inner semiconductor

Inner semiconductor (1000 Ω.m as per IEC 840):
=#

#=
### Main insulation

Add the insulation layer:
=#

#=
### Outer semiconductor

Outer semiconductor (500 Ω.m as per IEC 840):
=#

cable_id = "525kV_1600mm2"
datasheet_info = DatasheetInfo(
    designation_code = "(N)2XH(F)RK2Y",
    U0 = 500.0,                        # Phase (pole)-to-ground voltage [kV]
    U = 525.0,                         # Phase (pole)-to-phase (pole) voltage [kV]
    conductor_cross_section = 1600.0,  # [mm²]
    screen_cross_section = 1000.0,     # [mm²]
    resistance = nothing,              # DC resistance [Ω/km]
    capacitance = nothing,             # Capacitance [μF/km]
    inductance = nothing               # Inductance in trefoil [mH/km]
)

#=
### Lead screen/sheath

The lead sheath is an absolute annulus. The following PE inner sheath and PP
bedding are contextual layers and therefore state only their thickness.
=#

#=
### Armor and outer jacket regions

=#

# Declare and complete the cable in physical order. Contextual conductor and
# insulation layers obtain their inner boundary from the preceding region, so
# the datasheet thicknesses are stated directly.
cable_design = @cable cable_id begin
    @terminal :core begin
        stranded_core
        screen(semicon1; t = t_sc_in)
        insulation(pe; t = t_ins)
        screen(semicon2; t = t_sc_out)
        screen(polyacrylate; t = t_wbt)
    end
    @terminal :sheath begin
        sheath(lead; t = t_sc)
    end
    jacket(pe; t = t_pe)
    bedding(pp; t = t_bed)
    @terminal :armor begin
        armor(
            steel;
            shape = Disk(d_wa / 2),
            n = num_ar_wires,
            lay = LayRatio(10),
            tag = :armor_wire
        )
    end
    jacket(pp; t = t_jac)
end;
cable_library = CablesLibrary()
add!(cable_library, cable_design; catalogue = datasheet_info);

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
constants = CableConstants(cable_design);
constants

# Publish detached scientific observations without reaching into result fields:
published_constants = observables(constants, (R, L, C));

# Construct a table from the explicit detached publication:
constants_table = DataFrame(published_constants)

# Request the homogeneous design explicitly when an equivalent
# cable, rather than a solver input, is the desired result:
equivalent_design = homogenize(cable_design; new_id = cable_id * "_equivalent")
equivalent_summary = equivalent_design

# Inspect the completed physical design through its bounded Base display:
cable_design

#=
## Saving the cable design

Load an existing [`CablesLibrary`](@ref) file or create a new one:
=#

library = CablesLibrary()
library_file = fullfile("cables_library.json")
isfile(library_file) && load!(library, file_name = library_file);
add!(library, cable_design);
library

# Save to file for later use:
save(library, file_name = library_file);

# Verify that the saved cable can be recovered in a fresh session:
loaded_library = CablesLibrary()
load!(loaded_library, file_name = library_file)
loaded_design = get(loaded_library, cable_id);

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
earth = Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0);

#=
### Underground bipole configuration

=#

# Place the two poles and state their terminal assignments at the same surface:
positive_pole = @at loaded_design (-0.5, -1.0) connections = (
    core = 1, sheath = 0, armor = 0)
negative_pole = @at loaded_design (0.5, -1.0) connections = (
    core = 2, sheath = 0, armor = 0)
placements = [positive_pole, negative_pole]
cable_system = build(
    LineCableSystem,
    placements;
    environment = earth,
    system_id = "525kV_1600mm2_bipole",
    line_length = 1000.0
)

# Attach calculation-only state to the completed physical system:
problem = LineParametersProblem(
    cable_system;
    temperature = 20.0,
    earth_props = earth,
    frequencies = f
)
earth_params = problem.earth_props;

# Inspect the static earth declaration attached to the problem:
earth_params

#=
### Cable system preview

In this section the complete bipole cable system is examined.
=#

# Display system details:
cable_system

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
native formulation uses the scaled-Bessel internal-impedance method, lossless
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
rlgc_table = DataFrame(observables(
    line_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
    length_unit = :kilo
));

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
sequence_parameters = compute(
    ModalTransformationProblem(line_parameters),
    ModalTransformationFormulation(:Fortescue; tolerance = 1e-5)
);
Tv = operators(sequence_parameters).voltage;

# Read one transformed Z/Y term through the same observation boundary:
sequence_impedance = @observe sequence_parameters Z[1, 1, :]
sequence_admittance = @observe sequence_parameters Y[1, 1, :]

# Publish and tabulate the complete transformed quantities:
sequence_table = DataFrame(observables(
    sequence_parameters,
    (
        @observe(R[:, :, :]),
        @observe(L[:, :, :]),
        @observe(G[:, :, :]),
        @observe(C[:, :, :])
    );
    length_unit = :kilo
));

# Display a compact slice with ordinary DataFrames transformations:
first_sequence_term = subset(
    sequence_table,
    :row => ByRow(==(1)),
    :column => ByRow(==(1))
)
first(first_sequence_term, 12)

# Plot the two modal R/L and G/C responses. Diagonal observation is explicit:
sequence_plots = CairoMakie.plot(
    sequence_parameters,
    (
        @observe((R, diag)[:, :]),
        @observe((L, diag)[:, :]),
        @observe((G, diag)[:, :]),
        @observe((C, diag)[:, :])
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
