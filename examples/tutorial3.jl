#=
# Tutorial 3 - Computing line parameters

Compute frequency-dependent line parameters for a 525 kV cable with a
1600 mm² copper conductor, a 3.3 mm lead sheath, and 68 galvanized steel armor
wires of 5.827 mm diameter. The construction is based on [Karmokar2025](@cite);
the dimensions used in this calculation are listed below.

The example covers cable construction, cable constants, an underground bipole,
external-tool export, frequency-dependent parameters, and modal transformation.
Tables and plots are requested directly from the computed results.
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

HVDC cables have a central conductor enclosed by an insulation system comprising
inner and outer semiconductive layers and the main insulation. Metallic screens
and protective sheaths surround the insulation. Subsea designs also include
steel wire armor for mechanical protection and tensile strength. A reference
525 kV design is available in the
[manufacturer's datasheet](https://nkt.widen.net/content/pnwgwjfudf/pdf/Extruded_DC_525kV_DS_EN_DEHV_HV_DS_DE-EN.pdf).

The reference construction is described with XLPE main insulation. The numerical
example supplied with this tutorial uses the library material `:pe` for the main
insulation and the PE inner sheath. That material choice is retained below; the
reference's XLPE designation does not change the properties selected by the code.
=#

#=
## Getting started

Load the modeling and reporting API, CairoMakie for figures in the documentation,
and DataFrames for the construction-dimensions table.
=#

using LineCableModels
import LineCableModels: homogenize
import CairoMakie
using DataFrames: DataFrame

fullfile(filename) = joinpath(@__DIR__, filename); #hide

# Initialize and inspect the material library:
materials = MaterialsLibrary(add_defaults = true)

#=
## Cable dimensions

The cable consists of a stranded copper conductor, semiconductive layers and
main insulation, water-blocking tape, a lead sheath, a PE inner sheath, PP
bedding, steel armor, and a PP jacket. All dimensions in the declarations below
are in metres.
=#

num_ar_wires = 68  # Number of armor wires.
d_core = 0.0463    # Finished overall core diameter.
d_w = 3.6649e-3    # Source strand diameter used to match the specified core.
t_sc_in = 2e-3     # Inner semiconductor thickness.
t_ins = 26e-3      # Main insulation thickness.
t_sc_out = 1.8e-3  # Outer semiconductor thickness.
t_wbt = 0.3e-3     # Water-blocking tape thickness.
t_sc = 3.3e-3      # Lead sheath thickness.
t_pe = 3e-3        # PE inner sheath thickness.
t_bed = 3e-3       # PP bedding thickness.
d_wa = 5.827e-3    # Armor wire diameter.
t_jac = 10e-3;     # PP outer jacket thickness.

layer_names = ( #hide
    "Conductor", "Inner semiconductor", "Main insulation", #hide
    "Outer semiconductor", "Swellable tape", "Lead sheath", #hide
    "PE inner sheath", "PP bedding", "Stranded wire armor", "PP jacket" #hide
) #hide
layer_thicknesses = ( #hide
    missing, t_sc_in, t_ins, t_sc_out, t_wbt, t_sc, t_pe, t_bed, d_wa, t_jac #hide
) #hide
radial_increments = ( #hide
    0.0, t_sc_in, t_ins, t_sc_out, t_wbt, t_sc, t_pe, t_bed, d_wa, t_jac #hide
) #hide
layer_diameters = d_core .+ 2 .* cumsum(radial_increments) #hide

# Summarize the nominal radial dimensions in millimetres:
cable_dimensions = DataFrame(
    "layer" => collect(layer_names),
    "thickness [mm]" => [ismissing(t) ? missing : round(1000t, sigdigits = 2)
     for t in layer_thicknesses],
    "diameter [mm]" => collect(round.(1000 .* layer_diameters, digits = 2))
)

#=
## Core and main insulation

The source wire diameter and finished core boundary define the stranded core.
[`stranded`](@ref) preserves the source wire areas during compaction and retains
the common lay specification used for helical corrections.
=#

# Select the materials used throughout the cable:
copper = Material(materials, :copper)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe)
polyacrylate = Material(materials, :polyacrylate)
lead = Material(materials, :lead)
pp = Material(materials, :pp)
steel = Material(materials, :steel);

stranded_core = stranded(
    copper;
    shape = Disk(d_w / 2),
    lay = LayRatio(11.0),
    compact = true,
    boundary = Disk(d_core / 2)
);

#=
### Semiconductive and insulating layers

The inner semiconductor, main insulation, and outer semiconductor are specified
by `t_sc_in`, `t_ins`, and `t_sc_out`. The material records are `semicon1`, `pe`,
and `semicon2`, respectively. The source example lists nominal semiconductor
resistivities of 1000 Ω·m and 500 Ω·m, with a reference to IEC 840; the calculation
uses the properties in the selected library records.

Water-blocking tape follows the outer semiconductor. The complete declaration
below states these layers in radial order rather than constructing intermediate
cables for each layer.
=#

#=
### Catalogue information

Catalogue information is retained separately from the physical construction.
The following values are those supplied for this example.
=#

cable_id = "525kV_1600mm2"
datasheet_info = DatasheetInfo(
    designation_code = "(N)2XH(F)RK2Y",
    U0 = 500.0,                        # Pole-to-ground rating [kV].
    U = 525.0,                         # Pole-to-pole rating [kV].
    conductor_cross_section = 1600.0,  # [mm²].
    screen_cross_section = 1000.0,     # [mm²].
    resistance = nothing,              # DC resistance [Ω/km].
    capacitance = nothing,             # Capacitance [μF/km].
    inductance = nothing               # Inductance [mH/km].
)

#=
### Lead sheath, armor, and outer jacket

The lead sheath is followed by the PE inner sheath and PP bedding. The armor
contains 68 steel wires with `LayRatio(10)`; the PP jacket encloses the armor.
The thickness-based layers take their inner boundary from the preceding region.

The core, lead sheath, and armor are separate terminals. Their connections are
assigned when the cable is placed in a system.
=#

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

# Inspect the completed physical design:
cable_design

# Display the cross-section with its material scales:
cable_preview = preview(
    cable_design;
    backend = :cairo,
    legend_overflow = :show_all,
    display_plot = false, #hide
    controls = false #hide
)
cable_preview.figure #hide

#=
## Cable constants and equivalent design

Calculate the cable constants and select R/L and G/C with `values`.
`report(constants)` uses the available quantities and default display units.

The four separate quantity tables are grouped under `constants`.
Each cable-constant table contains the operating frequency and named assembly
columns; quantities with different units are not combined into one table.
=#

constants = CableConstants(cable_design);
rlgc = (R, L, G, C)

constants_report = report(constants; values = rlgc)

# The report exposes the grouped ordinary DataFrames:
constants_report.tables

#=
An equivalent physical design is a separate operation from either computation
or reporting. Request it explicitly with `homogenize`:
=#

equivalent_design = homogenize(cable_design; new_id = cable_id * "_equivalent")

#=
## Saving the cable design

Load an existing [`CablesLibrary`](@ref) or create one, then add the design and
its catalogue information. This saves the physical declaration for reuse; it
does not save a plotting recipe or a computed frequency scan.
=#

library = CablesLibrary()
library_file = fullfile("cables_library.json")
isfile(library_file) && load!(library; file_name = library_file);
add!(library, cable_design; catalogue = datasheet_info);
library

# Write the library:
save(library; file_name = library_file);

# Recover the cable through a fresh library object:
loaded_library = CablesLibrary()
load!(loaded_library; file_name = library_file)
loaded_design = loaded_library[cable_id];

#=
## Defining a cable system

### Earth model and frequency scan

Use a homogeneous earth with resistivity 100 Ω·m, relative permittivity 10,
and relative permeability 1. Earth properties are declared independently of
the scan. The problem supplies 61 logarithmically spaced frequencies from
1 Hz to 1 MHz.
=#

f = collect(10.0 .^ range(0, stop = 6, length = 61))
earth = homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 1.0);

#=
### Underground bipole configuration

Place the pole cables at horizontal positions −0.5 m and 0.5 m, both at depth
1 m. The core connections are numbered 1 and 2; the sheath and armor connections
are assigned 0. The system length is 1000 m.
=#

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

# Attach temperature, earth properties, and frequencies to the physical system:
problem = LineParametersProblem(
    cable_system;
    temperature = 20.0,
    earth_props = earth,
    frequencies = f
)
earth_params = problem.earth_props;

# Inspect the earth parameters retained by the problem:
earth_params

#=
### Cable system preview

Inspect the completed bipole and its cross-section. The earth model supplies the
background layers; `zoom_factor` controls the initial view around the cables.
=#

cable_system

system_preview = preview(
    cable_system;
    earth_model = earth_params,
    zoom_factor = 2.0,
    backend = :cairo,
    legend_overflow = :show_all,
    display_plot = false, #hide
    controls = false #hide
)
system_preview.figure #hide

#=
## PSCAD and ATPDraw export

Export the physical system and earth parameters using the public exporters.
These calls write input files; they do not launch either external solver.
=#

pscad_file = export_data(
    :pscad, cable_system, earth_params;
    file_name = fullfile("pscad_export.pscx")
);

atp_system_file = export_data(
    :atp, cable_system, earth_params;
    file_name = fullfile("atp_export.xml")
);

#=
## Frequency-dependent line parameters

[`Formulation`](@ref) selects the physical and numerical methods. The native
default uses scaled-Bessel internal impedance, lossless insulation models, and
the default earth-return equations. The formulation value can also be used by
parametric and uncertainty-analysis workflows.
=#

formulation = Formulation()

# Run the frequency scan:
@time line_parameters = compute(
    problem,
    formulation;
    options = (verbosity = (default = 0,),)
);

#=
### Numerical access

The immediate form of `@observe` extracts numerical values for further
calculation. It does not construct a table or figure. Inspect the range of the
computed first self-conductance coefficient before reporting resolution is
applied:
=#

conductance_residual = extrema(@observe line_parameters G[1, 1, :])

#=
### Quantity tables

Use the same R/L and G/C selection declared for the cable constants. The report
contains four separate quantity tables. Each full matrix table has one frequency
column followed by all ordered matrix coefficients, including both off-diagonals.

`length_unit=:kilo` selects per-kilometre reporting units. It does not change
the numerical result or the physical system length.
=#

phase_report = report(
    line_parameters;
    values = rlgc,
    length_unit = :kilo
)

#=
To tabulate a particular coefficient or frequency subset, express that request
with `@observe` in `values`. The following selects the first
self-resistance coefficient at the first twelve frequencies. No knowledge of
the generated DataFrame column names is needed to make the selection.
=#

first_term_report = report(
    line_parameters;
    values = @observe(R[1, 1, 1:12]),
    length_unit = :kilo
);

# Retrieve the already-selected resistance DataFrame from the report:
first_term_report.tables.Z.R

#=
### R/L and G/C plots

Pass the computed result directly to `LineCableModels.plot`. Plotting uses
`ydata` for the same quantity selection that reporting names `values`.

Each requested quantity has its own matrix dashboard. Matrix coordinates
identify subplots; each trace follows that coefficient over frequency. The
plots retain the complete matrices rather than assuming symmetry or discarding
small entries in the renderer.
=#

phase_plots = LineCableModels.plot(
    line_parameters;
    ydata = rlgc,
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 750),
    backend = :cairo,
    display_plot = false, #hide
    controls = false #hide
)
phase_plots[1].figure #hide
phase_plots[2].figure #hide
phase_plots[3].figure #hide
phase_plots[4].figure #hide

#=
### Real and imaginary components of Z and Y

With no `ydata` override, the line-result convenience selects R/X and G/B:
the real and imaginary components of Z and Y. This gives four separate matrix
dashboards, using the same public call and the same reporting-unit choice.
=#

zy_plots = LineCableModels.plot(
    line_parameters;
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 750),
    backend = :cairo,
    display_plot = false, #hide
    controls = false #hide
)
zy_plots[1].figure #hide
zy_plots[2].figure #hide
zy_plots[3].figure #hide
zy_plots[4].figure #hide

#=
### Exporting the computed Z/Y matrices

The ATP exporter also accepts the completed line parameters. Supply the cable
system alongside the numerical result:
=#

atp_parameters_file = export_data(
    :atp, line_parameters;
    file_name = fullfile("ZY_export.xml"),
    cable_system = cable_system
);

#=
## Modal transformation

Compute the default frequency-dependent modal transformation as a separate
operation on the completed phase-domain result. Retain the transformation
operators for numerical inspection.
=#

modal_parameters = compute(
    ModalTransformationProblem(line_parameters),
    ModalTransformationFormulation(:default);
    options = (offdiagonal_tolerance = 1e-5,)
);

Tv = operators(modal_parameters).voltage;

# Read transformed coefficients through the same immediate observation API:
modal_impedance = @observe modal_parameters Z[1, 1, :]
modal_admittance = @observe modal_parameters Y[1, 1, :]

#=
### Modal quantity tables

The transformed result uses the same reporting API and the same full-matrix
R/L and G/C selection. Modal-domain coordinates change the interpretation of
the matrix indices, not the organization of the quantity tables.
=#

modal_report = report(
    modal_parameters;
    values = rlgc,
    length_unit = :kilo
)

# Apply the same coefficient and sample selection used for the phase result:
first_modal_term_report = report(
    modal_parameters;
    values = @observe(R[1, 1, 1:12]),
    length_unit = :kilo
);
first_modal_term_report.tables.Z.R

#=
### Full modal matrix plots

Keep the off-diagonal coefficients in the modal view. Entries reduced to zero
by the observation layer remain visible as zero traces; residual coupling that
survives reporting resolution remains visible in its original matrix position.
The plot does not select the diagonal merely because the result is in
`ModalDomain`.

This uses the same quantity selection and plotting call as the phase-domain
view, so the full transformed matrices remain available for inspection.
=#

modal_plots = LineCableModels.plot(
    modal_parameters;
    ydata = rlgc,
    xscale = :log10,
    length_unit = :kilo,
    fig_size = (1100, 750),
    backend = :cairo,
    display_plot = false, #hide
    controls = false #hide
)
modal_plots[1].figure #hide
modal_plots[2].figure #hide
modal_plots[3].figure #hide
modal_plots[4].figure #hide
