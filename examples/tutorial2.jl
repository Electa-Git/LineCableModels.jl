#=
# Tutorial 2 - Building a cable design

Build an 18/30 kV single-core cable with a nominal 1000 mm² aluminum conductor
and a 35 mm² copper screen. Inspect its geometry, calculate and report its
cable constants, save the design, and arrange three cables in a buried trefoil
system for external-tool export.
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
regions. Their dimensions and material properties determine the cable
constants calculated below.

[`CableConstants`](@ref) evaluates earth-free resistance, inductance,
conductance, and capacitance for each independent concentric assembly at a
specified frequency and temperature. These local quantities are distinct from
the frequency-dependent matrices of an installed cable system, which are the
subject of Tutorial 3. In particular, a resistance evaluated at 50 Hz is not
the same measurement as a datasheet DC resistance.

The main operations in this tutorial are physical construction, preview,
calculation, reporting, library persistence, and system export. The completed
numerical result is passed directly to `report`; the reporting convenience
constructs its detached observation internally.
=#

#=
## Getting started

Load the public modeling and reporting API, CairoMakie for the documentation
figures, and DataFrames for the construction-dimensions table.
=#

using LineCableModels
import LineCableModels: homogenize
import CairoMakie
using DataFrames: DataFrame

fullfile(filename) = joinpath(@__DIR__, filename); #hide

# Initialize and inspect the material library:
materials = MaterialsLibrary(add_defaults = true)
materials

#=
The material library saved in Tutorial 1 can also be loaded:

```julia
load!(materials; file_name = "materials_library.json")
```
=#

#=
## Cable dimensions

The reference construction is an 18/30 kV cable with a stranded aluminum
conductor, XLPE insulation, a concentric copper wire screen, water-blocking
layers, an aluminum moisture barrier, and a PE jacket. Its designation is
associated with HD 620 10C [CENELEC_HD620_S3_2023](@cite) and DIN VDE 0276-620
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

The numerical example uses the library's `:pe` material for the main
insulation as well as the outer PE layers. It also uses `:polyacrylate` for the
tapes identified below as semiconductive and water-blocking tapes. These are
the material assignments used by this example; the construction labels do
not substitute different material properties.

All dimensions in the declarations below are in metres.
=#

num_sc_wires = 49  # Number of copper screen wires.
d_core = 38.1e-3   # Finished overall core diameter.
d_w = 4.7e-3       # Source strand diameter.
t_sc_in = 0.6e-3   # Inner semiconductor thickness.
t_ins = 8e-3       # Main insulation thickness.
t_sc_out = 0.3e-3  # Outer semiconductor thickness.
d_ws = 0.95e-3     # Screen-wire diameter.
t_cut = 0.1e-3     # Copper tape thickness.
w_cut = 10e-3      # Copper tape width.
t_wbt = 0.3e-3     # Water-blocking tape thickness.
t_sct = 0.3e-3     # Semiconductive tape thickness.
t_alt = 0.15e-3    # Aluminum tape thickness.
t_pet = 0.05e-3    # PE facing thickness over the aluminum tape.
t_jac = 2.4e-3;    # Outer PE jacket thickness.

layer_names = ( #hide
    "Conductor", "Inner semiconductive tape", "Inner semiconductor", #hide
    "Main insulation", "Outer semiconductor", "Outer semiconductive tape", #hide
    "Wire screen", "Copper tape", "Water-blocking tape", "Aluminum tape", #hide
    "PE facing", "PE jacket" #hide
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

# Summarize the nominal radial schedule in millimetres. The cable constructor
# resolves the detailed strand and tape geometry from the declarations below.
cable_dimensions = DataFrame(
    "layer" => collect(layer_names),
    "thickness [mm]" => [ismissing(t) ? missing : round(1000t, sigdigits = 2)
     for t in layer_thicknesses],
    "diameter [mm]" => collect(round.(1000 .* layer_diameters, digits = 2))
)

#=
## Cable construction

Declare regions from the center outward. [`@cable`](@ref) completes the physical
design, while [`@terminal`](@ref) marks its electrical terminal groups.

### Materials and stranded core

The wire diameter and the 38.1 mm finished boundary define the stranded core.
The four successive outer strand layers use the lay ratios specified below.
Compaction retains the area of each source wire rather than changing its
conductive area to fill the boundary.
=#

aluminum = Material(materials, :aluminum)
copper = Material(materials, :copper)
polyacrylate = Material(materials, :polyacrylate)
semicon1 = Material(materials, :semicon1)
semicon2 = Material(materials, :semicon2)
pe = Material(materials, :pe);

stranded_core = stranded(
    aluminum;
    shape = Disk(d_w / 2),
    lay = LayRatio(15.0, 13.5, 12.5, 11.0),
    compact = true,
    boundary = Disk(d_core / 2)
);

#=
### Insulation, screen, and jacket

The inner semiconducting layer provides the conductor-to-insulation interface;
a tape lies between it and the stranded conductor. The main insulation is
followed by the outer semiconductor and a second tape. The metallic screen
provides shielding and a fault-current path. Outside the screen, the
water-blocking layer, aluminum barrier, PE facing, and outer jacket complete
the construction.

[`screen`](@ref), [`insulation`](@ref), [`sheath`](@ref), and [`jacket`](@ref)
can state their thickness relative to the preceding outward boundary. The
wire screen instead needs the physical radius of its wire-center locus. The
copper tape retains its rectangular width and thickness when placed around
the preceding boundary.
=#

conductor_outer = d_core / 2
screen_wire_locus = conductor_outer + t_sct + t_sc_in + t_ins +
                    t_sc_out + t_sct + d_ws / 2

# Catalogue information is associated with the library entry. It does not
# replace the physical geometry or the selected numerical material properties.
cable_id = "18kV_1000mm2"
datasheet_info = DatasheetInfo(
    designation_code = "NA2XS(FL)2Y",
    U0 = 18.0,                        # Phase-to-ground voltage [kV].
    U = 30.0,                         # Phase-to-phase voltage [kV].
    conductor_cross_section = 1000.0, # [mm²].
    screen_cross_section = 35.0,      # [mm²].
    resistance = 0.0291,              # DC resistance [Ω/km].
    capacitance = 0.39,               # Capacitance [μF/km].
    inductance = 0.3                  # Inductance in trefoil [mH/km].
)

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

#=
The retained terminal names are `:core`, `:sheath`, and `:jacket`. In this
example, `:sheath` groups the copper screen wires and copper tape, while
`:jacket` identifies the conductive aluminum barrier. The outer PE jacket is
an insulating region, not an additional electrical terminal.
=#

# Inspect the completed design and its cross-section:
cable_design

cable_preview = preview(
    cable_design;
    backend = :cairo,
    legend_overflow = :show_all,
    display_plot = false, #hide
    controls = false #hide
)
cable_preview.figure #hide

#=
## Cable constants

Calculate the local cable constants at the defaults of 50 Hz and 20 °C.
The cable-constant calculation treats the innermost terminal of each concentric
assembly as active and the outward terminals as grounded. It does not yet
calculate the three-cable installation introduced later.
=#

constants = CableConstants(cable_design);
constants

#=
### Quantity tables

Pass the numerical result directly to `report`. The R/L and G/C selection
produces four separate quantity tables. Each table contains the operating
frequency and a named column for each assembly; these are not mutual-coupling
matrices between different cables.

The observation owner performs unit conversion and reporting-resolution
handling before tabulation. Choose Ω/km, mH/km, μS/km, and μF/km explicitly so
the R/L/C units match the catalogue entries used below.
=#

rlgc = (R, L, G, C)

constants_report = report(
    TableReportDefinition(rlgc),
    constants;
    length_unit = :kilo,
    quantity_units = (R = :base, L = :milli, G = :micro, C = :micro)
)

# The report exposes the separate ordinary DataFrames for these cable constants:
constants_report.tables

#=
### Assembly selection

For cable constants, an observation index selects an assembly, not a matrix
coefficient or a frequency sample. To report only the resistance of the first
assembly, state that intent with `@observe` in the report request. The example
has one concentric assembly, named `:core`.
=#

core_resistance_request = @observe R[1]
core_resistance_report = report(
    TableReportDefinition((core_resistance_request,)),
    constants;
    length_unit = :kilo,
    quantity_units = (R = :base,)
)

#=
### Catalogue comparison

The following displays keep the three quantities separate. They inspect the
DataFrames already produced by the report; they do not read numerical result
fields or repeat unit conversion.

The catalogue resistance is a DC value, whereas the calculated resistance is
at 50 Hz. Compare their magnitudes with those different conditions in mind:
this is not an equal-condition error calculation.
=#

# Catalogue DC resistance [Ω/km]:
datasheet_info.resistance

# Calculated resistance [Ω/km]; the frequency column records 50 Hz:
constants_report.tables.constants.R

#=
The catalogue inductance is specified for trefoil. `CableConstants` supplies
the earth-free concentric-assembly value, not a calculation of that trefoil
installation. The two are reference values with different stated scopes, not
a benchmark pair.
=#

# Catalogue trefoil inductance [mH/km]:
datasheet_info.inductance

# Calculated local inductance [mH/km]:
constants_report.tables.constants.L

#=
Inspect the catalogue capacitance alongside the calculated capacitance in the
same displayed units. The supplied catalogue record does not include a
capacitance test frequency or temperature; no such conditions are inferred.
The main insulation in this numerical model uses `:pe` as stated above.
=#

# Catalogue capacitance [μF/km]:
datasheet_info.capacitance

# Calculated local capacitance [μF/km]:
constants_report.tables.constants.C

# The conductance table remains a separate result; no catalogue G is supplied:
constants_report.tables.constants.G

#=
## Homogeneous equivalent design

A homogeneous equivalent is a separate physical design, not a report or an
observation of the cable constants. Request it explicitly when that geometry
is needed; the detailed design remains available for inspection and reuse.
=#

equivalent_design = homogenize(cable_design; new_id = cable_id * "_equivalent")

#=
## Saving the cable design

A [`CablesLibrary`](@ref) stores designs for reuse. Add the physical design and
its associated catalogue information, then save the library. Load it into a
fresh library object to retrieve the design by its identifier.
=#

library = CablesLibrary()
add!(library, cable_design; catalogue = datasheet_info);
library

library_file = fullfile("cables_library.json")
save(library; file_name = library_file);

loaded_library = CablesLibrary()
load!(loaded_library; file_name = library_file)
loaded_design = loaded_library[cable_id];
loaded_library

#=
## Defining a cable system

The installed arrangement is represented by [`LineCableSystem`](@ref). A
[`LineParametersProblem`](@ref LineCableModels.Engine.LineParametersProblem)
adds the operating temperature, earth properties, and analysis frequencies to
that completed physical system. This tutorial prepares the problem and exports
the system; it does not execute a frequency-dependent line-parameter scan.

### Earth model

Use homogeneous earth with resistivity 100 Ω·m, relative permittivity 10, and
relative permeability 1. The declared frequency scan contains ten
logarithmically spaced samples from 1 Hz to 1 MHz. The soil declaration is
independent of the frequency grid.
=#

f = collect(10.0 .^ range(0, stop = 6, length = 10))
earth = homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 1.0);

#=
### Three-phase trefoil

Place three copies of the saved design with 70 mm center-to-center spacing
and the trefoil center at `(0, -1)` m. The declared spacing is retained
explicitly; it is not inferred from the nominal cable diameter.

Core assignments `1`, `2`, and `3` identify the three phases. The copper-screen
and aluminum-barrier terminals are assigned `0`, marking them for grounded
terminal reduction when requested by the calculation. The system length is
1000 m.
=#

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

# Inspect the earth parameters retained by the problem:
earth_params

#=
!!! note "Terminal mapping"
    Connection schedules use the terminal names retained by the design.
    Reusing a nonzero assignment groups terminals for bundle reduction;
    zero assignments identify terminals for grounded-conductor reduction.
    The mapping declares those relationships rather than replacing the
    physical cable regions.

### Cable system preview

Inspect the completed three-phase system and its cross-section. The earth
model supplies the background, and `zoom_factor` controls the initial view.
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

Export the physical system and earth parameters through the public exporters.
These calls write input files for electromagnetic-transient analysis; they do
not launch PSCAD or ATPDraw and do not compute the line-parameter frequency
scan. Tutorial 3 covers the frequency-dependent calculation and its reports
and plots.
=#

pscad_file = export_data(
    :pscad, cable_system, earth_params;
    file_name = fullfile("pscad_export.pscx")
);

atp_file = export_data(
    :atp, cable_system, earth_params;
    file_name = fullfile("atp_export.xml")
);
