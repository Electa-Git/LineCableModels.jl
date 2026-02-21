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

# Load the package and set up the environment:
using LineCableModels
using LineCableModels.Engine.FEM
using LineCableModels.Engine.Transforms: Fortescue
using DataFrames
using Printf
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_verbosity!(0); #hide
set_backend!(:gl); #hide

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
	diameter = [round(d, digits = 2) for d in getindex.(layers, 3)], #hide
) #hide

#=
## Core and main insulation

Initialize the conductor object and assign the central wire:
=#

material = get(materials, "copper")
core = ConductorGroup(CircStrands(0.0, Diameter(d_w), 1, 0.0, material))

# Add the subsequent layers of wires and inspect the object:
n_strands = 6 # Strands per layer
n_layers = 6 # Layers of strands
for i in 1:n_layers
	add!(core, CircStrands, Diameter(d_w), i * n_strands, 11.0, material)
end
core

#=
### Inner semiconductor

Inner semiconductor (1000 Ω.m as per IEC 840):
=#

material = get(materials, "semicon1")
main_insu = InsulatorGroup(Semicon(core, Thickness(t_sc_in), material))

#=
### Main insulation

Add the insulation layer:
=#

material = get(materials, "pe")
add!(main_insu, Insulator, Thickness(t_ins), material)

#=
### Outer semiconductor

Outer semiconductor (500 Ω.m as per IEC 840):
=#

material = get(materials, "semicon2")
add!(main_insu, Semicon, Thickness(t_sc_out), material)

# Water blocking (swellable) tape:
material = get(materials, "polyacrylate")
add!(main_insu, Semicon, Thickness(t_wbt), material)

# Group core-related components:
core_cc = CableComponent("core", core, main_insu)

cable_id = "525kV_1600mm2"
datasheet_info = NominalData(
	designation_code = "(N)2XH(F)RK2Y",
	U0 = 500.0,                        # Phase (pole)-to-ground voltage [kV]
	U = 525.0,                         # Phase (pole)-to-phase (pole) voltage [kV]
	conductor_cross_section = 1600.0,  # [mm²]
	screen_cross_section = 1000.0,     # [mm²]
	resistance = nothing,              # DC resistance [Ω/km]
	capacitance = nothing,             # Capacitance [μF/km]
	inductance = nothing,              # Inductance in trifoil [mH/km]
)
cable_design = CableDesign(cable_id, core_cc, nominal_data = datasheet_info)

#=
### Lead screen/sheath

Build the wire screens on top of the previous layer:
=#

material = get(materials, "lead")
screen_con = ConductorGroup(Tubular(main_insu, Thickness(t_sc), material))

# PE inner sheath:
material = get(materials, "pe")
screen_insu = InsulatorGroup(Insulator(screen_con, Thickness(t_pe), material))

# PP bedding:
material = get(materials, "pp")
add!(screen_insu, Insulator, Thickness(t_bed), material)

# Group sheath components and assign to design:
sheath_cc = CableComponent("sheath", screen_con, screen_insu)
add!(cable_design, sheath_cc)

#=
### Armor and outer jacket components

=#

# Add the armor wires on top of the previous layer:
lay_ratio = 10.0 # typical value for wire screens
material = get(materials, "steel")
armor_con = ConductorGroup(
	CircStrands(screen_insu, Diameter(d_wa), num_ar_wires, lay_ratio, material))

# PP layer after armor:
material = get(materials, "pp")
armor_insu = InsulatorGroup(Insulator(armor_con, Thickness(t_jac), material))

# Assign the armor parts directly to the design:
add!(cable_design, "armor", armor_con, armor_insu)

# Inspect the finished cable design:
plt1, _ = preview(cable_design)
plt1 #hide

#=
## Examining the cable parameters (RLC)

=#

# Summarize DC lumped parameters (R, L, C):
core_df = DataFrame(cable_design, :baseparams)

# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)

#=
## Saving the cable design

Load an existing [`CablesLibrary`](@ref) file or create a new one:
=#


library = CablesLibrary()
library_file = fullfile("cables_library.json")
load!(library, file_name = library_file)
add!(library, cable_design)
library_df = DataFrame(library)

# Save to file for later use:
save(library, file_name = library_file);

#=
## Defining a cable system

=#

#=
### Earth model 

Define a constant frequency earth model:
=#

f = [1e-3] # Near DC frequency for the analysis
earth_params = EarthModel(f, 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1

# Earth model base (DC) properties:
earthmodel_df = DataFrame(earth_params)

#=
### Underground bipole configuration

=#

# Define the coordinates for both cables:
xp, xn, y0 = -0.5, 0.5, -1.0;

# Initialize the `LineCableSystem` with positive pole:
cablepos = CablePosition(cable_design, xp, y0,
	Dict("core" => 1, "sheath" => 0, "armor" => 0))
cable_system = LineCableSystem("525kV_1600mm2_bipole", 1000.0, cablepos)

# Add the other pole (negative) to the system:
add!(cable_system, cable_design, xn, y0,
	Dict("core" => 2, "sheath" => 0, "armor" => 0))

#=
### Cable system preview

In this section the complete bipole cable system is examined.
=#

# Display system details:
system_df = DataFrame(cable_system)

# Visualize the cross-section of the three-phase system:
plt2, _ = preview(cable_system, earth_model = earth_params, zoom_factor = 2.0)
plt2 #hide

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
## FEM calculations
=#

# Define a LineParametersProblem with the cable system and earth model
problem = LineParametersProblem(
	cable_system,
	temperature = 20.0,  # Operating temperature
	earth_props = earth_params,
	frequencies = f,   # Frequency for the analysis
);

# Estimate domain size based on skin depth in the earth
domain_radius = 10.0; #calc_domain_size(earth_params, f);

# Define custom mesh transitions around each cable
mesh_transition1 = MeshTransition(
	cable_system,
	[1],
	r_min = 0.08,
	r_length = 0.25,
	mesh_factor_min = 0.01 / (domain_radius / 5),
	mesh_factor_max = 0.25 / (domain_radius / 5),
	n_regions = 5)

mesh_transition2 = MeshTransition(
	cable_system,
	[2],
	r_min = 0.08,
	r_length = 0.25,
	mesh_factor_min = 0.01 / (domain_radius / 5),
	mesh_factor_max = 0.25 / (domain_radius / 5),
	n_regions = 5);

# Define runtime options 
opts = (
	force_remesh = true,                # Force remeshing
	force_overwrite = true,             # Overwrite existing files
	plot_field_maps = false,            # Do not compute/ plot field maps
	mesh_only = true,                  # Preview the mesh
	save_path = fullfile("fem_output"), # Results directory
	keep_run_files = true,              # Archive files after each run
	verbosity = 1,                      # Verbosity
);

# Define the FEM formulation with the specified parameters
F = FormulationSet(:FEM,
	impedance = Darwin(),
	admittance = Electrodynamics(),
	domain_radius = domain_radius,
	domain_radius_inf = domain_radius * 1.25,
	elements_per_length_conductor = 1,
	elements_per_length_insulator = 2,
	elements_per_length_semicon = 1,
	elements_per_length_interfaces = 5,
	points_per_circumference = 16,
	mesh_size_min = 1e-6,
	mesh_size_max = domain_radius / 5,
	# mesh_transitions = [mesh_transition1,
	# 	mesh_transition2],
	mesh_size_default = domain_radius / 10,
	mesh_algorithm = 5,
	mesh_max_retries = 20,
	materials = materials,
	options = opts,
);

# # Run the FEM solver
# @time ws, p = compute!(problem, F);

# # Display computation results
# per_km(p, 1; mode = :RLCG, tol = 1e-9)

# # Export ZY matrices to ATPDraw
# output_file = fullfile("ZY_export.xml")
# export_file = export_data(:atp, p; file_name = output_file, cable_system = cable_system);

# # Obtain the symmetrical components via Fortescue transformation
# Tv, p012 = Fortescue(tol = 1e-5)(p);

# # Inspect the transformed matrices
# per_km(p012, 1; mode = :ZY, tol = 1e-9)

# # Or the corresponding lumped circuit quantities
# per_km(p012, 1; mode = :RLCG, tol = 1e-9)

using LineCableModels.Utils: to_nominal
using LineCableModels.DataModel: calc_circstrands_coords

function collect_layer_geometry(design::CableDesign;
	x_offset::Real = 0.0,
	y_offset::Real = 0.0,
)
	layer_no = Ref(0)
	out = Vector{
		NamedTuple{(:layer, :type, :element, :x, :y, :radius_out, :rin, :rext),
			Tuple{Int, String, Int, Float64, Float64, Float64, Float64, Float64}},
	}()

	function push_layer!(layer, x0, y0)
		if layer isa CircStrands
			layer_no[] += 1
			rwire = to_nominal(layer.radius_wire)
			lay_r = layer.num_wires == 1 ? 0.0 : to_nominal(layer.radius_in)
			coords = calc_circstrands_coords(layer.num_wires, rwire, lay_r; C = (x0, y0))
			for (i, (x, y)) in enumerate(coords)
				push!(out, (layer_no[], "circstrands", i, x, y, rwire, NaN, NaN))
			end
		elseif layer isa Strip || layer isa Tubular || layer isa Semicon ||
			   layer isa Insulator
			layer_no[] += 1
			rin = to_nominal(layer.radius_in);
			rex = to_nominal(layer.radius_ext)
			typ = lowercase(string(nameof(typeof(layer))))
			push!(out, (layer_no[], typ, 1, x0, y0, rex, rin, rex))
		elseif layer isa ConductorGroup
			for sub in layer.layers
				push_layer!(sub, x0, y0)
			end
		else
			@warn "Skipping layer $(typeof(layer)) in geometry collection"
		end
	end

	for comp in design.components
		for layer in comp.conductor_group.layers
			push_layer!(layer, x_offset, y_offset)
		end
		for layer in comp.insulator_group.layers
			push_layer!(layer, x_offset, y_offset)
		end
	end
	return out
end

function dump_layer_geometry_csv(design::CableDesign;
	x_offset::Real = 0.0,
	y_offset::Real = 0.0,
	path::AbstractString = "layer_dump.csv",
)
	layer_no = Ref(0)

	function emit!(io, layer, x0, y0)
		if layer isa CircStrands
			layer_no[] += 1
			rwire = to_nominal(layer.radius_wire)
			lay_r = layer.num_wires == 1 ? 0.0 : to_nominal(layer.radius_in)
			coords = calc_circstrands_coords(layer.num_wires, rwire, lay_r; C = (x0, y0))
			for (i, (x, y)) in enumerate(coords)
				println(io, "$(layer_no[]),circstrands,$i,$x,$y,$rwire,,")
			end
		elseif layer isa Strip || layer isa Tubular || layer isa Semicon ||
			   layer isa Insulator
			layer_no[] += 1
			rin = to_nominal(layer.radius_in);
			rex = to_nominal(layer.radius_ext)
			typ = lowercase(string(nameof(typeof(layer))))
			println(io, "$(layer_no[]),$typ,1,$x0,$y0,$rex,$rin,$rex")
		elseif layer isa ConductorGroup
			for sub in layer.layers
				emit!(io, sub, x0, y0)
			end
		else
			@warn "Skipping layer $(typeof(layer)) in CSV dump"
		end
	end

	open(path, "w") do io
		println(io, "layer,type,element,x0,y0,radius_out,rin,rext")
		for comp in design.components
			for layer in comp.conductor_group.layers
				emit!(io, layer, x_offset, y_offset)
			end
			for layer in comp.insulator_group.layers
				emit!(io, layer, x_offset, y_offset)
			end
		end
	end
	return path
end


# dump_layer_geometry_csv(cable_design; path = fullfile("cable_layers.csv"))


all_wires = collect_layer_geometry(cable_design)

using Gmsh: gmsh

gmsh.initialize()

system_id = "test"
gmsh.model.add(system_id)

gmsh.option.set_number("General.InitialModule", 0)
gmsh.option.set_string("General.DefaultFileName", system_id * ".geo")

# Define verbosity level
gmsh.option.set_number("General.Verbosity", 1)

# Set OCC model healing options
gmsh.option.set_number("Geometry.AutoCoherence", 1)
gmsh.option.set_number("Geometry.OCCFixDegenerated", 1)
gmsh.option.set_number("Geometry.OCCFixSmallEdges", 1)
gmsh.option.set_number("Geometry.OCCFixSmallFaces", 1)
gmsh.option.set_number("Geometry.OCCSewFaces", 1)
gmsh.option.set_number("Geometry.OCCMakeSolids", 1)

# Log settings based on verbosity
@info "Initialized Gmsh model: $system_id"

function my_disk!(x, y, r, lcp, wires, llwires)
	cen = gmsh.model.geo.addPoint(x, y, 0.0, lcp)
	p1 = gmsh.model.geo.addPoint(x+r, y, 0.0, lcp)
	p2 = gmsh.model.geo.addPoint(x, y+r, 0.0, lcp)
	p3 = gmsh.model.geo.addPoint(x-r, y, 0, lcp)
	p4 = gmsh.model.geo.addPoint(x, y-r, 0, lcp)

	c1 = gmsh.model.geo.addCircleArc(p1, cen, p2)
	c2 = gmsh.model.geo.addCircleArc(p2, cen, p3)
	c3 = gmsh.model.geo.addCircleArc(p3, cen, p4)
	c4 = gmsh.model.geo.addCircleArc(p4, cen, p1)

	ll = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
	s = gmsh.model.geo.addPlaneSurface([ll])
	push!(wires, s)
	push!(llwires, ll)
	return s
end

lc = 0.0
z0 = 0.0
i=1
wires = []
llwires = []
for wire in all_wires
	x = wire.x
	y = wire.y
	r = wire.radius_out
	my_disk!(x, y, r, lc, wires, llwires)
end

gmsh.model.geo.synchronize()

gmsh.option.set_number("Geometry.SurfaceLabels", 0)  # Show surface labels
gmsh.option.set_number("Geometry.PointNumbers", 0)
gmsh.option.set_number("Geometry.CurveNumbers", 0)
gmsh.option.set_number("Geometry.SurfaceNumbers", 0)
gmsh.option.set_number("Geometry.NumSubEdges", 160)
gmsh.option.set_number("Geometry.Points", 1)
gmsh.option.set_number("Geometry.Curves", 1)
gmsh.option.set_number("Geometry.Surfaces", 0)
gmsh.option.set_number("Mesh.ColorCarousel", 2)  # Colors by physical group
gmsh.option.set_number("Mesh.LineWidth", 1)
gmsh.option.set_number("Mesh.SurfaceFaces", 1)

gmsh.fltk.initialize()

# Define event check function
function check_for_event()
	action = gmsh.onelab.get_string("ONELAB/Action")
	if length(action) > 0 && action[1] == "check"
		gmsh.onelab.set_string("ONELAB/Action", [""])
		@debug "UI interaction detected"
		gmsh.graphics.draw()
	end
	return true
end

# Wait for user to close the window
while gmsh.fltk.is_available() == 1 && check_for_event()
	gmsh.fltk.wait()
end

gmsh.finalize()