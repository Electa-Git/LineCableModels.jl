using Revise
using LineCableModels
using LineCableModels.Engine
using LineCableModels.Engine.Transforms: Fortescue, Levenberg
using LineCableModels.Utils
using LineCableModels.ParametricBuilder: CableBuilder, build, Conductor, Insulator, Material, Earth, SystemBuilder, at, determinize, trifoil
using LineCableModels.UQ: sample, mc, plot
using Measurements
using JSON
using Plots
using Random
using Dates
using JLD2

include("plot_runner.jl")  # Include the plotting functions
using .plot_runner  # Use the module to access its functions

# @load "System_Matrices/results_NA2XS2Y 1x50-16 18-30kV RM.jld2" run_seed
run_seed = abs(rand(Int))
Random.seed!(run_seed)

# --- load dataset ---
filename = "dataset.json"
dataset = [Dict(d) for d in JSON.parsefile(joinpath(@__DIR__, filename))]

# --- helper functions ---
sg(d::Dict, key) = get(d, key, nothing)
scn(x, factor) = (x isa Number) ? x * factor : nothing

set_verbosity!(1);
fullfile(filename) = joinpath(@__DIR__, filename);

# ---Initialize MaterialsLibrary and material uncertainties ---
materials = MaterialsLibrary(add_defaults = true)
material_uncertainties = Dict(
	"aluminum"      => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"copper"        => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"pe"            => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"xlpe"          => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"semicon1"      => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"semicon2"      => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"polyacrylate"  => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
)

# --- corresponding uncertainty vector for cable parameters ---
cable_uncertainties = [
	# --- Cable Dimensions ---
	rand(1.0:20.0), # 1. conductor_wire_diameter (co_d)
	rand(1.0:20.0), # 2. inner_semicon_thickness (t_sc_in)
	rand(1.0:20.0), # 3. insulation_thickness (t_ins)
	rand(1.0:20.0), # 4. outer_semicon_thickness (t_sc_out)
	rand(1.0:20.0), # 5. metallic_screen_wire_diameter (sc_d)
	rand(1.0:20.0), # 6. metallic_screen_tape_thickness (t_cut)
	rand(1.0:20.0), # 7. metallic_screen_tape_width (w_cut)
	nothing, # 8. radial_waterblocking_thickness (t_wbt)
	nothing, # 9. longitudinal_waterblocking_thickness (t_alt)
	rand(1.0:20.0), # 10. outer_sheath_thickness (t_jac)
	nothing, # 11. semicon_tape_thickness (t_sct)
	nothing, # 12. pe_face_thickness (t_pet)
	nothing, # 13. core lay ratio (co_lay)
	nothing, # 14. screen lay ratio (sc_lay)
]

# --- corresponding uncertainty vector for system parameters ---
system_uncertainties = [
	nothing,  # Distance between cables in trifoil formation
	nothing,  # resistivity of earth
	nothing,  # relative permittivity of earth
	nothing,  # relative permeability of earth
	nothing,  # length
	nothing,  # temperature
]

# --- Main Loop over Cable Designs ---
cable_series = [
	"NA2XS2Y 1x50/16 12/20kV RM",
	"NA2XS2Y 1x70/16 12/20kV RM",
	"NA2XS2Y 1x95/16 12/20kV RM",
	"NA2XS2Y 1x120/16 12/20kV RM",
	"NA2XS2Y 1x150/25 12/20kV RM",
	"NA2XS2Y 1x185/25 12/20kV RM",
	"NA2XS2Y 1x240/25 12/20kV RM",
	"NA2XS2Y 1x300/25 12/20kV RM",
	"NA2XS2Y 1x400/35 12/20kV RM",
	"NA2XS2Y 1x500/35 12/20kV RM",
	"NA2XS2Y 1x630/35 12/20kV RM",
	"NA2XS2Y 1x50/16 18/30kV RM",
	"NA2XS2Y 1x70/16 18/30kV RM",
	"NA2XS2Y 1x95/16 18/30kV RM",
	"NA2XS2Y 1x120/16 18/30kV RM",
	"NA2XS2Y 1x150/25 18/30kV RM",
	"NA2XS2Y 1x185/25 18/30kV RM",
	"NA2XS2Y 1x240/25 18/30kV RM",
	"NA2XS2Y 1x300/25 18/30kV RM",
	"NA2XS2Y 1x400/35 18/30kV RM",
	"NA2XS2Y 1x500/35 18/30kV RM",
	"NA2XS2Y 1x630/35 18/30kV RM",
	"NA2XS2Y 1x50/16 6/10kV RM",
	"NA2XS2Y 1x70/16 6/10kV RM",
	"NA2XS2Y 1x95/16 6/10kV RM",
	"NA2XS2Y 1x120/16 6/10kV RM",
	"NA2XS2Y 1x150/25 6/10kV RM",
	"NA2XS2Y 1x185/25 6/10kV RM",
	"NA2XS2Y 1x240/25 6/10kV RM",
	"NA2XS2Y 1x300/25 6/10kV RM",
	"NA2XS2Y 1x400/35 6/10kV RM",
	"NA2XS2Y 1x500/35 6/10kV RM",
	"NA2XS2Y 1x630/35 6/10kV RM",
	# Add more cable codes as needed
]

for code in cable_series
	# --- Nominal data from dataset for sbsific cable ---
	code_map = Dict(d["code"] => d for d in dataset)
	cable_data = code_map[code]
	datasheet_info = NominalData(
		designation_code = cable_data["code"],
		U0 = cable_data["Uo"], U = cable_data["U"],
		conductor_cross_section = cable_data["conductor_cross_section"], screen_cross_section = cable_data["metallic_screen_cross_section"], armor_cross_section = sg(cable_data, "armor_cross_section"),
		resistance = sg(cable_data, "resistance20C"), capacitance = sg(cable_data, "capacitance"), inductance = sg(cable_data, "inductance_trifoil"),
	)

	# --- Cable Parameters (deterministic) ---
	cable_parameters = [
		# --- Cable Dimensions ---
		scn(sg(cable_data, "conductor_wire_diameter"), 1e-3),                 # 1.conductor_wire_diameter (co_d)
		scn(sg(cable_data, "inner_semicon_thickness"), 1e-3),                 # 2.inner_semicon_thickness (t_sc_in)
		scn(sg(cable_data, "insulation_thickness"), 1e-3),                    # 3.insulation_thickness (t_ins)
		scn(sg(cable_data, "outer_semicon_thickness"), 1e-3),                 # 4.outer_semicon_thickness (t_sc_out)
		scn(sg(cable_data, "metallic_screen_wire_diameter"), 1e-3),           # 5.metallic_screen_wire_diameter (sc_d)
		scn(sg(cable_data, "metallic_screen_tape_thickness"), 1e-3),          # 6.metallic_screen_tape_thickness (t_cut)
		scn(sg(cable_data, "metallic_screen_tape_width"), 1e-3),              # 7.metallic_screen_tape_width (w_cut)
		scn(sg(cable_data, "radial_waterblocking_thickness"), 1e-3),          # 8.radial_waterblocking_thickness (t_wbt)
		scn(sg(cable_data, "longitudinal_waterblocking_thickness"), 1e-3),    # 9.longitudinal_waterblocking_thickness (t_alt)
		scn(sg(cable_data, "outer_sheath_thickness"), 1e-3),                  # 10.outer_sheath_thickness (t_jac)
		scn(sg(cable_data, "semicon_tape_thickness"), 1e-3),                  # 11.semicon_tape_thickness (t_sct)
		scn(sg(cable_data, "pe_face_thickness"), 1e-3),                       # 12.pe_face_thickness (t_pet)
		# --- Core and Screen Lay ---
		13.0,                                                                 # 13.core lay ratio (co_lay)
		10.0,                                                                 # 14.screen lay ratio (sc_lay)
		# --- Number of wires and layers ---
		cable_data["conductor_num_wires"],                                    # 15.core number of wires (co_w)
		Int(round(0.5 + sqrt((cable_data["conductor_num_wires"] - 0.25) / 3))),                        # 16.number of layers in hexa pattern (co_n)
		cable_data["metallic_screen_num_wires"],                              # 17.screen number of wires (sc_w)
		1,                                                                    # 18.screen number of layers (sc_n)
	]

	# --- System Parameters (deterministic) --- 
	system_parameters = [
		0.0,                # 1. Distance between cables in trifoil formation
		100.0,            # 2. resistivity of earth
		10.0,            # 3. relative permittivity of earth
		1.0,           # 4. relative permeability of earth
		1000.0,            # 5. length
		20.0,                 # 6. temperature
	]

	# Dynamic Construction Function --- based on layer definitions
	function create_layer(layer_str::String, cable_parameters::Vector, cable_uncertainties::Vector, materials::MaterialsLibrary, material_uncertainties::Dict)
		# Parse layer info
		layer_info = split(layer_str, ',')
		part_name = Symbol(layer_info[2])
		sbsific_class = layer_info[4]
		material_id = String(layer_info[5])

		# Get the material and its uncertainties
		unc = get(material_uncertainties, material_id, (rho = nothing, mu_r = nothing, eps_r = nothing))
		mat = Material(materials, material_id, rho = unc.rho, mu_r = unc.mu_r, eps_r = unc.eps_r)
		isnothing(mat) && (@warn "Unknown material ID: $material_id"; return nothing)

		# --- CONDUCTOR CONSTRUCTORS ---
		if sbsific_class == "Stranded"
			return Conductor.Stranded(
				part_name;
				layers = Int(cable_parameters[16]),                 # co_n (index 16)
				d = (cable_parameters[1], cable_uncertainties[1]),        # co_d (index 1)
				n = 6,
				lay = (cable_parameters[13], cable_uncertainties[13]),    # co_lay (index 13)
				m = mat,
			)
		elseif sbsific_class == "Wires"
			return Conductor.Wires(
				part_name;
				layers = Int(cable_parameters[18]),                 # sc_n (index 18)
				d = (cable_parameters[5], cable_uncertainties[5]),        # sc_d (index 5)
				n = Int(cable_parameters[17]),                      # sc_w (index 17)
				lay = (cable_parameters[14], cable_uncertainties[14]),    # sc_lay (index 14)
				m = mat,
			)
		elseif sbsific_class == "Strip"
			return Conductor.Strip(
				part_name;
				layers = 1,
				t = (cable_parameters[6], cable_uncertainties[6]),        # t_cut (index 6)
				w = (cable_parameters[7], cable_uncertainties[7]),        # w_cut (index 7)
				lay = (cable_parameters[14], cable_uncertainties[14]),    # sc_lay (index 14)
				m = mat,
			)

			# --- INSULATOR CONSTRUCTORS ---
		elseif sbsific_class == "Semicon" || sbsific_class == "Tubular"
			# Determine the index based on material_id
			if material_id == "semicon1"
				idx = 2  # t_sc_in
			elseif material_id == "semicon2"
				idx = 4  # t_sc_out
			elseif material_id == "xlpe"
				idx = 3  # t_ins
			elseif material_id == "pe"
				idx = 10 # t_jac
			else
				@warn "Unknown material ID for insulator: $material_id"
				return nothing
			end

			t_tuple = (cable_parameters[idx], cable_uncertainties[idx])

			if sbsific_class == "Semicon"
				return Insulator.Semicon(part_name; layers = 1, t = t_tuple, m = mat)
			else # sbsific_class == "Tubular"
				return Insulator.Tubular(part_name; layers = 1, t = t_tuple, m = mat)
			end
		end

		@warn "Unknown sbsific class: $sbsific_class"
		return nothing
	end

	# Build the parts list from dataset
	@info("-------------------------------------------")
	@info("Building the cablemodel for cable type: $code")
	layers_to_build = cable_data["layers"]
	parts = [
		create_layer(layer_str, cable_parameters, cable_uncertainties, materials, material_uncertainties)
		for layer_str in layers_to_build
	]

	# Build the CableModel
	cbs = CableBuilder(code, parts; nominal = datasheet_info)
	# Preview the built cable model
	designs = build(cbs)
	design = first(designs)
	# preview(design)

	# frequency grid is deterministic
	f = 50:50:1000

	# Define system center point (underground at 1 m depth) and the trifoil positions
	positions = [
		trifoil(
			y0 = -1.0,
			d = (system_parameters[1], system_uncertainties[1]),  # distance between cables
			phases = (
				:core   => (1, 2, 3),
				:sheath => (0, 0, 0),
				# :jacket => (0, 0, 0), # Ideally this could take notihing if no jacket is present
			),
		),
	]

	# Earth model (uniform, 100 Ωm with uncertainty)
	earth = Earth(rho = (system_parameters[2], system_uncertainties[2]), eps_r = (system_parameters[3], system_uncertainties[3]), mu_r = (system_parameters[4], system_uncertainties[4]))

	# System specification
	@info("Building the system model in trifoil formation...")
	sbs = SystemBuilder("trifoil_case", cbs, positions;
		length = (system_parameters[5], system_uncertainties[5]),
		temperature = (system_parameters[6], system_uncertainties[6]),
		earth = earth,
		f = f,
	)
	# Preview the built system
	sbs_det = determinize(sbs)
	sbs_sample = sample(sbs_det, distribution = :uniform)
	# preview(sbs_sample.system)

	# Define runtime options 
	opts = (
		force_overwrite = true,                    # Overwrite existing files
		save_path = fullfile("lineparams_output"), # Results directory
		verbosity = 0,                             # Verbosity
	);

	# Define the EMT-type model with the sbsified formulations
	F = FormulationSet(:EMT,
		internal_impedance = InternalImpedance.ScaledBessel(),
		insulation_impedance = InsulationImpedance.Lossless(),
		earth_impedance = EarthImpedance.Papadopoulos(),
		insulation_admittance = InsulationAdmittance.Lossless(),
		earth_admittance = EarthAdmittance.Papadopoulos(),
		modal_transform = Transforms.Fortescue(),
		equivalent_earth = EHEM.EnforceLayer(layer = -1),  # Use the last layer as effective earth
		options = opts,
	);

	###################################################################
	# --- MC SIMULATIONS + VISUALISATION ---
	###################################################################
	@info("Running Monte Carlo simulation on system parameters...")
	# --- Monte Carlo on system parameters ---
	@time res_sys = mc(
		sbs, F,
		return_pdf = true,
		return_samples = true,
		tol = 0.02,
		distribution = :normal,
		per_length = true,
	)

	# --- SAVE TO DISK ---
	# Create a descriptive filename, perhaps including a cable ID or timestamp
	timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM")
	safe_cable_id = replace(code, "/" => "-")
	save_path = "System_Matrices/results_$(safe_cable_id)_$(timestamp).jld2"

	@info "Saving results to $(save_path)..."
	# The macro @save automatically saves the variable `res_sys` into the file
	@save save_path res_sys run_seed material_uncertainties cable_uncertainties system_uncertainties

end
