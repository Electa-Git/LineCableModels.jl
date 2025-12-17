using Revise
using LineCableModels
using DataFrames
using LineCableModels.Engine
using LineCableModels.Engine.Transforms: Fortescue, Levenberg
using LineCableModels.Utils
using LineCableModels.ParametricBuilder: CableBuilder, build, Conductor, Insulator, Material, Earth, SystemBuilder, at, determinize, trifoil
using LineCableModels.UQ: sample, mc, plot
using Measurements
using JSON
using Plots
using Random
using Test

include("plot_runner.jl")  # Include the plotting functions
using .plot_runner  # Use the module to access its functions

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
	"aluminum" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"copper" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"pe" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"xlpe" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"semicon1" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"semicon2" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
	"polyacrylate" => (rho = rand(1.0:20.0), mu_r = rand(1.0:20.0), eps_r = rand(1.0:20.0)),
)

# --- Nominal data from dataset for sbsific cable ---
code_map = Dict(d["code"] => d for d in dataset)
code = "NA2XS2Y 1x50/16 18/30kV RM"
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
	scn(sg(cable_data, "conductor_wire_diameter"), 1e-3),                 	# 1.conductor_wire_diameter (co_d)
	scn(sg(cable_data, "inner_semicon_thickness"), 1e-3),                 	# 2.inner_semicon_thickness (t_sc_in)
	scn(sg(cable_data, "insulation_thickness"), 1e-3),                    	# 3.insulation_thickness (t_ins)
	scn(sg(cable_data, "outer_semicon_thickness"), 1e-3),                 	# 4.outer_semicon_thickness (t_sc_out)
	scn(sg(cable_data, "metallic_screen_wire_diameter"), 1e-3),           	# 5.metallic_screen_wire_diameter (sc_d)
	scn(sg(cable_data, "metallic_screen_tape_thickness"), 1e-3),          	# 6.metallic_screen_tape_thickness (t_cut)
	scn(sg(cable_data, "metallic_screen_tape_width"), 1e-3),              	# 7.metallic_screen_tape_width (w_cut)
	scn(sg(cable_data, "radial_waterblocking_thickness"), 1e-3),          	# 8.radial_waterblocking_thickness (t_wbt)
	scn(sg(cable_data, "longitudinal_waterblocking_thickness"), 1e-3),    	# 9.longitudinal_waterblocking_thickness (t_alt)
	scn(sg(cable_data, "outer_sheath_thickness"), 1e-3),                  	# 10.outer_sheath_thickness (t_jac)
	scn(sg(cable_data, "semicon_tape_thickness"), 1e-3),                  	# 11.semicon_tape_thickness (t_sct)
	scn(sg(cable_data, "pe_face_thickness"), 1e-3),                       	# 12.pe_face_thickness (t_pet)
	# --- Core and Screen Lay ---
	13.0,                                                                 	# 13.core lay ratio (co_lay)
	10.0,                                                                 	# 14.screen lay ratio (sc_lay)
	# --- Number of wires and layers ---
	cable_data["conductor_num_wires"],                                    	# 15.core number of wires (co_w)
	Int(round(0.5 + sqrt((cable_data["conductor_num_wires"] - 0.25) / 3))),	# 16.number of layers in hexa pattern (co_n)
	cable_data["metallic_screen_num_wires"],                              	# 17.screen number of wires (sc_w)
	1,                                                                    	# 18.screen number of layers (sc_n)
]

# --- corresponding uncertainty vector for cable parameters ---
cable_uncertainties = [
	# --- Cable Dimensions ---
	rand(1.0:20.0),   # 1. conductor_wire_diameter (co_d)
	rand(1.0:20.0),   # 2. inner_semicon_thickness (t_sc_in)
	rand(1.0:20.0),   # 3. insulation_thickness (t_ins)
	rand(1.0:20.0),   # 4. outer_semicon_thickness (t_sc_out)
	rand(1.0:20.0),   # 5. metallic_screen_wire_diameter (sc_d)
	rand(1.0:20.0),   # 6. metallic_screen_tape_thickness (t_cut)
	rand(1.0:20.0),   # 7. metallic_screen_tape_width (w_cut)
	nothing,   # 8. radial_waterblocking_thickness (t_wbt)
	nothing,   # 9. longitudinal_waterblocking_thickness (t_alt)
	rand(1.0:20.0),   # 10. outer_sheath_thickness (t_jac)
	nothing,   # 11. semicon_tape_thickness (t_sct)
	nothing,   # 12. pe_face_thickness (t_pet)
	nothing,   # 13. core lay ratio (co_lay)
	nothing,   # 14. screen lay ratio (sc_lay)
]

# --- System Parameters (deterministic) ---
system_parameters = [
	0.025,  # Distance between cables in trifoil formation
	100.0,  # resistivity of earth
	10.0,  # relative permittivity of earth
	1.0, # relative permeability of earth
	1000.0,  # length
	20.0,   # temperature
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
@info("Building the cablemodel for code: $code")
layers_to_build = cable_data["layers"]
parts = [
	create_layer(layer_str, cable_parameters, cable_uncertainties, materials, material_uncertainties)
	for layer_str in layers_to_build
]

# Build the CableModel
cbs = CableBuilder(code, parts; nominal = datasheet_info)
# Determinize the CableBuilder
cbs_det = determinize(cbs)
# Sample a design from the determinized range of options
cbs_sample = sample(cbs_det, distribution = :uniform)

# Build designs
designs = build(cbs)
# DataFrames for easy comparison
df1 = DataFrame(designs[1], :baseparams) # DataFrame of the nominal design
df_random = DataFrame(cbs_sample, :baseparams) # DataFrame of a random sampled design

# frequency grid is deterministic
# f = 10.0 .^ range(0, stop = 6, length = 10)
f = 50:50:1000

# Define system center point (underground at 1 m depth) and the trifoil positions
positions = [
	trifoil(
		y0 = -1.0,
		d = (system_parameters[1], system_uncertainties[1]),  # distance between cables
		phases = (
			:core   => (1, 2, 3),
			:sheath => (0, 0, 0),
			# :jacket => (0, 0, 0), # Ideally this could take nothing if no jacket is present
		),
	),
]

# Earth model (uniform, 100 Ωm with uncertainty)
earth = Earth(rho = (system_parameters[2], system_uncertainties[2]), eps_r = (system_parameters[3], system_uncertainties[3]), mu_r = (system_parameters[4], system_uncertainties[4]))

# System specification
sbs = SystemBuilder("trifoil_case", cbs, positions;
	length = (system_parameters[5], system_uncertainties[5]),
	temperature = (system_parameters[6], system_uncertainties[6]),
	earth = earth,
	f = f,
)
# Determinize and sample the system
sbs_det = determinize(sbs)
sbs_sample = sample(sbs_det, distribution = :uniform)

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

Ps = Vector{LineParameters}()
for problem in sbs
	_, p = compute!(problem, F)  #  FormulationSet as before
	push!(Ps, p)
end

###################################################################
# --- MC SIMULATIONS + VISUALISATION ---
###################################################################

@info("Running Monte Carlo simulation on cable parameters...")
# --- Monte Carlo on cable arameters ---
@time res = mc(
	cbs,
	return_pdf = true,
	return_samples = true,
	tol = 0.02,
	distribution = :normal,
)

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

# --- Plots ---
# --- Plots of sampled cable parameters ---
p1 = Plots.histogram(
	res.samples.R,
	bins = 80,
	normalize = :pdf,
	xlabel = "R",
	ylabel = "density",
	label = "R",
)
p2 = Plots.histogram(
	res.samples.L,
	bins = 80,
	normalize = :pdf,
	xlabel = "L",
	ylabel = "density",
	label = "L",
)
p3 = Plots.histogram(
	res.samples.C,
	bins = 80,
	normalize = :pdf,
	xlabel = "C",
	ylabel = "density",
	label = "C",
)

# Single Frequency PDFs of system parameters
freq_indices = [1, 2, 3, 7, 13, 20] # corresponding to entries in f
seq = 1 # 1 = zero sequence, 2 = positive sequence
for freq_index in freq_indices
	freq = freq_index * 50.0
	# PDFs of specific frequency samples from system MC
	p4 = Plots.histogram(
		res_sys.samples.R[seq, seq, freq_index, :],
		bins = 80,
		normalize = :pdf,
		xlabel = "R [Ω/km]",
		ylabel = "density",
		label = "System R₀ at $freq Hz",
	)
	display(p4)
	# fname = "Figures/R_zero_$(Int(freq))Hz.png" # .png extension is defined here
	# savefig(p4, fname)
end

# Fan Chart and Ridgeplots over frequency
# Evolution of zero sequence Resistance for selected frequencies
pr_0 = plot_runner.plot_frequency_ridgeline(
	res_sys,
	:R,
	freq_indices;
	sequence = :zero,
	overlap = 0.5,
	length_unit = :kilo,
	per_length = true,
	nbins = 80,
	step_style = false,
);
# Fan Chart of quantity over full frequency range
fR_0 = plot_runner.plot_frequency_fanchart(
	res_sys,
	:R,
	frequencies = f;
	sequence = :zero,
);

# Evolution of pos. sequence Resistance for selected frequencies
pr_1 = plot_runner.plot_frequency_ridgeline(
	res_sys,
	:R,
	freq_indices;
	sequence = :positive,
	overlap = 0.5,
	length_unit = :kilo,
	per_length = true,
	nbins = 80,
	step_style = false,
);
# Fan Chart of quantity over full frequency range
fR_1 = plot_runner.plot_frequency_fanchart(
	res_sys,
	:R,
	frequencies = f;
	sequence = :positive,
);

# Evolution of zero sequence Inductance for selected frequencies
pl_0 = plot_runner.plot_frequency_ridgeline(
	res_sys,
	:L,
	freq_indices;
	sequence = :zero,
	overlap = 0.5,
	length_unit = :kilo,
	per_length = true,
	nbins = 80,
	step_style = false,
);
# Fan Chart of quantity over full frequency range
fL_0 = plot_runner.plot_frequency_fanchart(
	res_sys,
	:L,
	frequencies = f;
	sequence = :zero,
);

# Evolution of pos. sequence Inductance for selected frequencies
pl_1 = plot_runner.plot_frequency_ridgeline(
	res_sys,
	:L,
	freq_indices;
	sequence = :positive,
	overlap = 0.5,
	length_unit = :kilo,
	per_length = true,
	nbins = 80,
	step_style = false,
);
# Fan Chart of quantity over full frequency range
fL_1 = plot_runner.plot_frequency_fanchart(
	res_sys,
	:L,
	frequencies = f;
	sequence = :positive,
);
