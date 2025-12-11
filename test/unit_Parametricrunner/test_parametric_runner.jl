using Revise
using LineCableModels
using DataFrames
using LineCableModels.Engine
using LineCableModels.Engine.Transforms: Fortescue, Levenberg
using LineCableModels.Utils
using LineCableModels.ParametricBuilder: CableBuilder, build, Conductor, Insulator, Material, Earth, SystemBuilder, at, determinize, trifoil
using LineCableModels.UQ: sample, mc
using Measurements
using JSON
using Plots
using Random
using Test

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
    "aluminum"    	=> (rho = 10.0, mu_r = 12.0, eps_r = 8.0),
    "copper"      	=> (rho = nothing, mu_r = 10.0, eps_r = nothing),
	"pe"          	=> (rho = 10.0,  mu_r = nothing, eps_r = nothing),
    "xlpe"        	=> (rho = 10.0,  mu_r = 10.0, eps_r = nothing),
    "semicon1"   	=> (rho = nothing, mu_r = nothing, eps_r = 10.0),
	"semicon2"    	=> (rho = nothing, mu_r = 10.0, eps_r = nothing),
	"polyacrylate" 	=> (rho = 10.0, mu_r = nothing, eps_r = nothing),
)

# --- Nominal data from dataset for specific cable ---
code_map = Dict(d["code"] => d for d in dataset)
code = "NA2XS2Y 1x50/16 18/30kV RM"
cable_data = code_map[code]
datasheet_info = NominalData(
	designation_code = cable_data["code"],
	U0 = cable_data["Uo"], U = cable_data["U"],
	conductor_cross_section = cable_data["conductor_cross_section"], screen_cross_section = cable_data["metallic_screen_cross_section"], armor_cross_section = sg(cable_data, "armor_cross_section"),
	resistance = sg(cable_data, "resistance20C"), capacitance = sg(cable_data, "capacitance"), inductance = sg(cable_data, "inductance_trifoil")
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
	(cable_data["conductor_num_wires"] - 1)/6 + 1,                        # 16.number of layers in hexa pattern (co_n)
	cable_data["metallic_screen_num_wires"],                              # 17.screen number of wires (sc_w)
	1,                                                                    # 18.screen number of layers (sc_n)
]

# --- corresponding uncertainty vector for cable parameters ---
cable_uncertainties = [
	# --- Cable Dimensions ---
    rand(1.0:20.0), 													  # 1. conductor_wire_diameter (co_d)
    rand(1.0:20.0), 															  # 2. inner_semicon_thickness (t_sc_in)
    rand(1.0:20.0), 															  # 3. insulation_thickness (t_ins)
    rand(1.0:20.0), 															  # 4. outer_semicon_thickness (t_sc_out)
    rand(1.0:20.0), 															  # 5. metallic_screen_wire_diameter (sc_d)
    rand(1.0:20.0), 															  # 6. metallic_screen_tape_thickness (t_cut)
    rand(1.0:20.0), 															  # 7. metallic_screen_tape_width (w_cut)
    nothing, 															  # 8. radial_waterblocking_thickness (t_wbt)
    nothing, 															  # 9. longitudinal_waterblocking_thickness (t_alt)
    rand(1.0:20.0), 															  # 10. outer_sheath_thickness (t_jac)
    nothing, 															  # 11. semicon_tape_thickness (t_sct)
    nothing, 															  # 12. pe_face_thickness (t_pet)
    nothing, 															  # 13. core lay ratio (co_lay)
    nothing, 															  # 14. screen lay ratio (sc_lay)
]

# --- System Parameters (deterministic) ---
system_parameters = [
	0.0,  														# Distance between cables in trifoil formation
	100.0,  													# resistivity of earth
	10.0,  														# relative permittivity of earth
	1.0, 														# relative permeability of earth
	1000.0,  													# length
	20.0,   													# temperature
]

# --- corresponding uncertainty vector for system parameters ---
system_uncertainties = [
	nothing,  													# Distance between cables in trifoil formation
	nothing,  													# resistivity of earth
	nothing,  													# relative permittivity of earth
	nothing,  													# relative permeability of earth
	10.0,  													# length
	10.0,   													# temperature
]


# Dynamic Construction Function --- based on layer definitions
function create_layer(layer_str::String, cable_parameters::Vector, cable_uncertainties::Vector, materials::MaterialsLibrary, material_uncertainties::Dict)
    # Parse layer info
    layer_info = split(layer_str, ',')
    part_name = Symbol(layer_info[2])
    specific_class = layer_info[4]
    material_id = String(layer_info[5])
    
    # Get the material constant
	unc = get(material_uncertainties, material_id, (rho=nothing, mu_r=nothing, eps_r=nothing))
	mat = Material(materials, material_id, rho = unc.rho, mu_r = unc.mu_r, eps_r = unc.eps_r)
    isnothing(mat) && (@warn "Unknown material ID: $material_id"; return nothing)

    # --- CONDUCTOR CONSTRUCTORS ---
    if specific_class == "Stranded"
        # 1,core,Conductor,Stranded,aluminum
        # d uses index 1, lay uses index 13, counts use 15 and 16
        return Conductor.Stranded(
            part_name;
            layers = Int(cable_parameters[16]),                 # co_n (index 16)
            d = (cable_parameters[1], cable_uncertainties[1]),        # co_d (index 1)
            n = 6,											
            lay = (cable_parameters[13], cable_uncertainties[13]),    # co_lay (index 13)
            m = mat,
        )
    elseif specific_class == "Wires"
        # 5,sheath,Conductor,Wires,copper
        # d uses index 5, lay uses index 14, counts use 17 and 18
        return Conductor.Wires(
            part_name;
            layers = Int(cable_parameters[18]),                 # sc_n (index 18)
            d = (cable_parameters[5], cable_uncertainties[5]),        # sc_d (index 5)
            n = Int(cable_parameters[17]),                      # sc_w (index 17)
            lay = (cable_parameters[14], cable_uncertainties[14]),    # sc_lay (index 14)
            m = mat,
        )
    elseif specific_class == "Strip"
        # 6,sheath,Conductor,Strip,copper
        # t uses index 6, w uses index 7, lay uses index 14
        return Conductor.Strip(
            part_name;
            layers = 1,
            t = (cable_parameters[6], cable_uncertainties[6]),        # t_cut (index 6)
            w = (cable_parameters[7], cable_uncertainties[7]),        # w_cut (index 7)
            lay = (cable_parameters[14], cable_uncertainties[14]),    # sc_lay (index 14)
            m = mat,
        )

    # --- INSULATOR CONSTRUCTORS ---
    elseif specific_class == "Semicon" || specific_class == "Tubular"
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
        
        if specific_class == "Semicon"
            return Insulator.Semicon(part_name; layers = 1, t = t_tuple, m = mat)
        else # specific_class == "Tubular"
            return Insulator.Tubular(part_name; layers = 1, t = t_tuple, m = mat)
        end
    end

	@warn "Unknown specific class: $specific_class"
    return nothing
end

# Build the parts list from dataset
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
sample_des = sample(cbs_det, distribution = :uniform)

# Build designs
designs = build(cbs)
# DataFrames for easy comparison
df1 = DataFrame(designs[1], :baseparams) # DataFrame of the nominal design
df_random = DataFrame(sample_des, :baseparams) # DataFrame of a random sampled design

# --- TEST ALL PARAMETERS WITH UNCERTAINTY ---
@testset "Uncertainty propagation" begin
	@info("Testing uncertainty propagation for all cable parameters with assigned uncertainty...")
    for eachindex in 1:length(cable_uncertainties)
		idx = eachindex
        applied_unc = cable_uncertainties[idx]
		if applied_unc !== nothing
			@info("-------------------------------")
			@info("Testing uncertainty of parameter index $idx with uncertainty $applied_unc...")
		end
        isnothing(applied_unc) && continue   # Skip unused indices
        deterministic_value = cable_parameters[idx]
        
        # ============================================================
        # 1. TEST parameter presence in cbs.parts
        # ============================================================
        found = false
        for (part_id, p) in enumerate(cbs.parts)
            for field in fieldnames(typeof(p))
                val = getfield(p, field)
                if val isa Tuple{Float64, Union{Nothing,Float64}} # thickness or diameter case (checks dim)
                    if val[1] == deterministic_value && val[2] == applied_unc
                        found = true
						@info("Found uncertainty in the correct place for parameter index $idx in cbs.parts $part_id, field $field")
                        break
                    end
				elseif val isa Tuple{Int64, Tuple{Float64, Union{Nothing,Float64}}} # lay_ratio case (checks args[2])
					if val[2][1] == deterministic_value && val[2][2] == applied_unc
						found = true
						break
					end
				elseif val isa Tuple{Tuple{Float64, Union{Nothing,Float64}}, Tuple{Float64, Union{Nothing,Float64}}} # width case (checks args[1])
					if val[1][1] == deterministic_value && val[1][2] == applied_unc
						found = true
						@info("Found uncertainty in the correct place for parameter index $idx in cbs.parts $part_id, field $field")
						break
					end
				else
					continue  # skip anything else
                end
            end
            found && break
        end

        @test found || error("Parameter index $idx not found in any cbs.parts field")

        # ============================================================
        # 2. TEST determinized interval
        # ============================================================
		# Wirearray core has multiple layers_to_build
		lower_expected_core = cable_parameters[1] * (1 - cable_uncertainties[1] / 100)
		upper_expected_core = cable_parameters[1] * (1 + cable_uncertainties[1] / 100)

        lower_expected = deterministic_value * (1 - applied_unc / 100)
        upper_expected = deterministic_value * (1 + applied_unc / 100)

        det_found = false
		for (part_id, p) in enumerate(cbs_det.parts)
			for field in fieldnames(typeof(p))
				val = getfield(p, field)
				if val isa Tuple{Tuple{Float64, Float64, Int64}, Nothing} # thickness or diameter case (checks dim)
					lo = val[1][1]
					hi = val[1][2]
					# Check if this field matches the known Core bounds
    				is_core_values = (lo == lower_expected_core && hi == upper_expected_core)
    				# If we found values that look like the Core, but we are NOT currently testing the Core parameter, skip this part.
					if idx != 1 && is_core_values
						continue
					end
					if lo == lower_expected && hi == upper_expected
						det_found = true
						@info("Found determinized bounds [$lower_expected, $upper_expected] for parameter index $idx in cbs_det.parts $part_id, field $field")
						break
					end
				elseif val isa Tuple{Int64, Tuple{Float64, Float64, Int64}} # lay_ratio case (checks args[2])
					lo = val[2][1]
					hi = val[2][2]
					if lo == lower_expected && hi == upper_expected
						det_found = true
						@info("Found determinized bounds [$lower_expected, $upper_expected] for parameter index $idx in cbs_det.parts $part_id, field $field")
						break
					end
				elseif val isa Tuple{Tuple{Tuple{Float64, Float64, Int64}, Nothing}, Tuple{Float64, Nothing}} # width case (checks args[1])
					lo = val[1][1][1]
					hi = val[1][1][2]
					if lo == lower_expected && hi == upper_expected
						det_found = true
						@info("Found determinized bounds [$lower_expected, $upper_expected] for parameter index $idx in cbs_det.parts $part_id, field $field")
						break
					end
				else
					continue  # skip anything else
				end
			end
			det_found && break
		end

        @test det_found || error("Determinized interval for parameter idx=$idx not found")

        # ============================================================
        # 3. TEST sampled designs stay within expected bounds
        # ============================================================
        lower_sample = deterministic_value * (1 - applied_unc * sqrt(3) / 100)
        upper_sample = deterministic_value * (1 + applied_unc * sqrt(3) / 100)

        # --- Mapping idx → function to extract (r_in, r_ext) pair ---
		layer_accessor = Dict(
			# Core conductor
			1  => s -> begin	# co_d
				layer = s.components[1].conductor_group.layers[1]
				(layer.radius_in, layer.radius_ext)
			end,
			13 => s -> begin	# co_lay
				layer = s.components[1].conductor_group.layers[2]
				(layer.lay_ratio)
			end,

			# Inner semicon
			2  => s -> begin	# t_sc_in
				layer = s.components[1].insulator_group.layers[1]
				(layer.radius_in, layer.radius_ext)
			end,

			# Insulator
			3  => s -> begin	# t_ins
				layer = s.components[1].insulator_group.layers[2]
				(layer.radius_in, layer.radius_ext)
			end,

			# Outer semicon
			4  => s -> begin	# t_sc_out
				layer = s.components[1].insulator_group.layers[3]
				(layer.radius_in, layer.radius_ext)
			end,

			# Wire diameters
			5  => s -> begin	# sc_d
				layer = s.components[2].conductor_group.layers[1]
				(layer.radius_in, layer.radius_ext)
			end,
			14 => s -> begin	# sc_lay
				layer = s.components[2].conductor_group.layers[1]
				(layer.lay_ratio)
			end,

			# Cut (strip)
			6  => s -> begin	# t_cut
				layer = s.components[2].conductor_group.layers[2]
				(layer.radius_in, layer.radius_ext)
			end,
			7  => s -> begin	# w_cut	
				layer = s.components[2].conductor_group.layers[2]
				(layer.width)
			end,

			# Jacket
			10 => s -> begin	# t_jac
				layer = s.components[2].insulator_group.layers[1]
				(layer.radius_in, layer.radius_ext)
			end
		)
		# Ensure the index is mapped
		haskey(layer_accessor, idx) || error("No layer mapping for idx=$idx")

		for i in 1:10
			s = sample(cbs_det, distribution = :uniform)
			
			if idx == 13 || idx == 14  # lay ratio special case
				lay_ratio = layer_accessor[idx](s)
				ok = lower_sample ≤ lay_ratio ≤ upper_sample	
				@test ok || error("Sample #$i out of range for idx=$idx: lay_ratio=$lay_ratio but expected ∈ [$lower_sample, $upper_sample]")
				continue
			elseif idx == 7  # width special case
				width = layer_accessor[idx](s)
				ok = lower_sample ≤ width ≤ upper_sample	
				@test ok || error("Sample #$i out of range for idx=$idx: width=$width but expected ∈ [$lower_sample, $upper_sample]")
				continue
			else
				# Regular thickness-based cable_parameters
				r_in, r_ext = layer_accessor[idx](s)
				if idx == 1  # diameter case core conductor
					lower_bound = deterministic_value * (1 - applied_unc * sqrt(3) / 100) / 2
					upper_bound = deterministic_value * (1 + applied_unc * sqrt(3) / 100) / 2
					lo = r_in + lower_bound
					hi = r_in + upper_bound
					ok = lo ≤ r_ext ≤ hi
					@test ok || error("Sample #$i out of range for idx=$idx: r_ext=$r_ext but expected ∈ [$lo, $hi]")
				else  # thickness case
					lo = r_in + lower_sample
					hi = r_in + upper_sample
					ok = lo ≤ r_ext ≤ hi
					@test ok || error("Sample #$i out of range for idx=$idx: r_ext=$r_ext but expected ∈ [$lo, $hi]")
				end				
			end
		end
    end

	@info("Testing uncertainty propagation for all material parameters with assigned uncertainty...")
	logical_materials = Dict{Int, String}()
	for layer_str in layers_to_build
		parts = split(layer_str, ',')
		idx = parse(Int, parts[1])
		logical_materials[idx] = String(parts[5]) # Store material ID
	end

	# Calculate how many layers the Core has.
	num_core_layers = length(cbs.parts) - length(layers_to_build) + 1

	# ============================================================
	# TEST LOOP OVER LOGICAL LAYERS
	# ============================================================
	cbs_cursor = 1

	for logical_idx in sort(collect(keys(logical_materials)))
		
		mat_id = logical_materials[logical_idx]
		
		# Retrieve uncertainty for this material
		expected_uncs = get(material_uncertainties, mat_id, (rho=nothing, mu_r=nothing, eps_r=nothing))

		# Determine how many physical parts this logical layer consumes
		is_core = (logical_idx == 1)
		parts_to_check = is_core ? num_core_layers : 1
		
		# Loop through the physical parts belonging to this logical layer
		for offset in 0:(parts_to_check - 1)
			current_part_idx = cbs_cursor + offset
			
			# Determine the sub-layer index for the sample accessor
			sub_layer_idx = offset + 1 
			
			# --- DEFINE ACCESSOR FOR 'SAMPLE' (s) ---
			get_sample_mat = if is_core
				# Core: components[1] -> conductor -> layers[1..N]
				s -> s.components[1].conductor_group.layers[sub_layer_idx].material_props
			elseif logical_idx == 2
				s -> s.components[1].insulator_group.layers[1].material_props # Semicon In
			elseif logical_idx == 3
				s -> s.components[1].insulator_group.layers[2].material_props # XLPE
			elseif logical_idx == 4
				s -> s.components[1].insulator_group.layers[3].material_props # Semicon Out
			elseif logical_idx == 5
				s -> s.components[2].conductor_group.layers[1].material_props # Screen Wires
			elseif logical_idx == 6
				s -> s.components[2].conductor_group.layers[2].material_props # Screen Strip
			elseif logical_idx == 7
				s -> s.components[2].insulator_group.layers[1].material_props # Jacket
			else
				nothing
			end

			if isnothing(get_sample_mat)
				continue
			end

			# --- TEST PROPERTIES ---
			for prop in [:rho, :mu_r, :eps_r]
				applied_unc = getproperty(expected_uncs, prop)
				(isnothing(applied_unc) || applied_unc == 0.0) && continue
				@info("-------------------------------")
				@info "Checking Logical Layer $logical_idx (Part $current_part_idx) - $mat_id - $prop"

				# 1. CHECK SPECIFICATION (cbs.parts[current_part_idx])
				spec_part = cbs.parts[current_part_idx]
				val_struct = getproperty(spec_part.material, prop)
				
				nominal_val = val_struct isa Tuple ? val_struct[1] : val_struct
				is_tuple = val_struct isa Tuple{Float64, Union{Nothing, Float64}}
				@test is_tuple && (val_struct[2] ≈ applied_unc) || error("Part $current_part_idx ($prop) not (val, unc). Found: $val_struct")

				# 2. CHECK DETERMINIZED
				det_part = cbs_det.parts[current_part_idx]
				det_val = getproperty(det_part.material, prop)

				lower_ex = nominal_val * (1 - applied_unc/100)
				upper_ex = nominal_val * (1 + applied_unc/100)
				@info("Expected bounds for Part $current_part_idx ($prop): [$lower_ex, $upper_ex]")
				# Handle Tuple/Vector outputs from cbs_det
				if det_val isa Tuple{Tuple{Float64, Float64, Int64}, Nothing}
					lo = det_val[1][1]
					hi = det_val[1][2]
					@info("Extracted bounds for Part $current_part_idx ($prop): [$lo, $hi]")
					# Perform the check on the extracted bounds
					is_correct = (lo ≈ lower_ex) && (hi ≈ upper_ex)
					@test is_correct || error("Bounds mismatch Part $current_part_idx ($prop). Expected [$lower_ex, $upper_ex], got [$lo, $hi]")
				end

				# 3. CHECK SAMPLING
				lower_samp = nominal_val * (1 - applied_unc * sqrt(3) / 100)
				upper_samp = nominal_val * (1 + applied_unc * sqrt(3) / 100)

				for i in 1:10
					s = sample(cbs_det, distribution = :uniform)
					samp_mat = get_sample_mat(s)
					val = getproperty(samp_mat, prop)
					@test lower_samp <= val <= upper_samp || error("Sample out of bounds")
				end
			end
		end
		cbs_cursor += parts_to_check
	end
end

# # --- Monte Carlo on cable arameters ---
# @time res = mc(
#    cbs,
#    return_pdf = true,
#    return_samples = true,
#    dkw_eps = 0.01,
#    distribution = :gaussian,
# )
# show(res.summary, allrows = true, allcols = true)

# # --- Plots of sampled parameters ---
# p1 = Plots.histogram(
#    res.samples.μ_R,
#    bins = 80,
#    normalize = :pdf,
#    xlabel = "R",
#    ylabel = "density",
#    label = "R",
# )
# p2 = Plots.histogram(
#    res.samples.μ_L,
#    bins = 80,
#    normalize = :pdf,
#    xlabel = "L",
#    ylabel = "density",
#    label = "L",
# ) 
# display(p1)
# display(p2)

# frequency grid is deterministic
f = 10.0 .^ range(0, stop = 6, length = 10)

# Define system center point (underground at 1 m depth) and the trifoil positions
positions = [
    trifoil(
        y0 = -1.0,
        d = (system_parameters[1], system_uncertainties[1]),  # distance between cables
        phases = (
            :core   => (1, 2, 3),
            :sheath => (0, 0, 0),
            # :jacket => (0, 0, 0),
        ),
    ),
]

# Earth model (uniform, 100 Ωm with uncertainty)
earth = Earth(rho = (system_parameters[2], system_uncertainties[2]), eps_r = (system_parameters[3], system_uncertainties[3]), mu_r = (system_parameters[4], system_uncertainties[4]))

# System spec
spec = SystemBuilder("trifoil_case", cbs, positions;
	length = (system_parameters[5], system_uncertainties[5]),
	temperature = (system_parameters[6], system_uncertainties[6]),
	earth = earth,
	f = f,
)
# Determinize and sample the system
det_sys = determinize(spec)
sample_sys = sample(det_sys, distribution = :uniform)

# Define runtime options 
opts = (
	force_overwrite = true,                    # Overwrite existing files
	save_path = fullfile("lineparams_output"), # Results directory
	verbosity = 0,                             # Verbosity
);

# Define the EMT-type model with the specified formulations
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
for problem in spec
	_, p = compute!(problem, F)  #  FormulationSet as before
	push!(Ps, p)
end