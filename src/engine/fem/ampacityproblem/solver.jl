
# function make_fem_problem!(
# 	fem_formulation::AmpacityFormulation,
# 	frequency::Float64,
# 	workspace::ElectroThermalFEMWorkspace,
# )

# 	fem_formulation.problem = GetDP.Problem()
# 	define_jacobian!(fem_formulation.problem, workspace)
# 	define_integration!(fem_formulation.problem)
# 	define_material_props!(fem_formulation.problem, workspace)
# 	define_constants!(fem_formulation.problem, fem_formulation, frequency, workspace )
# 	define_domain_groups!(fem_formulation.problem, fem_formulation, workspace)
# 	define_constraint!(fem_formulation.problem, fem_formulation, workspace)
# 	define_resolution!(fem_formulation.problem, fem_formulation, workspace)

# 	make_problem!(fem_formulation.problem)

# 	index = findfirst(item -> item isa typeof(fem_formulation), workspace.formulation.analysis_type)
# 	fem_formulation.problem.filename = workspace.paths[:multiphysics_file][index]
# 	@info "Writing multiphysics problem to $(fem_formulation.problem.filename)"

# 	write_file(fem_formulation.problem)
# end

# function define_constants!(
# 	problem::GetDP.Problem,
# 	fem_formulation::AmpacityFormulation,
# 	frequency::Float64,
# 	workspace::ElectroThermalFEMWorkspace,
# )
# 	func = GetDP.Function()

# 	add_constant!(func, "Freq", frequency)
# 	add_constant!(func, "Tambient[]", workspace.temperature+273.15) # Kelvin
# 	add_constant!(func, "V_wind", workspace.wind_velocity) # m/s
# 	add_constant!(func, "h[]", "7.371 + 6.43*V_wind^0.75") # Convective coefficient [W/(m^2 K)]

# 	push!(problem.function_obj, func)
# end

# function define_domain_groups!(
# 	problem::GetDP.Problem,
# 	fem_formulation::AmpacityFormulation,
# 	workspace::ElectroThermalFEMWorkspace,
# )

# 	material_reg = Dict{Symbol, Vector{Int}}(
# 		:DomainC => Int[],
# 		:DomainCC => Int[],
# 		:DomainInf => Int[],
# 	)
# 	inds_reg = Int[]
# 	cables_reg = Dict{Int, Vector{Int}}()
# 	boundary_reg = Int[]
# 	is_magneto_thermal = fem_formulation isa MagnetoThermal

# 	for tag in keys(workspace.physical_groups)
# 		if tag > 10^8
# 			# Decode tag information
# 			surface_type, entity_num, component_num, material_group, _ =
# 				decode_physical_group_tag(tag)

# 			# Categorize regions 
# 			if surface_type == 1
# 				push!(get!(cables_reg, entity_num, Int[]), tag)
# 				if material_group == 1 && (!is_magneto_thermal || component_num == 1)
# 					push!(inds_reg, tag)
# 				end
# 			end
# 			if material_group == 1
# 				push!(material_reg[:DomainC], tag)
# 			elseif material_group == 2
# 				push!(material_reg[:DomainCC], tag)
# 			end

# 			surface_type == 3 && push!(material_reg[:DomainInf], tag)

# 		else
# 			decode_boundary_tag(tag)[1] == 2 && push!(boundary_reg, tag)
# 		end
# 	end
# 	inds_reg = sort(inds_reg)
# 	material_reg[:DomainC] = sort(material_reg[:DomainC])
# 	material_reg[:DomainCC] = sort(material_reg[:DomainCC])

# 	# Create and configure groups
# 	group = GetDP.Group()

# 	# Add common domains
# 	add!(
# 		group,
# 		"DomainInf",
# 		material_reg[:DomainInf],
# 		"Region",
# 		comment = "Domain transformation to infinity",
# 	)

# 	for (key, tag) in enumerate(inds_reg)
# 		add!(group, "Con_$key", [tag], "Region";
# 			comment = "$(create_physical_group_name(workspace, tag))")
# 	end

# 	add!(group, "Conductors", inds_reg, "Region")

# 	# Add standard FEM domains
# 	domain_configs = [
# 		("DomainC", Int[], "All conductor materials"),
# 		("DomainCC", Int[], "All non-conductor materials"),
# 	]

# 	for (name, regions, comment) in domain_configs
# 		add!(group, name, regions, "Region"; comment = comment)
# 	end

# 	for tag in material_reg[:DomainC]
# 		add!(group, "DomainC", [tag], "Region";
# 			operation = "+=",
# 			comment = "$(create_physical_group_name(workspace, tag))")
# 	end

# 	for tag in material_reg[:DomainCC]
# 		add!(group, "DomainCC", [tag], "Region";
# 			operation = "+=",
# 			comment = "$(create_physical_group_name(workspace, tag))")
# 	end


# 	if fem_formulation isa MagnetoThermal

# 		add!(group, "Domain_Mag", ["DomainCC", "DomainC"], "Region")
# 		add!(group, "Sur_Dirichlet_Mag", boundary_reg, "Region")
# 		add!(group, "Sur_Dirichlet_The", boundary_reg, "Region")
# 		add!(group, "Sur_Convection_Thermal", [], "Region")
# 	end

# 	problem.group = group
# end

# function define_constraint!(
# 	problem::GetDP.Problem,
# 	fem_formulation::MagnetoThermal,
# 	workspace::ElectroThermalFEMWorkspace,
# )
# 	constraint = GetDP.Constraint()

# 	# MagneticVectorPotential_2D
# 	mvp = assign!(constraint, "MagneticVectorPotential_2D")
# 	case!(mvp, "Sur_Dirichlet_Mag", value = "0.0")

# 	# Voltage_2D (placeholder)
# 	voltage = assign!(constraint, "Voltage_2D")
# 	case!(voltage, "")

# 	# Current_2D
# 	current = assign!(constraint, "Current_2D")
# 	for (idx, curr) in enumerate(workspace.energizations)
# 		case!(current, "Con_$idx", value = "Complex[$(real(curr)), $(imag(curr))]")
# 	end

# 	temp = assign!(constraint, "DirichletTemp")
# 	case!(temp, "Sur_Dirichlet_The", value = "Tambient[]") # Ambient temperature

# 	problem.constraint = constraint

# end


# function define_resolution!(
# 	problem::GetDP.Problem,
# 	formulation::MagnetoThermal,
# 	workspace::ElectroThermalFEMWorkspace,
# )

# 	resolution_name = formulation.resolution_name

# 	# Create a new Problem instance
# 	functionspace = FunctionSpace()

# 	# FunctionSpace section
# 	fs1 = add!(functionspace, "Hcurl_a_Mag_2D", nothing, nothing, Type = "Form1P")
# 	add_basis_function!(
# 		functionspace,
# 		"se",
# 		"ae",
# 		"BF_PerpendicularEdge";
# 		Support = "Domain_Mag",
# 		Entity = "NodesOf[ All ]",
# 	)

# 	add_constraint!(functionspace, "ae", "NodesOf", "MagneticVectorPotential_2D")

# 	fs3 = add!(functionspace, "Hregion_u_Mag_2D", nothing, nothing, Type = "Form1P")
# 	add_basis_function!(
# 		functionspace,
# 		"sr",
# 		"ur",
# 		"BF_RegionZ";
# 		Support = "DomainC",
# 		Entity = "DomainC",
# 	)
# 	add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef = "ur")
# 	add_global_quantity!(functionspace, "I", "AssociatedWith"; NameOfCoef = "ur")
# 	add_constraint!(functionspace, "U", "Region", "Voltage_2D")
# 	add_constraint!(functionspace, "I", "Region", "Current_2D")

# 	fs1 = add!(functionspace, "Hgrad_Thermal", nothing, nothing, Type = "Form0")
# 	add_basis_function!(
# 		functionspace,
# 		"sn",
# 		"t",
# 		"BF_Node";
# 		Support = "Domain_Mag",
# 		Entity = "NodesOf[ All ]",
# 	)

# 	add_constraint!(functionspace, "t", "NodesOf", "DirichletTemp")

# 	problem.functionspace = functionspace

# 	# Define Formulation
# 	formulation = GetDP.Formulation()

# 	form = add!(formulation, "Darwin_a_2D", "FemEquation")
# 	add_quantity!(form, "a", Type = "Local", NameOfSpace = "Hcurl_a_Mag_2D")
# 	add_quantity!(form, "ur", Type = "Local", NameOfSpace = "Hregion_u_Mag_2D")
# 	add_quantity!(form, "T", Type = "Local", NameOfSpace = "Hgrad_Thermal")
# 	add_quantity!(form, "I", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [I]")
# 	add_quantity!(form, "U", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [U]")

# 	eq = add_equation!(form)

# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ nu[] * Dof{d a} , {d a} ]",
# 		In = "Domain_Mag",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDof [ sigma[{T}] * Dof{a} , {a} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ sigma[{T}] * Dof{ur}, {a} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDof [ sigma[{T}] * Dof{a} , {ur} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ sigma[{T}] * Dof{ur}, {ur}]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDtDof [ epsilon[] * Dof{a} , {a}]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 		comment = " Darwin approximation term",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDof[ epsilon[] * Dof{ur}, {a} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDtDof [ epsilon[] * Dof{a} , {ur}]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(
# 		eq,
# 		"Galerkin",
# 		"DtDof[ epsilon[] * Dof{ur}, {ur} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)
# 	add!(eq, "GlobalTerm", "[ Dof{I} , {U} ]", In = "Conductors") #DomainActive

# 	form = add!(formulation, "ThermalSta", "FemEquation")
# 	add_quantity!(form, "T", Type = "Local", NameOfSpace = "Hgrad_Thermal")
# 	add_quantity!(form, "a", Type = "Local", NameOfSpace = "Hcurl_a_Mag_2D")
# 	add_quantity!(form, "ur", Type = "Local", NameOfSpace = "Hregion_u_Mag_2D")

# 	eq = add_equation!(form)

# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ k[] * Dof{d T} , {d T} ]",
# 		In = "Domain_Mag",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)

# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ -0.5*sigma[{T}] * <a>[ SquNorm[Dt[{a}]+{ur}] ], {T} ]",
# 		In = "DomainC",
# 		Jacobian = "Vol",
# 		Integration = "I1",
# 	)

# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[ h[] * Dof{T} , {T} ]",
# 		In = "Sur_Convection_Thermal",
# 		Jacobian = "Sur",
# 		Integration = "I1",
# 		comment = " Convection boundary condition",
# 	)

# 	add!(
# 		eq,
# 		"Galerkin",
# 		"[-h[] * Tambient[] , {T} ]",
# 		In = "Sur_Convection_Thermal",
# 		Jacobian = "Sur",
# 		Integration = "I1",
# 	)
# 	# Add the formulation to the problem
# 	problem.formulation = formulation

# 	# Define Resolution
# 	resolution = Resolution()

# 	sys_mag = SystemItem("Sys_Mag", "Darwin_a_2D"; Type="Complex", Frequency="Freq")
# 	sys_the = SystemItem("Sys_The", "ThermalSta")

# 	# Add a resolution
# 	output_dir = joinpath("results", lowercase(resolution_name))
# 	output_dir = replace(output_dir, "\\" => "/")     # for compatibility with Windows paths

# 	# Construct the final Operation vector
# 	add!(resolution, resolution_name, [sys_mag, sys_the],
# 		Operation = [
# 			"CreateDir[\"$(output_dir)\"]",
# 			"InitSolution[Sys_Mag]",
# 			"InitSolution[Sys_The]",
# 			"Generate[Sys_Mag]",
# 			"Solve[Sys_Mag]",
# 			"Generate[Sys_The]",
# 			"Solve[Sys_The]",
# 			"SaveSolution[Sys_Mag]",
# 			"SaveSolution[Sys_The]",
# 			"PostOperation[LineParams]",
# 		]
# 	)
# 	# Add the resolution to the problem
# 	problem.resolution = resolution

# 	# PostProcessing section
# 	postprocessing = PostProcessing()

# 	pp = add!(postprocessing, "Darwin_a_2D", "Darwin_a_2D")

# 	q = add_post_quantity_term!(pp, "a")
# 	add!(q, "Term", "{a}"; In = "Domain_Mag", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "az")
# 	add!(q, "Term", "CompZ[{a}]"; In = "Domain_Mag", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "b")
# 	add!(q, "Term", "{d a}"; In = "Domain_Mag", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "bm")
# 	add!(q, "Term", "Norm[{d a}]"; In = "Domain_Mag", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "j")
# 	add!(q, "Term", "-sigma[]*(Dt[{a}]+{ur})"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "jz")
# 	add!(q, "Term", "CompZ[-sigma[{T}]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "jm")
# 	add!(q, "Term", "Norm[-sigma[{T}]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "d")
# 	add!(q, "Term", "epsilon[] * Dt[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "dz")
# 	add!(q, "Term", "CompZ[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "dm")
# 	add!(q, "Term", "Norm[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
# 	q = add_post_quantity_term!(pp, "rhoj2")
# 	add!(q, "Term", "0.5*sigma[{T}]*SquNorm[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")

# 	q = add_post_quantity_term!(pp, "U")
# 	add!(q, "Term", "{U}"; In = "DomainC")
# 	q = add_post_quantity_term!(pp, "I")
# 	add!(q, "Term", "{I}"; In = "DomainC")
# 	q = add_post_quantity_term!(pp, "Z")
# 	add!(q, "Term", "-{U}"; In = "DomainC")


# 	pp_the = add!(postprocessing, "ThermalSta", "ThermalSta")

# 	q = add_quantity_term!(pp_the, "T")
# 	add!(q, "Local", "{T}"; In = "Domain_Mag", Jacobian = "Vol")

# 	q = add_quantity_term!(pp_the, "TinC")
# 	add!(q, "Local", "{T}-273.15"; In = "Domain_Mag", Jacobian = "Vol")

# 	q = add_quantity_term!(pp_the, "q")
# 	add!(q, "Local", "-k[]*{d T}"; In = "Domain_Mag", Jacobian = "Vol")

# 	problem.postprocessing = postprocessing

# 	# PostOperation section
# 	postoperation = PostOperation()

# 	# Add post-operation items
# 	po1 = add!(postoperation, "Field_Maps", "ThermalSta")
# 	op1 = add_operation!(po1)
# 	add_operation!(op1, "Print[ TinC, OnElementsOf Domain_Mag, Name \"T [°C] around cable\", File StrCat[ \"$(joinpath(output_dir,"TinC"))\", \".pos\" ] ];")
# 	add_operation!(op1, "Print[ q , OnElementsOf Domain_Mag, Name \"heat flux [W/m²] around cable\", File StrCat[ \"$(joinpath(output_dir,"q"))\", \".pos\" ] ];")

# 	po2 = add!(postoperation, "LineParams", "ThermalSta")
# 	op2 = add_operation!(po2)
# 	add_operation!(op2, "Print[ TinC, OnElementsOf Domain_Mag, Name \"T [°C] around cable\", File StrCat[ \"$(joinpath(output_dir,"TinC"))\", \".pos\" ] ];")
# 	add_operation!(op2, "Print[ q , OnElementsOf Domain_Mag, Name \"heat flux [W/m²] around cable\", File StrCat[ \"$(joinpath(output_dir,"q"))\", \".pos\" ] ];")
# 	# add_operation!(op2, "Print[ Z, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Z.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];")


# 	# Add the post-operation to the problem
# 	problem.postoperation = postoperation

# end

# function run_getdp(workspace::FEMWorkspace, fem_formulation::AmpacityFormulation)
# 	# Initialize Gmsh if not already initialized
# 	if gmsh.is_initialized() == 0
# 		gmsh.initialize()
# 	end

# 	# Number of iterations (from the original function)
# 	n_phases =
# 		sum([length(c.design_data.components) for c in workspace.core.system.cables])

# 	# Flag to track if all solves are successful
# 	all_success = true

# 	# Map verbosity to Gmsh/GetDP level
# 	gmsh_verbosity = map_verbosity_to_gmsh(workspace.opts.verbosity)
# 	gmsh.option.set_number("General.Verbosity", gmsh_verbosity)


# 	getdp_verbosity = map_verbosity_to_getdp(workspace.opts.verbosity)

# 	# Construct solver command with -setnumber active_ind i
# 	solve_cmd = "$(workspace.opts.getdp_executable) $(fem_formulation.problem.filename) -msh $(workspace.paths[:mesh_file]) -solve $(fem_formulation.resolution_name) -v2 -verbose $(getdp_verbosity)"

# 	# Log the current solve attempt
# 	@info "Solving multiphysics... (Resolution = $(fem_formulation.resolution_name))"

# 	# Attempt to run the solver
# 	try
# 		gmsh.onelab.run("GetDP", solve_cmd)

# 		if workspace.opts.plot_field_maps
# 			@info "Building field maps multiphysics... (Resolution = $(fem_formulation.resolution_name))"

# 			post_cmd = "$(workspace.opts.getdp_executable) $(fem_formulation.problem.filename) -msh $(workspace.paths[:mesh_file]) -pos Field_Maps -v2 -verbose $(getdp_verbosity)"

# 			gmsh.onelab.run("GetDP", post_cmd)
# 		end

# 		@info "Solve successful!"
# 	catch e
# 		# Log the error and update the success flag
# 		@error "Solver failed: $e"
# 		all_success = false
# 		# Continue to the next iteration even if this one fails
# 	end

# 	# Return true only if all solves were successful
# 	return all_success
# end

# using LinearAlgebra: BLAS, BlasFloat

# # function run_solver!(::Val{:AmpacityProblem}, workspace::ElectroThermalFEMWorkspace)

# # 	n_frequencies = workspace.n_frequencies

# # 	for (k, frequency) in enumerate(workspace.freq)
# # 		@info "Solving frequency $k/$n_frequencies: $frequency Hz"

# # 		# Fill Z,Y (original ordering) for this slice
# # 		_do_run_solver!(k, workspace)

# # 		# Archive if requested
# # 		if workspace.opts.keep_run_files
# # 			archive_frequency_results(workspace, frequency)
# # 		end
# # 	end

# # end


# # function _do_run_solver!(freq_idx::Int,
# # 	workspace::ElectroThermalFEMWorkspace) # Z::Array{ComplexF64, 3}, Y::Array{ComplexF64, 3})

# # 	# Get formulation from workspace
# # 	formulation = workspace.formulation
# # 	# Z, Y = workspace.Z, workspace.Y
# # 	frequency = workspace.freq[freq_idx]

# # 	# Build and solve both formulations
# # 	for fem_formulation in formulation.analysis_type
# # 		@debug "Processing $(fem_formulation.resolution_name) formulation"

# # 		make_fem_problem!(fem_formulation, frequency, workspace)

# # 		if !run_getdp(workspace, fem_formulation)
# # 			Base.error("$(fem_formulation.resolution_name) solver failed")
# # 		end
# # 	end

# # end

