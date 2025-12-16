

@inline function get_output_filename(::AbstractImpedanceFormulation, workspace::FEMWorkspace)
    return workspace.core.paths[:impedance_file]
end

@inline function get_output_filename(::AbstractAdmittanceFormulation, workspace::FEMWorkspace)
    return workspace.core.paths[:admittance_file]
end

@inline function get_output_filename(::AmpacityFormulation, workspace::FEMWorkspace)
    return workspace.core.paths[:analysis_file]
end

struct DefineJacobian end
struct DefineIntegration end
struct DefineMaterialProps end
struct DefineConstants end
struct DefineDomainGroups end
struct DefineConstraint end
struct DefineResolution end

@inline function (f::Union{AbstractImpedanceFormulation, AbstractAdmittanceFormulation})(
    frequency::Float64,
    workspace::FEMWorkspace,
	active_cond::Int,
)
    # Create and build the problem
    getdp_problem = GetDP.Problem()

    DefineJacobian()(getdp_problem, workspace)
    DefineIntegration()(getdp_problem)
    DefineMaterialProps()(getdp_problem, workspace)
    DefineConstants()(getdp_problem, frequency)
    DefineDomainGroups()(getdp_problem, f, workspace, active_cond)
    DefineConstraint()(getdp_problem, f)
    DefineResolution()(getdp_problem, f)

    make_problem!(getdp_problem)
	
	# Set the filename on the problem object based on formulation type
    getdp_problem.filename = get_output_filename(f, workspace)

	write_file(getdp_problem)

	return getdp_problem
end

@inline function (f::AmpacityFormulation)(
    frequency::Float64,
    workspace::FEMWorkspace,
)
    # Create and build the problem
    getdp_problem = GetDP.Problem()
    
    DefineJacobian()(getdp_problem, workspace)
    DefineIntegration()(getdp_problem)
    DefineMaterialProps()(getdp_problem, workspace)
    DefineConstants()(getdp_problem, frequency, workspace)
    DefineDomainGroups()(getdp_problem, f, workspace)
    DefineConstraint()(getdp_problem, workspace)
    DefineResolution()(getdp_problem, f)

    make_problem!(getdp_problem)
	
	# Set the filename on the problem object based on formulation type
    getdp_problem.filename = get_output_filename(f, workspace)

	write_file(getdp_problem)

	return getdp_problem
end


@inline function (f::DefineJacobian)(problem::GetDP.Problem, workspace::FEMWorkspace)
	# Initialize Jacobian
	jac = GetDP.Jacobian()

	Rint = workspace.core.formulation.domain_radius
	Rext = workspace.core.formulation.domain_radius_inf

	# Add Vol Jacobian
	vol = add!(jac, "Vol")
	add!(vol;
		Region = "DomainInf",
		Jacobian = GetDP.VolSphShell(
			Rint = Rint,
			Rext = Rext,
			center_X = 0.0,
			center_Y = 0.0,
			center_Z = 0.0,
		),
	)
	add!(vol; Region = "All", Jacobian = "Vol")

	# Add Sur Jacobian
	sur = add!(jac, "Sur")
	add!(sur;
		Region = "All",
		Jacobian = "Sur",
	)

	# Add Jacobian to problem
	problem.jacobian = jac
end

@inline function (f::DefineIntegration)(problem::GetDP.Problem)
	# Initialize Integration
	integ = GetDP.Integration()
	i1 = add!(integ, "I1")
	case = add!(i1)
	geo_case = GetDP.add_nested_case!(case; type = "Gauss")
	add!(geo_case; GeoElement = "Point", NumberOfPoints = 1)
	add!(geo_case; GeoElement = "Line", NumberOfPoints = 4)
	add!(geo_case; GeoElement = "Triangle", NumberOfPoints = 4)
	add!(geo_case; GeoElement = "Quadrangle", NumberOfPoints = 4)
	problem.integration = integ

end

@inline function (f::DefineMaterialProps)(problem::GetDP.Problem, workspace::FEMWorkspace)
	# Create material properties function
	func = GetDP.Function()

	for (tag, mat) in workspace.core.physical_groups
		if tag > 10^8
			# Add material properties for this region
			add_comment!(
				func,
				"Material properties for region $(tag): $(create_physical_group_name(workspace, tag))",
				false,
			)
			add_space!(func)
			add!(func, "nu", expression = 1 / (mat.mu_r * μ₀), region = [tag])
			add!(
				func,
				"sigma",
				expression = isinf(mat.rho) ? 0.0 : 1 / mat.rho,
				region = [tag],
			)
			add!(func, "epsilon", expression = mat.eps_r * ε₀, region = [tag])
			add!(func, "k", expression = mat.kappa, region = [tag])
		end
	end

	push!(problem.function_obj, func)
end


@inline function (f::DefineConstants)(
	problem::GetDP.Problem,
	frequency::Float64,
)
	func = GetDP.Function()

	add_constant!(func, "Freq", frequency)
	add_constant!(func, "UnitAmplitude", 1.0)
	
	push!(problem.function_obj, func)
end

@inline function (f::DefineConstants)(
	problem::GetDP.Problem,
	frequency::Float64,
	workspace::FEMWorkspace,
)
	func = GetDP.Function()

	add_constant!(func, "Freq", frequency)
	add_constant!(func, "Tambient[]", workspace.core.temp+273.15) # Kelvin
	add_constant!(func, "V_wind", workspace.wind_velocity) # m/s
	add_constant!(func, "h[]", "7.371 + 6.43*V_wind^0.75") # Convective coefficient [W/(m^2 K)]

	push!(problem.function_obj, func)
end

@inline function (f::DefineDomainGroups)(
	problem::GetDP.Problem,
	fem_formulation::Union{AbstractImpedanceFormulation, AbstractAdmittanceFormulation},
	workspace::FEMWorkspace,
	active_cond::Int,
)

	material_reg = Dict{Symbol, Vector{Int}}(
		:DomainC => Int[],
		:DomainCC => Int[],
		:DomainInf => Int[],
	)
	inds_reg = Int[]
	cables_reg = Dict{Int, Vector{Int}}()
	boundary_reg = Int[]
	add_raw_code!(problem,
		"""
		active_con = $active_cond;
		""")
	for tag in keys(workspace.core.physical_groups)
		if tag > 10^8
			# Decode tag information
			surface_type, entity_num, component_num, material_group, _ =
				decode_physical_group_tag(tag)

			# Categorize regions 
			if surface_type == 1
				push!(get!(cables_reg, entity_num, Int[]), tag)
				if material_group == 1
					push!(inds_reg, tag)
				end
			end
			if material_group == 1
				push!(material_reg[:DomainC], tag)
			elseif material_group == 2
				push!(material_reg[:DomainCC], tag)
			end

			surface_type == 3 && push!(material_reg[:DomainInf], tag)

		else
			decode_boundary_tag(tag)[1] == 2 && push!(boundary_reg, tag)
		end
	end
	inds_reg = sort(inds_reg)
	material_reg[:DomainC] = sort(material_reg[:DomainC])
	material_reg[:DomainCC] = sort(material_reg[:DomainCC])

	# Create and configure groups
	group = GetDP.Group()

	# Add common domains
	add!(
		group,
		"DomainInf",
		material_reg[:DomainInf],
		"Region",
		comment = "Domain transformation to infinity",
	)

	for (key, tag) in enumerate(inds_reg)
		add!(group, "Con_$key", [tag], "Region";
			comment = "$(create_physical_group_name(workspace, tag))")
	end

	add!(group, "Conductors", inds_reg, "Region")

	# Add standard FEM domains
	domain_configs = [
		("DomainC", Int[], "All conductor materials"),
		("DomainCC", Int[], "All non-conductor materials"),
		("DomainActive", ["Con~{active_con}"], "Sources"),
		(
			"DomainInactive",
			["Conductors - Con~{active_con}"],
			"Conductors set to zero energization",
		),
	]

	for (name, regions, comment) in domain_configs
		add!(group, name, regions, "Region"; comment = comment)
	end

	for tag in material_reg[:DomainC]
		add!(group, "DomainC", [tag], "Region";
			operation = "+=",
			comment = "$(create_physical_group_name(workspace, tag))")
	end

	for tag in material_reg[:DomainCC]
		add!(group, "DomainCC", [tag], "Region";
			operation = "+=",
			comment = "$(create_physical_group_name(workspace, tag))")
	end

	if fem_formulation isa AbstractAdmittanceFormulation
		add!(group, "Domain_Ele", ["DomainCC", "DomainC"], "Region")
		add!(group, "Sur_Dirichlet_Ele", boundary_reg, "Region")

	else
		# Add domain groups
		add!(group, "Domain_Mag", ["DomainCC", "DomainC"], "Region")
		add!(group, "Sur_Dirichlet_Mag", boundary_reg, "Region")
	end

	problem.group = group
end

@inline function (f::DefineDomainGroups)(
	problem::GetDP.Problem,
	fem_formulation::MagnetoThermal,
	workspace::FEMWorkspace,
)

	material_reg = Dict{Symbol, Vector{Int}}(
		:DomainC => Int[],
		:DomainCC => Int[],
		:DomainInf => Int[],
	)
	inds_reg = Int[]
	cables_reg = Dict{Int, Vector{Int}}()
	boundary_reg = Int[]
	is_magneto_thermal = fem_formulation isa MagnetoThermal

	for tag in keys(workspace.core.physical_groups)
		if tag > 10^8
			# Decode tag information
			surface_type, entity_num, component_num, material_group, _ =
				decode_physical_group_tag(tag)

			# Categorize regions 
			if surface_type == 1
				push!(get!(cables_reg, entity_num, Int[]), tag)
				if material_group == 1 && (!is_magneto_thermal || component_num == 1)
					push!(inds_reg, tag)
				end
			end
			if material_group == 1
				push!(material_reg[:DomainC], tag)
			elseif material_group == 2
				push!(material_reg[:DomainCC], tag)
			end

			surface_type == 3 && push!(material_reg[:DomainInf], tag)

		else
			decode_boundary_tag(tag)[1] == 2 && push!(boundary_reg, tag)
		end
	end
	inds_reg = sort(inds_reg)
	material_reg[:DomainC] = sort(material_reg[:DomainC])
	material_reg[:DomainCC] = sort(material_reg[:DomainCC])

	# Create and configure groups
	group = GetDP.Group()

	# Add common domains
	add!(
		group,
		"DomainInf",
		material_reg[:DomainInf],
		"Region",
		comment = "Domain transformation to infinity",
	)

	for (key, tag) in enumerate(inds_reg)
		add!(group, "Con_$key", [tag], "Region";
			comment = "$(create_physical_group_name(workspace, tag))")
	end

	add!(group, "Conductors", inds_reg, "Region")

	# Add standard FEM domains
	domain_configs = [
		("DomainC", Int[], "All conductor materials"),
		("DomainCC", Int[], "All non-conductor materials"),
	]

	for (name, regions, comment) in domain_configs
		add!(group, name, regions, "Region"; comment = comment)
	end

	for tag in material_reg[:DomainC]
		add!(group, "DomainC", [tag], "Region";
			operation = "+=",
			comment = "$(create_physical_group_name(workspace, tag))")
	end

	for tag in material_reg[:DomainCC]
		add!(group, "DomainCC", [tag], "Region";
			operation = "+=",
			comment = "$(create_physical_group_name(workspace, tag))")
	end


	add!(group, "Domain_Mag", ["DomainCC", "DomainC"], "Region")
	add!(group, "Sur_Dirichlet_Mag", boundary_reg, "Region")
	add!(group, "Sur_Dirichlet_The", boundary_reg, "Region")
	add!(group, "Sur_Convection_Thermal", [], "Region")
	
	problem.group = group
end


@inline function (f::DefineConstraint)(
	problem::GetDP.Problem,
	fem_formulation::Union{AbstractImpedanceFormulation, AbstractAdmittanceFormulation},
)
	constraint = GetDP.Constraint()

	# num_cores = workspace.core.system.num_cables

	if fem_formulation isa AbstractAdmittanceFormulation
		# ScalarPotential_2D
		esp = assign!(constraint, "ScalarPotential_2D")
		case!(esp, "DomainInactive", value = "0.0")
		case!(esp, "Con~{active_con}", value = "UnitAmplitude")
		case!(esp, "Sur_Dirichlet_Ele", value = "0.0")

		charge = assign!(constraint, "Charge_2D")
	else
		# MagneticVectorPotential_2D
		mvp = assign!(constraint, "MagneticVectorPotential_2D")
		case!(mvp, "Sur_Dirichlet_Mag", value = "0.0")

		# Voltage_2D (placeholder)
		voltage = assign!(constraint, "Voltage_2D")
		case!(voltage, "")

		# Current_2D
		current = assign!(constraint, "Current_2D")

		case!(current, "DomainInactive", value = "0.0")
		case!(current, "Con~{active_con}", value = "UnitAmplitude")
	end

	problem.constraint = constraint

end

@inline function (f::DefineConstraint)(
	problem::GetDP.Problem,
	workspace::FEMWorkspace,
)
	constraint = GetDP.Constraint()

	# MagneticVectorPotential_2D
	mvp = assign!(constraint, "MagneticVectorPotential_2D")
	case!(mvp, "Sur_Dirichlet_Mag", value = "0.0")

	# Voltage_2D (placeholder)
	voltage = assign!(constraint, "Voltage_2D")
	case!(voltage, "")

	# Current_2D
	current = assign!(constraint, "Current_2D")
	for (idx, curr) in enumerate(workspace.energizations)
		case!(current, "Con_$idx", value = "Complex[$(real(curr)), $(imag(curr))]")
	end

	temp = assign!(constraint, "DirichletTemp")
	case!(temp, "Sur_Dirichlet_The", value = "Tambient[]") # Ambient temperature

	problem.constraint = constraint

end


@inline function (f::DefineResolution)(
	problem::GetDP.Problem,
	formulation_type::Electrodynamics,
)

	# FunctionSpace section
	functionspace = FunctionSpace()
	fs1 = add!(functionspace, "Hgrad_v_Ele", nothing, nothing, Type = "Form0")
	add_basis_function!(
		functionspace,
		"sn",
		"vn",
		"BF_Node";
		Support = "Domain_Ele",
		Entity = "NodesOf[ All, Not Conductors ]",
	)
	add_basis_function!(
		functionspace,
		"sf",
		"vf",
		"BF_GroupOfNodes";
		Support = "Domain_Ele",
		Entity = "GroupsOfNodesOf[ Conductors ]",
	)
	add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef = "vf")
	add_global_quantity!(functionspace, "Q", "AssociatedWith"; NameOfCoef = "vf")
	add_constraint!(functionspace, "U", "Region", "ScalarPotential_2D")
	add_constraint!(functionspace, "Q", "Region", "Charge_2D")
	add_constraint!(functionspace, "vn", "NodesOf", "ScalarPotential_2D")

	problem.functionspace = functionspace

	# Formulation section
	formulation = Formulation()
	form = add!(formulation, "Electrodynamics_v", "FemEquation")
	add_quantity!(form, "v", Type = "Local", NameOfSpace = "Hgrad_v_Ele")
	add_quantity!(form, "U", Type = "Global", NameOfSpace = "Hgrad_v_Ele [U]")
	add_quantity!(form, "Q", Type = "Global", NameOfSpace = "Hgrad_v_Ele [Q]")

	eq = add_equation!(form)
	add!(
		eq,
		"Galerkin",
		"[ sigma[] * Dof{d v} , {d v} ]",
		In = "Domain_Ele",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof[ epsilon[] * Dof{d v} , {d v} ]",
		In = "DomainCC",
		Jacobian = "Vol",
		Integration = "I1",
	) #CHECKME
	add!(eq, "GlobalTerm", "[ Dof{Q} , {U} ]", In = "Conductors")

	problem.formulation = formulation

	# Resolution section
	resolution_name = get_resolution_name(formulation_type)
	output_dir = joinpath("results", lowercase(resolution_name))
	output_dir = replace(output_dir, "\\" => "/") # for compatibility with Windows paths
	resolution = Resolution()
	sys_ele = SystemItem("Sys_Ele", "Electrodynamics_v"; 
		Type="Complex", 
		Frequency="Freq"
	)
	add!(resolution, resolution_name, [sys_ele],
		Operation = [
			"CreateDir[\"$(output_dir)\"]",
			"Generate[Sys_Ele]",
			"Solve[Sys_Ele]",
			"SaveSolution[Sys_Ele]",
			"PostOperation[LineParams]",
		])

	problem.resolution = resolution

	# PostProcessing section
	postprocessing = PostProcessing()
	pp = add!(postprocessing, "EleDyn_v", "Electrodynamics_v")

	# Add field maps quantities
	for (name, expr, options) in [
		("v", "{v}", Dict()),
		("e", "-{d v}", Dict()),
		("em", "Norm[-{d v}]", Dict()),
		("d", "-epsilon[] * {d v}", Dict()),
		("dm", "Norm[-epsilon[] * {d v}]", Dict()),
		("j", "-sigma[] * {d v}", Dict()),
		("jm", "Norm[-sigma[] * {d v}]", Dict()),
	]
		q = add_post_quantity_term!(pp, name)
		add!(q, "Term", expr; In = "Domain_Ele", Jacobian = "Vol", options...)
	end

	# Add jtot (combination of j and d)
	q = add_post_quantity_term!(pp, "jtot")
	add!(
		q,
		"Term",
		"-sigma[] * {d v}";
		Type = "Global",
		In = "Domain_Ele",
		Jacobian = "Vol",
	)
	add!(
		q,
		"Term",
		"-epsilon[] * Dt[{d v}]";
		Type = "Global",
		In = "Domain_Ele",
		Jacobian = "Vol",
	)

	q = add_post_quantity_term!(pp, "U")
	add!(q, "Term", "{U}"; In = "Domain_Ele")

	q = add_post_quantity_term!(pp, "Q")
	add!(q, "Term", "{Q}"; In = "Domain_Ele")

	q = add_post_quantity_term!(pp, "Y")
	add!(q, "Term", "-{Q}"; In = "Domain_Ele")

	problem.postprocessing = postprocessing

	# PostOperation section
	postoperation = PostOperation()

	# Field_Maps
	po1 = add!(postoperation, "Field_Maps", "EleDyn_v")
	op1 = add_operation!(po1)
	add_operation!(
		op1,
		"Print[ v, OnElementsOf Domain_Ele, File StrCat[ \"$(joinpath(output_dir,"v_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ em, OnElementsOf Domain_Ele, Name \"|E| [V/m]\", File StrCat[ \"$(joinpath(output_dir,"em_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ dm, OnElementsOf Domain_Ele, Name \"|D| [A/m²]\", File StrCat[ \"$(joinpath(output_dir,"dm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ e, OnElementsOf Domain_Ele, Name \"E [V/m]\", File StrCat[ \"$(joinpath(output_dir,"e_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)

	# LineParams
	po2 = add!(postoperation, "LineParams", "EleDyn_v")
	op2 = add_operation!(po2)
	add_operation!(
		op2,
		"Print[ Y, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Y.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];",
	)

	problem.postoperation = postoperation

end

@inline function (f::DefineResolution)(
	problem::GetDP.Problem,
	formulation_type::Darwin,
)

	# Create a new Problem instance
	functionspace = FunctionSpace()

	# FunctionSpace section
	fs1 = add!(functionspace, "Hcurl_a_Mag_2D", nothing, nothing, Type = "Form1P")
	add_basis_function!(
		functionspace,
		"se",
		"ae",
		"BF_PerpendicularEdge";
		Support = "Domain_Mag",
		Entity = "NodesOf[ All ]",
	)

	add_constraint!(functionspace, "ae", "NodesOf", "MagneticVectorPotential_2D")

	fs3 = add!(functionspace, "Hregion_u_Mag_2D", nothing, nothing, Type = "Form1P")
	add_basis_function!(
		functionspace,
		"sr",
		"ur",
		"BF_RegionZ";
		Support = "DomainC",
		Entity = "DomainC",
	)
	add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef = "ur")
	add_global_quantity!(functionspace, "I", "AssociatedWith"; NameOfCoef = "ur")
	add_constraint!(functionspace, "U", "Region", "Voltage_2D")
	add_constraint!(functionspace, "I", "Region", "Current_2D")

	problem.functionspace = functionspace

	# Define Formulation
	formulation = GetDP.Formulation()

	form = add!(formulation, "Darwin_a_2D", "FemEquation")
	add_quantity!(form, "a", Type = "Local", NameOfSpace = "Hcurl_a_Mag_2D")
	add_quantity!(form, "ur", Type = "Local", NameOfSpace = "Hregion_u_Mag_2D")
	add_quantity!(form, "I", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [I]")
	add_quantity!(form, "U", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [U]")

	eq = add_equation!(form)

	add!(
		eq,
		"Galerkin",
		"[ nu[] * Dof{d a} , {d a} ]",
		In = "Domain_Mag",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof [ sigma[] * Dof{a} , {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"[ sigma[] * Dof{ur}, {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof [ sigma[] * Dof{a} , {ur} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"[ sigma[] * Dof{ur}, {ur}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDtDof [ epsilon[] * Dof{a} , {a}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
		comment = " Darwin approximation term",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof[ epsilon[] * Dof{ur}, {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDtDof [ epsilon[] * Dof{a} , {ur}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof[ epsilon[] * Dof{ur}, {ur} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(eq, "GlobalTerm", "[ Dof{I} , {U} ]", In = "Conductors") #DomainActive

	# Add the formulation to the problem
	problem.formulation = formulation

	# Define Resolution
	resolution = Resolution()
	sys_mag = SystemItem("Sys_Mag", "Darwin_a_2D"; 
		Type="Complex", 
		Frequency="Freq"
	)
	# Add a resolution
	resolution_name = get_resolution_name(formulation_type)
	output_dir = joinpath("results", lowercase(resolution_name))
	output_dir = replace(output_dir, "\\" => "/")     # for compatibility with Windows paths
	add!(resolution, resolution_name, [sys_mag],
		Operation = [
			"CreateDir[\"$(output_dir)\"]",
			"InitSolution[Sys_Mag]",
			"Generate[Sys_Mag]",
			"Solve[Sys_Mag]",
			"SaveSolution[Sys_Mag]",
			"PostOperation[LineParams]",
		])

	# Add the resolution to the problem
	problem.resolution = resolution

	# PostProcessing section
	postprocessing = PostProcessing()

	pp = add!(postprocessing, "Darwin_a_2D", "Darwin_a_2D")
	q = add_post_quantity_term!(pp, "a")
	add!(q, "Term", "{a}"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "az")
	add!(q, "Term", "CompZ[{a}]"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "b")
	add!(q, "Term", "{d a}"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "bm")
	add!(q, "Term", "Norm[{d a}]"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "j")
	add!(q, "Term", "-sigma[]*(Dt[{a}]+{ur})"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "jz")
	add!(q, "Term", "CompZ[-sigma[]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "jm")
	add!(q, "Term", "Norm[-sigma[]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "d")
	add!(q, "Term", "epsilon[] * Dt[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "dz")
	add!(q, "Term", "CompZ[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "dm")
	add!(q, "Term", "Norm[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "rhoj2")
	add!(q, "Term", "0.5*sigma[]*SquNorm[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")

	q = add_post_quantity_term!(pp, "U")
	add!(q, "Term", "{U}"; In = "DomainC")
	q = add_post_quantity_term!(pp, "I")
	add!(q, "Term", "{I}"; In = "DomainC")
	q = add_post_quantity_term!(pp, "Z")
	add!(q, "Term", "-{U}"; In = "DomainC")

	problem.postprocessing = postprocessing

	# PostOperation section
	postoperation = PostOperation()

	# Add post-operation items
	po1 = add!(postoperation, "Field_Maps", "Darwin_a_2D")
	op1 = add_operation!(po1)

	add_operation!(
		op1,
		"Print[ az, OnElementsOf Domain_Mag, Smoothing 1, Name \"flux lines: Az [T m]\", File StrCat[ \"$(joinpath(output_dir,"az_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ b, OnElementsOf Domain_Mag, Smoothing 1, Name \"B [T]\", File StrCat[ \"$(joinpath(output_dir,"b_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ bm, OnElementsOf Domain_Mag, Smoothing 1, Name \"|B| [T]\", File StrCat[ \"$(joinpath(output_dir,"bm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ jz, OnElementsOf Region[{DomainC}], Smoothing 1, Name \"jz [A/m²] Conducting domain\", File StrCat[ \"$(joinpath(output_dir,"jz_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ rhoj2, OnElementsOf Region[{DomainC}], Smoothing 1, Name \"Power density\", File StrCat[ \"$(joinpath(output_dir,"rhoj2_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ jm, OnElementsOf DomainC, Smoothing 1, Name \"|j| [A/m²] Conducting domain\", File StrCat[ \"$(joinpath(output_dir,"jm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)
	add_operation!(
		op1,
		"Print[ dm, OnElementsOf DomainC, Smoothing 1, Name \"|D| [A/m²]\", File StrCat[ \"$(joinpath(output_dir,"dm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];",
	)

	po2 = add!(postoperation, "LineParams", "Darwin_a_2D")
	op2 = add_operation!(po2)
	add_operation!(
		op2,
		"Print[ Z, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Z.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];",
	)

	# Add the post-operation to the problem
	problem.postoperation = postoperation

end
@inline function (f::DefineResolution)(
	problem::GetDP.Problem,
	formulation_type::MagnetoThermal,
)
	# Create a new Problem instance
	functionspace = FunctionSpace()

	# FunctionSpace section
	fs1 = add!(functionspace, "Hcurl_a_Mag_2D", nothing, nothing, Type = "Form1P")
	add_basis_function!(
		functionspace,
		"se",
		"ae",
		"BF_PerpendicularEdge";
		Support = "Domain_Mag",
		Entity = "NodesOf[ All ]",
	)

	add_constraint!(functionspace, "ae", "NodesOf", "MagneticVectorPotential_2D")

	fs3 = add!(functionspace, "Hregion_u_Mag_2D", nothing, nothing, Type = "Form1P")
	add_basis_function!(
		functionspace,
		"sr",
		"ur",
		"BF_RegionZ";
		Support = "DomainC",
		Entity = "DomainC",
	)
	add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef = "ur")
	add_global_quantity!(functionspace, "I", "AssociatedWith"; NameOfCoef = "ur")
	add_constraint!(functionspace, "U", "Region", "Voltage_2D")
	add_constraint!(functionspace, "I", "Region", "Current_2D")

	fs1 = add!(functionspace, "Hgrad_Thermal", nothing, nothing, Type = "Form0")
	add_basis_function!(
		functionspace,
		"sn",
		"t",
		"BF_Node";
		Support = "Domain_Mag",
		Entity = "NodesOf[ All ]",
	)

	add_constraint!(functionspace, "t", "NodesOf", "DirichletTemp")

	problem.functionspace = functionspace

	# Define Formulation
	formulation = GetDP.Formulation()

	form = add!(formulation, "Darwin_a_2D", "FemEquation")
	add_quantity!(form, "a", Type = "Local", NameOfSpace = "Hcurl_a_Mag_2D")
	add_quantity!(form, "ur", Type = "Local", NameOfSpace = "Hregion_u_Mag_2D")
	add_quantity!(form, "T", Type = "Local", NameOfSpace = "Hgrad_Thermal")
	add_quantity!(form, "I", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [I]")
	add_quantity!(form, "U", Type = "Global", NameOfSpace = "Hregion_u_Mag_2D [U]")

	eq = add_equation!(form)

	add!(
		eq,
		"Galerkin",
		"[ nu[] * Dof{d a} , {d a} ]",
		In = "Domain_Mag",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof [ sigma[{T}] * Dof{a} , {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"[ sigma[{T}] * Dof{ur}, {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof [ sigma[{T}] * Dof{a} , {ur} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"[ sigma[{T}] * Dof{ur}, {ur}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDtDof [ epsilon[] * Dof{a} , {a}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
		comment = " Darwin approximation term",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof[ epsilon[] * Dof{ur}, {a} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDtDof [ epsilon[] * Dof{a} , {ur}]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(
		eq,
		"Galerkin",
		"DtDof[ epsilon[] * Dof{ur}, {ur} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)
	add!(eq, "GlobalTerm", "[ Dof{I} , {U} ]", In = "Conductors") #DomainActive

	form = add!(formulation, "ThermalSta", "FemEquation")
	add_quantity!(form, "T", Type = "Local", NameOfSpace = "Hgrad_Thermal")
	add_quantity!(form, "a", Type = "Local", NameOfSpace = "Hcurl_a_Mag_2D")
	add_quantity!(form, "ur", Type = "Local", NameOfSpace = "Hregion_u_Mag_2D")

	eq = add_equation!(form)

	add!(
		eq,
		"Galerkin",
		"[ k[] * Dof{d T} , {d T} ]",
		In = "Domain_Mag",
		Jacobian = "Vol",
		Integration = "I1",
	)

	add!(
		eq,
		"Galerkin",
		"[ -0.5*sigma[{T}] * <a>[ SquNorm[Dt[{a}]+{ur}] ], {T} ]",
		In = "DomainC",
		Jacobian = "Vol",
		Integration = "I1",
	)

	add!(
		eq,
		"Galerkin",
		"[ h[] * Dof{T} , {T} ]",
		In = "Sur_Convection_Thermal",
		Jacobian = "Sur",
		Integration = "I1",
		comment = " Convection boundary condition",
	)

	add!(
		eq,
		"Galerkin",
		"[-h[] * Tambient[] , {T} ]",
		In = "Sur_Convection_Thermal",
		Jacobian = "Sur",
		Integration = "I1",
	)
	# Add the formulation to the problem
	problem.formulation = formulation

	# Define Resolution
	resolution = Resolution()

	sys_mag = SystemItem("Sys_Mag", "Darwin_a_2D"; Type="Complex", Frequency="Freq")
	sys_the = SystemItem("Sys_The", "ThermalSta")

	# Add a resolution
	resolution_name = get_resolution_name(formulation_type)
	output_dir = joinpath("results", lowercase(resolution_name))
	output_dir = replace(output_dir, "\\" => "/")     # for compatibility with Windows paths

	# Construct the final Operation vector
	add!(resolution, resolution_name, [sys_mag, sys_the],
		Operation = [
			"CreateDir[\"$(output_dir)\"]",
			"InitSolution[Sys_Mag]",
			"InitSolution[Sys_The]",
			"Generate[Sys_Mag]",
			"Solve[Sys_Mag]",
			"Generate[Sys_The]",
			"Solve[Sys_The]",
			"SaveSolution[Sys_Mag]",
			"SaveSolution[Sys_The]",
			"PostOperation[LineParams]",
		]
	)
	# Add the resolution to the problem
	problem.resolution = resolution

	# PostProcessing section
	postprocessing = PostProcessing()

	pp = add!(postprocessing, "Darwin_a_2D", "Darwin_a_2D")

	q = add_post_quantity_term!(pp, "a")
	add!(q, "Term", "{a}"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "az")
	add!(q, "Term", "CompZ[{a}]"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "b")
	add!(q, "Term", "{d a}"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "bm")
	add!(q, "Term", "Norm[{d a}]"; In = "Domain_Mag", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "j")
	add!(q, "Term", "-sigma[]*(Dt[{a}]+{ur})"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "jz")
	add!(q, "Term", "CompZ[-sigma[{T}]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "jm")
	add!(q, "Term", "Norm[-sigma[{T}]*(Dt[{a}]+{ur})]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "d")
	add!(q, "Term", "epsilon[] * Dt[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "dz")
	add!(q, "Term", "CompZ[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "dm")
	add!(q, "Term", "Norm[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In = "DomainC", Jacobian = "Vol")
	q = add_post_quantity_term!(pp, "rhoj2")
	add!(q, "Term", "0.5*sigma[{T}]*SquNorm[Dt[{a}]+{ur}]"; In = "DomainC", Jacobian = "Vol")

	q = add_post_quantity_term!(pp, "U")
	add!(q, "Term", "{U}"; In = "DomainC")
	q = add_post_quantity_term!(pp, "I")
	add!(q, "Term", "{I}"; In = "DomainC")
	q = add_post_quantity_term!(pp, "Z")
	add!(q, "Term", "-{U}"; In = "DomainC")


	pp_the = add!(postprocessing, "ThermalSta", "ThermalSta")

	q = add_quantity_term!(pp_the, "T")
	add!(q, "Local", "{T}"; In = "Domain_Mag", Jacobian = "Vol")

	q = add_quantity_term!(pp_the, "TinC")
	add!(q, "Local", "{T}-273.15"; In = "Domain_Mag", Jacobian = "Vol")

	q = add_quantity_term!(pp_the, "q")
	add!(q, "Local", "-k[]*{d T}"; In = "Domain_Mag", Jacobian = "Vol")

	problem.postprocessing = postprocessing

	# PostOperation section
	postoperation = PostOperation()

	# Add post-operation items
	po1 = add!(postoperation, "Field_Maps", "ThermalSta")
	op1 = add_operation!(po1)
	add_operation!(op1, "Print[ TinC, OnElementsOf Domain_Mag, Name \"T [°C] around cable\", File StrCat[ \"$(joinpath(output_dir,"TinC"))\", \".pos\" ] ];")
	add_operation!(op1, "Print[ q , OnElementsOf Domain_Mag, Name \"heat flux [W/m²] around cable\", File StrCat[ \"$(joinpath(output_dir,"q"))\", \".pos\" ] ];")

	po2 = add!(postoperation, "LineParams", "ThermalSta")
	op2 = add_operation!(po2)
	add_operation!(op2, "Print[ TinC, OnElementsOf Domain_Mag, Name \"T [°C] around cable\", File StrCat[ \"$(joinpath(output_dir,"TinC"))\", \".pos\" ] ];")
	add_operation!(op2, "Print[ q , OnElementsOf Domain_Mag, Name \"heat flux [W/m²] around cable\", File StrCat[ \"$(joinpath(output_dir,"q"))\", \".pos\" ] ];")
	# add_operation!(op2, "Print[ Z, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Z.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];")


	# Add the post-operation to the problem
	problem.postoperation = postoperation

end

function run_getdp(workspace::FEMWorkspace, frequency::Float64,
	fem_formulation::Union{AbstractImpedanceFormulation, AbstractAdmittanceFormulation}, active_con::Int)

	resolution_name = get_resolution_name(fem_formulation)
    @debug "Processing $(resolution_name) formulation"
    getdp_problem = fem_formulation(frequency, workspace, active_con)
    getdp_verbosity = map_verbosity_to_getdp(workspace.core.opts.verbosity)
    solve_cmd = "$(workspace.core.opts.getdp_executable) $(getdp_problem.filename) -msh $(workspace.core.paths[:mesh_file]) -solve $(resolution_name) -v2 -verbose $(getdp_verbosity)"

    @info "Solving for source conductor $active_con... (Resolution = $(resolution_name))"

    # Attempt to run the solver
    success = true
    try
        gmsh.onelab.run("GetDP", solve_cmd)

        if workspace.core.opts.plot_field_maps
            @info "Building field maps for source conductor $active_con... (Resolution = $(resolution_name))"

            post_cmd = "$(workspace.core.opts.getdp_executable) $(getdp_problem.filename) -msh $(workspace.core.paths[:mesh_file]) -pos Field_Maps -v2 -verbose $(getdp_verbosity)"

            gmsh.onelab.run("GetDP", post_cmd)
        end

        @info "Solve successful for source conductor $(active_con)!"
    catch e
        # Log the error and update the success flag
        @error "Solver failed for source conductor $active_con: $e"
        success = false
    end
    
	return success

end
function run_getdp(workspace::FEMWorkspace, problem::GetDP.Problem, fem_formulation::AmpacityFormulation)
	# Initialize Gmsh if not already initialized
	if gmsh.is_initialized() == 0
		gmsh.initialize()
	end

	# Number of iterations (from the original function)
	n_phases =
		sum([length(c.design_data.components) for c in workspace.core.system.cables])

	# Flag to track if all solves are successful
	all_success = true

	# Map verbosity to Gmsh/GetDP level
	gmsh_verbosity = map_verbosity_to_gmsh(workspace.core.opts.verbosity)
	gmsh.option.set_number("General.Verbosity", gmsh_verbosity)


	getdp_verbosity = map_verbosity_to_getdp(workspace.core.opts.verbosity)
	resolution_name = get_resolution_name(fem_formulation)

	# Construct solver command with -setnumber active_ind i
	solve_cmd = "$(workspace.core.opts.getdp_executable) $(problem.filename) -msh $(workspace.core.paths[:mesh_file]) -solve $(resolution_name) -v2 -verbose $(getdp_verbosity)"

	# Log the current solve attempt
	@info "Solving multiphysics... (Resolution = $(resolution_name))"

	# Attempt to run the solver
	try
		gmsh.onelab.run("GetDP", solve_cmd)

		if workspace.core.opts.plot_field_maps
			@info "Building field maps multiphysics... (Resolution = $(resolution_name))"

			post_cmd = "$(workspace.core.opts.getdp_executable) $(problem.filename) -msh $(workspace.core.paths[:mesh_file]) -pos Field_Maps -v2 -verbose $(getdp_verbosity)"

			gmsh.onelab.run("GetDP", post_cmd)
		end

		@info "Solve successful!"
	catch e
		# Log the error and update the success flag
		@error "Solver failed: $e"
		all_success = false
		# Continue to the next iteration even if this one fails
	end

	# Return true only if all solves were successful
	return all_success
end


using LinearAlgebra: BLAS, BlasFloat

"""
$(TYPEDSIGNATURES)

Main function to run the FEM simulation workflow for a cable system.

# Arguments

- `cable_system`: Cable system to simulate.
- `formulation`: Problem definition parameters.
- `solver`: Solver parameters.
- `frequency`: Simulation frequency \\[Hz\\]. Default: 50.0.

# Returns

- A [`FEMWorkspace`](@ref) instance with the simulation results.

# Examples

```julia
# Run a FEM simulation
workspace = $(FUNCTIONNAME)(cable_system, formulation, solver)
```
"""

function compute!(problem::LineParametersProblem,
    formulation::FEMFormulation,
    workspace::Union{FEMWorkspace, Nothing} = nothing)

    opts = formulation.options

    # Initialize workspace
    workspace = init_workspace(problem, formulation, workspace)

    # Meshing phase: make_mesh! decides if it needs to run.
    # It returns true if the process should stop (e.g., mesh_only=true).
    mesh_only_flag = make_mesh!(workspace)
        
    # Only proceed with solver if not mesh_only
    if !mesh_only_flag
        @info "Starting FEM solver"
        
        # Extract necessary variables from workspace
        n_phases = workspace.core.n_phases
        n_frequencies = workspace.core.n_frequencies
        phase_map = workspace.core.phase_map

        # --- Index plan (once) ---
        perm = reorder_indices(phase_map)    # encounter-ordered: first of each phase, then tails, then zeros
        map_r = phase_map[perm]              # reordered map (constant across k)

        # --- Outputs: size decided by kron_map (here: map_r after merge_bundles! zeros tails) ---
        # Probe the keep-size once using a scratch (no heavy cost)
        _probe = Matrix{ComplexF64}(I, n_phases, n_phases)
        _, reduced_map = merge_bundles!(copy(_probe), map_r)
        n_keep = count(!=(0), reduced_map)

        Zr = zeros(ComplexF64, n_keep, n_keep, n_frequencies)
        Yr = zeros(ComplexF64, n_keep, n_keep, n_frequencies)

        # --- Scratch buffers (reused every k) ---
        Zbuf = Matrix{ComplexF64}(undef, n_phases, n_phases)   # reordered + merged target
        Ybuf = Matrix{ComplexF64}(undef, n_phases, n_phases)
        Pf = Matrix{ComplexF64}(undef, n_phases, n_phases)     # potentials (for Y path)

        # Tiny gather helper: reorder src[:,:,k] into dest without temp allocs
        @inline function _reorder_into!(dest::StridedMatrix{ComplexF64},
            src::Array{ComplexF64, 3},
            perm::Vector{Int}, k::Int)
            n = length(perm)
            @inbounds for j in 1:n, i in 1:n
                dest[i, j] = src[perm[i], perm[j], k]
            end
            dest
        end

        # --- Big loop over frequencies ---
        for (k, frequency) in enumerate(workspace.core.freq)
            @info "Solving frequency $k/$n_frequencies: $frequency Hz"

            # Build and solve both formulations
            for fem_formulation_item in formulation.analysis_type

				# Initialize Gmsh if not already initialized
                if gmsh.is_initialized() == 0
                    gmsh.initialize()
                end

                # Number of iterations
                n_phases_inner = sum([length(c.design_data.components) for c in workspace.core.system.cables])

                # Flag to track if all solves are successful
                all_success = true

                # Map verbosity to Gmsh/GetDP level
                gmsh_verbosity = map_verbosity_to_gmsh(workspace.core.opts.verbosity)
                gmsh.option.set_number("General.Verbosity", gmsh_verbosity)

                # Loop over each active_ind from 1 to n_phases
                for i in 1:n_phases_inner
                    # Construct solver command with -setnumber active_ind i
					all_success = run_getdp(workspace, frequency, fem_formulation_item, i)
                end

                # Check if solve was successful
                if !all_success
                    Base.error("$(get_resolution_name(fem_formulation_item)) solver failed")
                end
            end

            # Extract results into preallocated arrays
            workspace.Z[:, :, k] = read_results_file(formulation.analysis_type[1], workspace)
            workspace.Y[:, :, k] = read_results_file(formulation.analysis_type[2], workspace)

            # REORDER → Z
            _reorder_into!(Zbuf, workspace.Z, perm, k)

            # MERGE bundles (in-place on Zbuf) and get reduced map (tails → 0)
            Zm, reduced_map = merge_bundles!(Zbuf, map_r)

            # KRON on Z
            Zred = kronify(Zm, reduced_map)
            symtrans!(Zred)
            formulation.options.ideal_transposition || line_transpose!(Zred)
            @inbounds Zr[:, :, k] .= Zred

            # Y path goes via potentials: Pf = inv(Y/(jω))
            w = 2π * frequency
            # REORDER → Y
            _reorder_into!(Ybuf, workspace.Y, perm, k)

            # Pf = inv(Ybuf / (jω)) without extra temps
            @inbounds @views begin
                Pf .= Ybuf
                Pf ./= (1im * w)
            end
            Pf .= inv(Pf)

            # MERGE bundles for Pf (same reduced_map semantics)
            Pfm, reduced_map = merge_bundles!(Pf, map_r)

            # KRON on Pf, then invert back to Y
            Pr = kronify(Pfm, reduced_map)
            Yrk = (1im * w) * inv(Pr)
            symtrans!(Yrk)
            formulation.options.ideal_transposition || line_transpose!(Yrk)
            @inbounds Yr[:, :, k] .= Yrk

            # Archive if requested
            if workspace.core.opts.keep_run_files
                archive_frequency_results(workspace, frequency)
            end
        end

        ZY = LineParameters(Zr, Yr, workspace.core.freq)
        @info "FEM computation completed successfully"
    end

    return workspace, ZY
end
"""
$(TYPEDSIGNATURES)

Main function to run the FEM simulation workflow for a cable system.

# Arguments
- `problem`: Ampacity problem definition.
- `formulation`: Formulation parameters.
- `workspace`: (Optional) Pre-initialized [`FEMWorkspace`](@ref) instance.

# Returns
- A [`FEMWorkspace`](@ref) instance with the simulation results.

# Examples
```julia
# Run a FEM simulation
workspace = $(FUNCTIONNAME)(problem, formulation)
```
"""
function compute!(problem::AmpacityProblem,
	formulation::FEMFormulation,
	workspace::Union{FEMWorkspace, Nothing} = nothing)

	opts = formulation.options
	# Initialize workspace
	workspace = init_workspace(problem, formulation, workspace)

	# Meshing phase: make_mesh! decides if it needs to run.
	# It returns true if the process should stop (e.g., mesh_only=true).
	if make_mesh!(workspace)
		return workspace, nothing
	end

	fem_formulation = workspace.core.formulation.analysis_type

	# Solving phase - always runs unless mesh_only
	@info "Starting FEM solver"
	for freq in workspace.core.freq
		@debug "Processing $(workspace.core.formulation.resolution_name) formulation"

		getdp_problem = fem_formulation(freq, workspace)

		if !run_getdp(workspace, getdp_problem, fem_formulation)
			Base.error("$(fem_formulation.resolution_name) solver failed")
		end
		# Archive if requested
		if workspace.core.opts.keep_run_files
			archive_frequency_results(workspace, freq)
		end
	end

	@info "FEM computation completed successfully"
	return workspace
end

