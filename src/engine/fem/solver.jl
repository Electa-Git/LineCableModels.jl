
function make_fem_problem!(fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation}, frequency::Float64, workspace::FEMWorkspace)

    fem_formulation.problem = GetDP.Problem()
    define_jacobian!(fem_formulation.problem, workspace)
    define_integration!(fem_formulation.problem)
    define_material_props!(fem_formulation.problem, workspace)
    define_constants!(fem_formulation.problem, fem_formulation, frequency)
    define_domain_groups!(fem_formulation.problem, fem_formulation, workspace)
    define_constraint!(fem_formulation.problem, fem_formulation, workspace)
    define_resolution!(fem_formulation.problem, fem_formulation, workspace)

    make_problem!(fem_formulation.problem)
    fem_formulation.problem.filename = fem_formulation isa AbstractImpedanceFormulation ? workspace.paths[:impedance_file] : workspace.paths[:admittance_file]
    write_file(fem_formulation.problem)
end

function define_jacobian!(problem::GetDP.Problem, workspace::FEMWorkspace)
    # Initialize Jacobian
    jac = Jacobian()

    Rint = workspace.formulation.domain_radius
    Rext = workspace.formulation.domain_radius_inf

    # Add Vol Jacobian
    vol = add!(jac, "Vol")
    add!(vol;
        Region="DomainInf",
        Jacobian=VolSphShell(
            Rint=Rint,
            Rext=Rext,
            center_X=0.0,
            center_Y=0.0,
            center_Z=0.0
        )
    )
    add!(vol; Region="All", Jacobian="Vol")

    # Add Sur Jacobian
    sur = add!(jac, "Sur")
    add!(sur;
        Region="All",
        Jacobian="Sur",
    )

    # Add Jacobian to problem
    problem.jacobian = jac
end

function define_integration!(problem::GetDP.Problem)
    # Initialize Integration
    integ = Integration()
    i1 = add!(integ, "I1")
    case = add!(i1)
    geo_case = add_nested_case!(case; type="Gauss")
    add!(geo_case; GeoElement="Point", NumberOfPoints=1)
    add!(geo_case; GeoElement="Line", NumberOfPoints=4)
    add!(geo_case; GeoElement="Triangle", NumberOfPoints=4)
    add!(geo_case; GeoElement="Quadrangle", NumberOfPoints=4)
    problem.integration = integ

end

function define_material_props!(problem::GetDP.Problem, workspace::FEMWorkspace)
    # Create material properties function
    func = GetDP.Function()

    for (tag, mat) in workspace.physical_groups
        if tag > 10^8
            # Add material properties for this region
            add_comment!(func, "Material properties for region $(tag): $(create_physical_group_name(workspace, tag))", false)
            add_space!(func)
            add!(func, "nu", expression=1 / (mat.mu_r * μ₀), region=[tag])
            add!(func, "sigma", expression=isinf(mat.rho) ? 0.0 : 1 / mat.rho, region=[tag])
            add!(func, "epsilon", expression=mat.eps_r * ε₀, region=[tag])
        end
    end

    push!(problem.function_obj, func)
end

function define_constants!(problem::GetDP.Problem, fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation}, frequency::Float64)
    func = GetDP.Function()

    add_constant!(func, "Freq", frequency)
    add_constant!(func, "UnitAmplitude", 1.0)
    push!(problem.function_obj, func)
end

function define_domain_groups!(problem::GetDP.Problem, fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation}, workspace::FEMWorkspace)

    material_reg = Dict{Symbol,Vector{Int}}(
        :DomainC => Int[],
        :DomainCC => Int[],
        :DomainInf => Int[]
    )
    inds_reg = Int[]
    cables_reg = Dict{Int,Vector{Int}}()
    boundary_reg = Int[]
    add_raw_code!(problem,
        """
        DefineConstant[
        active_con = {1, Choices{1,9999}, Name "Input/Active conductor", Visible 1}];
        """)
    for tag in keys(workspace.physical_groups)
        if tag > 10^8
            # Decode tag information
            surface_type, entity_num, component_num, material_group, _ = decode_physical_group_tag(tag)

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
    add!(group, "DomainInf", material_reg[:DomainInf], "Region", comment="Domain transformation to infinity")

    for (key, tag) in enumerate(inds_reg)
        add!(group, "Con_$key", [tag], "Region";
            comment="$(create_physical_group_name(workspace, tag))")
    end

    add!(group, "Conductors", inds_reg, "Region")

    # Add standard FEM domains
    domain_configs = [
        ("DomainC", Int[], "All conductor materials"),
        ("DomainCC", Int[], "All non-conductor materials"),
        ("DomainActive", ["Con~{active_con}"], "Sources"),
        ("DomainInactive", ["Conductors - Con~{active_con}"], "Conductors set to zero energization")
    ]

    for (name, regions, comment) in domain_configs
        add!(group, name, regions, "Region"; comment=comment)
    end

    for tag in material_reg[:DomainC]
        add!(group, "DomainC", [tag], "Region";
            operation="+=",
            comment="$(create_physical_group_name(workspace, tag))")
    end

    for tag in material_reg[:DomainCC]
        add!(group, "DomainCC", [tag], "Region";
            operation="+=",
            comment="$(create_physical_group_name(workspace, tag))")
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

function define_constraint!(problem::GetDP.Problem, fem_formulation::Union{AbstractImpedanceFormulation,AbstractAdmittanceFormulation}, workspace::FEMWorkspace)
    constraint = GetDP.Constraint()

    # num_cores = workspace.problem_def.system.num_cables

    if fem_formulation isa AbstractAdmittanceFormulation
        # ScalarPotential_2D
        esp = assign!(constraint, "ScalarPotential_2D")
        case!(esp, "DomainInactive", value="0.0")
        case!(esp, "Con~{active_con}", value="UnitAmplitude")
        case!(esp, "Sur_Dirichlet_Ele", value="0.0")

        charge = assign!(constraint, "Charge_2D")
    else
        # MagneticVectorPotential_2D
        mvp = assign!(constraint, "MagneticVectorPotential_2D")
        case!(mvp, "Sur_Dirichlet_Mag", value="0.0")

        # Voltage_2D (placeholder)
        voltage = assign!(constraint, "Voltage_2D")
        case!(voltage, "")

        # Current_2D
        current = assign!(constraint, "Current_2D")

        case!(current, "DomainInactive", value="0.0")
        case!(current, "Con~{active_con}", value="UnitAmplitude")
    end

    problem.constraint = constraint

end

function define_resolution!(problem::GetDP.Problem, formulation::Electrodynamics, workspace::FEMWorkspace)
    resolution_name = formulation.resolution_name
    num_sources = workspace.problem_def.system.num_cables

    # FunctionSpace section
    functionspace = FunctionSpace()
    fs1 = add!(functionspace, "Hgrad_v_Ele", nothing, nothing, Type="Form0")
    add_basis_function!(functionspace, "sn", "vn", "BF_Node"; Support="Domain_Ele", Entity="NodesOf[ All, Not Conductors ]")
    add_basis_function!(functionspace, "sf", "vf", "BF_GroupOfNodes"; Support="Domain_Ele", Entity="GroupsOfNodesOf[ Conductors ]")
    add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef="vf")
    add_global_quantity!(functionspace, "Q", "AssociatedWith"; NameOfCoef="vf")
    add_constraint!(functionspace, "U", "Region", "ScalarPotential_2D")
    add_constraint!(functionspace, "Q", "Region", "Charge_2D")
    add_constraint!(functionspace, "vn", "NodesOf", "ScalarPotential_2D")

    problem.functionspace = functionspace

    # Formulation section
    formulation = Formulation()
    form = add!(formulation, "Electrodynamics_v", "FemEquation")
    add_quantity!(form, "v", Type="Local", NameOfSpace="Hgrad_v_Ele")
    add_quantity!(form, "U", Type="Global", NameOfSpace="Hgrad_v_Ele [U]")
    add_quantity!(form, "Q", Type="Global", NameOfSpace="Hgrad_v_Ele [Q]")

    eq = add_equation!(form)
    add!(eq, "Galerkin", "[ sigma[] * Dof{d v} , {d v} ]", In="Domain_Ele", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDof[ epsilon[] * Dof{d v} , {d v} ]", In="DomainCC", Jacobian="Vol", Integration="I1") #CHECKME
    add!(eq, "GlobalTerm", "[ Dof{Q} , {U} ]", In="Conductors")

    problem.formulation = formulation

    # Resolution section
    output_dir = joinpath("results", lowercase(resolution_name))
    output_dir = replace(output_dir, "\\" => "/") # for compatibility with Windows paths
    resolution = Resolution()
    add!(resolution, resolution_name, "Sys_Ele",
        NameOfFormulation="Electrodynamics_v",
        Type="Complex",
        Frequency="Freq",
        Operation=[
            "CreateDir[\"$(output_dir)\"]",
            "Generate[Sys_Ele]",
            "Solve[Sys_Ele]",
            "SaveSolution[Sys_Ele]",
            "PostOperation[LineParams]"
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
        ("jm", "Norm[-sigma[] * {d v}]", Dict())
    ]
        q = add!(pp, name)
        add!(q, "Term", expr; In="Domain_Ele", Jacobian="Vol", options...)
    end

    # Add jtot (combination of j and d)
    q = add!(pp, "jtot")
    add!(q, "Term", "-sigma[] * {d v}"; Type="Global", In="Domain_Ele", Jacobian="Vol")
    add!(q, "Term", "-epsilon[] * Dt[{d v}]"; Type="Global", In="Domain_Ele", Jacobian="Vol")

    q = add!(pp, "U")
    add!(q, "Term", "{U}"; In="Domain_Ele")

    q = add!(pp, "Q")
    add!(q, "Term", "{Q}"; In="Domain_Ele")

    q = add!(pp, "Y")
    add!(q, "Term", "-{Q}"; In="Domain_Ele")

    problem.postprocessing = postprocessing

    # PostOperation section
    postoperation = PostOperation()

    # Field_Maps
    po1 = add!(postoperation, "Field_Maps", "EleDyn_v")
    op1 = add_operation!(po1)
    add_operation!(op1, "Print[ v, OnElementsOf Domain_Ele, File StrCat[ \"$(joinpath(output_dir,"v_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ em, OnElementsOf Domain_Ele, Name \"|E| [V/m]\", File StrCat[ \"$(joinpath(output_dir,"em_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ dm, OnElementsOf Domain_Ele, Name \"|D| [A/m²]\", File StrCat[ \"$(joinpath(output_dir,"dm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ e, OnElementsOf Domain_Ele, Name \"E [V/m]\", File StrCat[ \"$(joinpath(output_dir,"e_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")

    # LineParams
    po2 = add!(postoperation, "LineParams", "EleDyn_v")
    op2 = add_operation!(po2)
    add_operation!(op2, "Print[ Y, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Y.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];")

    problem.postoperation = postoperation

end

function define_resolution!(problem::GetDP.Problem, formulation::Darwin, workspace::FEMWorkspace)

    resolution_name = formulation.resolution_name

    # Create a new Problem instance
    functionspace = FunctionSpace()

    # FunctionSpace section
    fs1 = add!(functionspace, "Hcurl_a_Mag_2D", nothing, nothing, Type="Form1P")
    add_basis_function!(functionspace, "se", "ae", "BF_PerpendicularEdge"; Support="Domain_Mag", Entity="NodesOf[ All ]")

    add_constraint!(functionspace, "ae", "NodesOf", "MagneticVectorPotential_2D")

    fs3 = add!(functionspace, "Hregion_u_Mag_2D", nothing, nothing, Type="Form1P")
    add_basis_function!(functionspace, "sr", "ur", "BF_RegionZ"; Support="DomainC", Entity="DomainC")
    add_global_quantity!(functionspace, "U", "AliasOf"; NameOfCoef="ur")
    add_global_quantity!(functionspace, "I", "AssociatedWith"; NameOfCoef="ur")
    add_constraint!(functionspace, "U", "Region", "Voltage_2D")
    add_constraint!(functionspace, "I", "Region", "Current_2D")

    problem.functionspace = functionspace

    # Define Formulation
    formulation = GetDP.Formulation()

    form = add!(formulation, "Darwin_a_2D", "FemEquation")
    add_quantity!(form, "a", Type="Local", NameOfSpace="Hcurl_a_Mag_2D")
    add_quantity!(form, "ur", Type="Local", NameOfSpace="Hregion_u_Mag_2D")
    add_quantity!(form, "I", Type="Global", NameOfSpace="Hregion_u_Mag_2D [I]")
    add_quantity!(form, "U", Type="Global", NameOfSpace="Hregion_u_Mag_2D [U]")

    eq = add_equation!(form)

    add!(eq, "Galerkin", "[ nu[] * Dof{d a} , {d a} ]", In="Domain_Mag", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDof [ sigma[] * Dof{a} , {a} ]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "[ sigma[] * Dof{ur}, {a} ]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDof [ sigma[] * Dof{a} , {ur} ]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "[ sigma[] * Dof{ur}, {ur}]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDtDof [ epsilon[] * Dof{a} , {a}]", In="DomainC", Jacobian="Vol", Integration="I1", comment=" Darwin approximation term")
    add!(eq, "Galerkin", "DtDof[ epsilon[] * Dof{ur}, {a} ]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDtDof [ epsilon[] * Dof{a} , {ur}]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "Galerkin", "DtDof[ epsilon[] * Dof{ur}, {ur} ]", In="DomainC", Jacobian="Vol", Integration="I1")
    add!(eq, "GlobalTerm", "[ Dof{I} , {U} ]", In="Conductors") #DomainActive

    # Add the formulation to the problem
    problem.formulation = formulation

    # Define Resolution
    resolution = Resolution()

    # Add a resolution
    output_dir = joinpath("results", lowercase(resolution_name))
    output_dir = replace(output_dir, "\\" => "/")     # for compatibility with Windows paths
    add!(resolution, resolution_name, "Sys_Mag",
        NameOfFormulation="Darwin_a_2D",
        Type="Complex", Frequency="Freq",
        Operation=[
            "CreateDir[\"$(output_dir)\"]",
            "InitSolution[Sys_Mag]",
            "Generate[Sys_Mag]",
            "Solve[Sys_Mag]",
            "SaveSolution[Sys_Mag]",
            "PostOperation[LineParams]"
        ])

    # Add the resolution to the problem
    problem.resolution = resolution

    # PostProcessing section
    postprocessing = PostProcessing()

    pp = add!(postprocessing, "Darwin_a_2D", "Darwin_a_2D")
    q = add!(pp, "a")
    add!(q, "Term", "{a}"; In="Domain_Mag", Jacobian="Vol")
    q = add!(pp, "az")
    add!(q, "Term", "CompZ[{a}]"; In="Domain_Mag", Jacobian="Vol")
    q = add!(pp, "b")
    add!(q, "Term", "{d a}"; In="Domain_Mag", Jacobian="Vol")
    q = add!(pp, "bm")
    add!(q, "Term", "Norm[{d a}]"; In="Domain_Mag", Jacobian="Vol")
    q = add!(pp, "j")
    add!(q, "Term", "-sigma[]*(Dt[{a}]+{ur})"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "jz")
    add!(q, "Term", "CompZ[-sigma[]*(Dt[{a}]+{ur})]"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "jm")
    add!(q, "Term", "Norm[-sigma[]*(Dt[{a}]+{ur})]"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "d")
    add!(q, "Term", "epsilon[] * Dt[Dt[{a}]+{ur}]"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "dz")
    add!(q, "Term", "CompZ[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "dm")
    add!(q, "Term", "Norm[epsilon[] * Dt[Dt[{a}]+{ur}]]"; In="DomainC", Jacobian="Vol")
    q = add!(pp, "rhoj2")
    add!(q, "Term", "0.5*sigma[]*SquNorm[Dt[{a}]+{ur}]"; In="DomainC", Jacobian="Vol")

    q = add!(pp, "U")
    add!(q, "Term", "{U}"; In="DomainC")
    q = add!(pp, "I")
    add!(q, "Term", "{I}"; In="DomainC")
    q = add!(pp, "Z")
    add!(q, "Term", "-{U}"; In="DomainC")

    problem.postprocessing = postprocessing

    # PostOperation section
    postoperation = PostOperation()

    # Add post-operation items
    po1 = add!(postoperation, "Field_Maps", "Darwin_a_2D")
    op1 = add_operation!(po1)

    add_operation!(op1, "Print[ az, OnElementsOf Domain_Mag, Smoothing 1, Name \"flux lines: Az [T m]\", File StrCat[ \"$(joinpath(output_dir,"az_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ b, OnElementsOf Domain_Mag, Smoothing 1, Name \"B [T]\", File StrCat[ \"$(joinpath(output_dir,"b_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ bm, OnElementsOf Domain_Mag, Smoothing 1, Name \"|B| [T]\", File StrCat[ \"$(joinpath(output_dir,"bm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ jz, OnElementsOf Region[{DomainC}], Smoothing 1, Name \"jz [A/m²] Conducting domain\", File StrCat[ \"$(joinpath(output_dir,"jz_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ rhoj2, OnElementsOf Region[{DomainC}], Smoothing 1, Name \"Power density\", File StrCat[ \"$(joinpath(output_dir,"rhoj2_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ jm, OnElementsOf DomainC, Smoothing 1, Name \"|j| [A/m²] Conducting domain\", File StrCat[ \"$(joinpath(output_dir,"jm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")
    add_operation!(op1, "Print[ dm, OnElementsOf DomainC, Smoothing 1, Name \"|D| [A/m²]\", File StrCat[ \"$(joinpath(output_dir,"dm_"))\", Sprintf(\"%g\",active_con), \".pos\" ] ];")

    po2 = add!(postoperation, "LineParams", "Darwin_a_2D")
    op2 = add_operation!(po2)
    add_operation!(op2, "Print[ Z, OnRegion Conductors, Format Table, File \"$(joinpath(output_dir,"Z.dat"))\", AppendToExistingFile (active_con > 1 ? 1 : 0) ];")

    # Add the post-operation to the problem
    problem.postoperation = postoperation

end


function run_getdp(workspace::FEMWorkspace, fem_formulation::AbstractFormulationSet)
    # Initialize Gmsh if not already initialized
    if gmsh.is_initialized() == 0
        gmsh.initialize()
    end

    # Number of iterations (from the original function)
    n_phases = sum([length(c.design_data.components) for c in workspace.problem_def.system.cables])

    # Flag to track if all solves are successful
    all_success = true

    # Map verbosity to Gmsh/GetDP level
    gmsh_verbosity = map_verbosity_to_gmsh(workspace.opts.verbosity)
    gmsh.option.set_number("General.Verbosity", gmsh_verbosity)


    getdp_verbosity = map_verbosity_to_getdp(workspace.opts.verbosity)

    # Loop over each active_ind from 1 to n_phases
    for i in 1:n_phases
        # Construct solver command with -setnumber active_ind i
        solve_cmd = "$(workspace.opts.getdp_executable) $(fem_formulation.problem.filename) -msh $(workspace.paths[:mesh_file]) -solve $(fem_formulation.resolution_name) -setnumber active_con $i -v2 -verbose $(getdp_verbosity)"

        # Log the current solve attempt
        @info "Solving for source conductor $i... (Resolution = $(fem_formulation.resolution_name))"

        # Attempt to run the solver
        try
            gmsh.onelab.run("GetDP", solve_cmd)

            if workspace.opts.plot_field_maps
                @info "Building field maps for source conductor $i... (Resolution = $(fem_formulation.resolution_name))"

                post_cmd = "$(workspace.opts.getdp_executable) $(fem_formulation.problem.filename) -msh $(workspace.paths[:mesh_file]) -pos Field_Maps -setnumber active_con $i -v2 -verbose $(getdp_verbosity)"

                gmsh.onelab.run("GetDP", post_cmd)
            end

            @info "Solve successful for source conductor $(i)!"
        catch e
            # Log the error and update the success flag
            @error "Solver failed for source conductor $i: $e"
            all_success = false
            # Continue to the next iteration even if this one fails
        end
    end

    # Return true only if all solves were successful
    return all_success
end

function run_solver!(workspace::FEMWorkspace)

    # Get problem and formulation from the workspace
    problem = workspace.problem_def
    formulation = workspace.formulation

    # Preallocate result matrices
    n_phases = sum(length(c.design_data.components) for c in problem.system.cables)
    n_frequencies = length(problem.frequencies)

    Z = zeros(ComplexF64, n_phases, n_phases, n_frequencies)
    Y = zeros(ComplexF64, n_phases, n_phases, n_frequencies)

    # Solve for each frequency
    for (freq_idx, frequency) in enumerate(problem.frequencies)
        @info "Solving frequency $freq_idx/$n_frequencies: $frequency Hz"

        try
            _do_run_solver!(frequency, freq_idx, workspace, Z, Y)
        catch e
            @error "Solver failed for frequency $frequency Hz" exception = e
            rethrow(e)
        end

        # Archive results if not cleaning up
        if workspace.opts.keep_run_files
            archive_frequency_results(workspace, frequency)
        end
    end

    return LineParameters(Z, Y)
end

function _do_run_solver!(frequency::Float64, freq_idx::Int,
    workspace::FEMWorkspace,
    Z::Array{ComplexF64,3}, Y::Array{ComplexF64,3})

    # Get formulation from workspace
    formulation = workspace.formulation

    # Build and solve both formulations
    for fem_formulation in formulation.analysis_type
        @debug "Processing $(fem_formulation.resolution_name) formulation"

        make_fem_problem!(fem_formulation, frequency, workspace)

        if !run_getdp(workspace, fem_formulation)
            Base.error("$(fem_formulation.resolution_name) solver failed")
        end
    end

    # Extract results into preallocated arrays
    Z[:, :, freq_idx] = read_results_file(formulation.analysis_type[1], workspace)
    Y[:, :, freq_idx] = read_results_file(formulation.analysis_type[2], workspace)
end

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
    workspace::Union{FEMWorkspace,Nothing}=nothing)

    opts = formulation.options

    # Initialize workspace
    workspace = init_workspace(problem, formulation, workspace)

    # Meshing phase: make_mesh! decides if it needs to run.
    # It returns true if the process should stop (e.g., mesh_only=true).
    if make_mesh!(workspace)
        return workspace, nothing
    end

    # Solving phase - always runs unless mesh_only
    @info "Starting FEM solver"
    line_params = run_solver!(workspace)
    # line_params = run_solver!(problem, formulation, workspace)

    @info "FEM computation completed successfully"
    return workspace, line_params
end

