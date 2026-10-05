const FEM_EXPORT_MESH_OPTIONS = (
    :domain_size_factor, :pml_thickness_factor, :pml_layers, :pml_grading,
    :pml_reflection,
    :mesh_size_factor, :exterior_mesh_size_factor,
    :interface_refinement_factor, :volume_quadrature,
    :physical_volume_quadrature, :pml_element_family, :pml_quadrature,
    :conductor_geometry_tolerance, :conductor_skin_depth_elements, :conductor_mesh_growth,
    :conductor_skin_depths, :conductor_thickness_elements)
const FEM_EXPORT_MESH_GLOBALS = (
    "Mesh.MshFileVersion", "Mesh.Binary", "Mesh.SaveAll", "Mesh.MeshSizeMin", "Mesh.MeshSizeMax",
    "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary")
const FEM_EXPORT_SOLVER_OPTIONS = (:mumps_ordering, :petsc_prealloc,
    :linear_solver, :gmres_iterations_max, :gmres_relative_tolerance, :gmres_absolute_tolerance,
    :mumps_error_analysis, :mumps_refinement_max,
    :mumps_backward_error_tolerance, :mumps_forward_error_tolerance)

_export_identifier(value) = replace(String(value), r"[^a-zA-Z0-9_]" => "_")
_export_material_name(index, material) = "Material_$(index)_$(_export_identifier(material.kind))"

function _write_export_data(path, model, formulation, controls)
    # Retain the shared input contract, with named native coefficients as its
    # only values. The arrays below are GetDP's frequency/region indexing.
    _write_model_data(path, model, controls)
    original = readlines(path)
    material_lines = ("MaterialRegionTags(", "MaterialIsConductor(",
        "MaterialHasLoss(", "MaterialMu(", "MaterialSigma_", "MaterialEpsilon_")
    open(path, "w") do io
        println(io, "// Authoritative native model data. SI units; exp(+j omega t).")
        println(io, "// Edit these coefficients directly. Native execution never rewrites this file.")
        native_controls = ("MeshSizeFactor", "ExteriorMeshSizeFactor", "InterfaceRefinementFactor", "PmlQuadrangles",
            "VolumeQuadrature", "PhysicalVolumeQuadrature", "PmlQuadrature",
            "ConductorGeometryTolerance", "ConductorSkinDepthElements", "ConductorMeshGrowth", "ConductorSkinDepths", "ConductorThicknessElements",
            "DomainSizeFactor", "PmlReflection", "PmlSideThicknessFactor", "PmlTopThicknessFactor", "PmlBottomThicknessFactor",
            "PmlSideLayers", "PmlTopLayers", "PmlBottomLayers", "PmlSideGrading", "PmlTopGrading", "PmlBottomGrading")
        for line in original
            any(prefix -> startswith(line,prefix),material_lines) && continue
            any(name -> startswith(line,name*" ="),native_controls) && continue
            println(io,line)
        end
        println(io,"// Indexed native values remain authoritative after selection and ONELAB edits.")
        println(io,"For FEMCase In {0:FrequencyCount-1}")
        for (name,label) in (("Frequencies","01Frequency [Hz]"),("GammaReValues","02Gamma real [1/m]"),
                ("GammaImValues","03Gamma imag [1/m]"),("EarthSigma","04Soil sigma [S/m]"),
                ("EarthEpsilon","05Soil epsilon [F/m]"),("EarthMu","06Soil mu [H/m]"))
            println(io,"  ",name,"(FEMCase) = DefineNumber[",name,"(FEMCase), Name Sprintf[\"Inputs/Cases/%04g/",label,"\",FEMCase+1]];")
        end
        println(io,"EndFor")
        for (name,label) in (("AirSigma","01Sigma [S/m]"),("AirEpsilon","02Epsilon [F/m]"),("AirMu","03Mu [H/m]"))
            println(io,name," = DefineNumber[",name,", Name \"Inputs/Air/",label,"\"];")
        end
        println(io,"FieldMapsHelmholtz() = Str[",
            join(_pro_string.(_field_quantities(formulation.options.data.physics)),", "),"];")
        println(io, "// Terminal encounter order and connections: positive phase IDs; 0 grounded.")
        for (i, name) in enumerate(model.terminal_ids)
            println(io, "// Terminal ", i, ": ", replace(name, '\n'=>' '))
            println(io, "// Terminal~{", i, "}: surface tag ", model.tags.terminal_base+i,
                "; TerminalContour~{", i, "}: curve tag ", model.tags.terminal_contour_base+i)
            println(io, "Connection_", i, " = ", model.problem.system.connection_order[i], ";")
        end
        println(io, "Connections() = {", join(["Connection_$i" for i in eachindex(model.terminal_ids)], ", "), "};")
        for (i, m) in enumerate(model.material_plans)
            name = _export_material_name(i, m)
            println(io, "\n// ", m.object_id, "; ", m.physical_name, "; selected law: ", m.field)
            println(io, name, "_Tag = ", m.physical_tag, ";")
            println(io, name, "_Mu = ", _pro_number(m.mu_r*4π*1e-7), "; // H/m")
            println(io, name, "_Sigma() = ", _pro_array(real.(m.admittivity)), "; // S/m, by frequency")
            println(io, name, "_Epsilon() = ", _pro_array(imag.(m.admittivity)./(2π.*model.problem.frequencies)), "; // F/m")
            println(io,name,"_Mu = DefineNumber[",name,"_Mu, Name \"Materials/",name,"/01Mu [H/m]\"];")
            println(io,"For FEMCase In {0:FrequencyCount-1}")
            for (coefficient,label) in (("Sigma","02Sigma [S/m]"),("Epsilon","03Epsilon [F/m]"))
                println(io,"  ",name,"_",coefficient,"(FEMCase) = DefineNumber[",name,"_",coefficient,"(FEMCase), Name Sprintf[\"Materials/",name,"/Cases/%04g/",label,"\",FEMCase+1]];")
            end
            println(io,"EndFor")
            println(io, "MaterialSigma_", i, "() = ", name, "_Sigma();")
            println(io, "MaterialEpsilon_", i, "() = ", name, "_Epsilon();")
            println(io, name, "_HasLoss = 0;")
            println(io, "For MaterialCase In {0:FrequencyCount-1}")
            println(io, "  If(",name,"_Sigma(MaterialCase) != 0) ",name,"_HasLoss = 1; EndIf")
            println(io, "EndFor")
        end
        names = [_export_material_name(i,m) for (i,m) in enumerate(model.material_plans)]
        println(io, "MaterialRegionTags() = {", join(names .* "_Tag", ", "), "};")
        println(io, "MaterialMu() = {", join(names .* "_Mu", ", "), "};")
        println(io, "MaterialIsConductor() = ", _pro_array([m.kind===:conductor for m in model.material_plans]), ";")
        println(io, "MaterialHasLoss() = {", join(names .* "_HasLoss", ", "), "};")

        println(io, "\n// Current amplitudes AND coefficient normalization used by the native equations.")
        println(io, "UnitSource = 1.; // axial A; constraints and normalized Z/P use this value")
        println(io, "LineLength = ", _pro_number(model.problem.system.line_length), "; // m")
        for (name, key) in (("ReduceBundle", :reduce_bundle), ("KronReduction", :kron_reduction), ("IdealTransposition", :ideal_transposition))
            println(io, name, " = ", Int(getproperty(formulation.options.data,key)), ";")
        end
        println(io, "DefineConstant[")
        frequency_choices = join(["$i=" * _pro_string(@sprintf("%.8g Hz", f))
            for (i,f) in enumerate(model.problem.frequencies)], ",")
        ui_controls = (
            "RunFrequencyScan = {0, Choices{0,1}, Name \"Inputs/00Run frequency scan\"}",
            "ScanFrequencyIndex = {1, Choices{1:FrequencyCount}, Min 1, Max FrequencyCount, Step 1, Loop RunFrequencyScan, Visible 0, Name \"Inputs/00Scan frequency index\"}",
            "FrequencyIndex = {1, Choices{$frequency_choices}, Visible !RunFrequencyScan, Name \"Inputs/01Frequency case\"}",
            "Physics = {$(_fem_physics_code(formulation)), Choices{1=\"Helmholtz\"}, Name \"Inputs/02Physics\"}",
            "BasisTerminal = {0, Min 0, Max NumTerminals, Step 1, Name \"Inputs/03Basis (0 = full matrix)\"}",
            "MeshSizeFactor = {$(_pro_number(controls.mesh_size_factor)), Name \"Mesh/01Physical size factor\"}",
            "ExteriorMeshSizeFactor = {$(_pro_number(controls.exterior_mesh_size_factor)), Min 1., Name \"Mesh/02Exterior size factor\"}",
            "InterfaceRefinementFactor = {$(_pro_number(controls.interface_refinement_factor)), Min 1., Name \"Mesh/09Interface footprint factor\"}",
            "ConductorGeometryTolerance = {$(_pro_number(controls.conductor_geometry_tolerance)), Name \"Mesh/03Conductor geometry tolerance\"}",
            "ConductorSkinDepthElements = {$(_pro_number(controls.conductor_skin_depth_elements)), Name \"Mesh/04Conductor skin depth elements\"}",
            "ConductorMeshGrowth = {$(_pro_number(controls.conductor_mesh_growth)), Name \"Mesh/05Conductor normal growth\"}",
            "ConductorSkinDepths = {$(_pro_number(controls.conductor_skin_depths)), Name \"Mesh/06Conductor graded skin depths\"}",
            "ConductorThicknessElements = {$(controls.conductor_thickness_elements), Min 1, Step 1, Name \"Mesh/07Conductor wall elements\"}",
            "RunAction = {1, Choices{0=\"Mesh only\",1=\"Mesh and solve\"}, Name \"Inputs/04Run action\"}",
            "PlotFieldMaps = {0, Choices{0,1}, Name \"Outputs/01Write field maps\"}",
            "GetDPThreads = {1, Min 1, Step 1, Name \"Numerics/02GetDP threads\"}",
            "MumpsOrdering = {$(something(controls.mumps_ordering,-1)), Choices{-1=\"solver default\",0=\"AMD\",2=\"AMF\",3=\"Scotch\",4=\"PORD\",5=\"METIS\",6=\"QAMD\",7=\"automatic\"}, Name \"Numerics/03MUMPS ordering\"}",
            "PetscPrealloc = {$(something(controls.petsc_prealloc,0)), Min 0, Step 1, Name \"Numerics/04Sparse row allocation (0 = default)\"}",
            "MumpsErrorAnalysis = {$(controls.mumps_error_analysis), Choices{0=\"off\",1=\"full sensitivity estimates\",2=\"backward errors\"}, Name \"Numerics/08MUMPS error analysis\"}",
            "MumpsRefinementMax = {$(controls.mumps_refinement_max), Min 0, Step 1, Name \"Numerics/09Maximum refinement steps\"}",
            "MumpsBackwardErrorTolerance = {$(_pro_number(controls.mumps_backward_error_tolerance)), Name \"Numerics/10Backward-error target\"}",
            "MumpsForwardErrorTolerance = {$(_pro_number(controls.mumps_forward_error_tolerance)), Name \"Numerics/11Forward-error comparison budget\"}",
            "LinearSolver = {$(Int(controls.linear_solver === :gmres)), Choices{0=\"Direct MUMPS\",1=\"GMRES with LU preconditioning\"}, Name \"Numerics/00Linear solver\"}",
            "GmresIterationsMax = {$(controls.gmres_iterations_max), Min 1, Step 1, Name \"Numerics/12GMRES maximum iterations\"}",
            "GmresRelativeTolerance = {$(_pro_number(controls.gmres_relative_tolerance)), Name \"Numerics/13GMRES relative residual tolerance\"}",
            "GmresAbsoluteTolerance = {$(_pro_number(controls.gmres_absolute_tolerance)), Min 0, Name \"Numerics/14GMRES absolute residual tolerance\"}",
            "VolumeQuadrature = {$(controls.volume_quadrature), Choices{4,7,12,13}, Name \"Numerics/05Triangle volume points\"}",
            "PhysicalVolumeQuadrature = {$(something(controls.physical_volume_quadrature,0)), Choices{0=\"inherit triangle rule\",3=\"3\",4=\"4\",7=\"7\",12=\"12\",13=\"13\"}, Name \"Numerics/06Physical triangle points\"}",
            "PmlQuadrangles = {$(Int(controls.pml_element_family === :quadrangle)), Choices{0=\"triangles\",1=\"quadrangles\"}, Name \"Mesh/10PML element family\"}",
            "PmlQuadrature = {$(controls.pml_quadrature), Choices{4,9,16}, Name \"Numerics/07PML quadrangle points\"}",
            "OutputTotal = {0, Choices{0=\"per metre\",1=\"total length\"}, Name \"Outputs/02Matrix basis\"}")
        pml_controls = String[
            "DomainSizeFactor = {$(_pro_number(controls.domain_size_factor)), Name \"Boundary/01Domain size factor [dimensionless]\"}",
            "PmlReflection = {$(_pro_number(controls.pml_reflection)), Name \"Boundary/02Nominal PML reflection\"}"]
        for (direction,name) in enumerate(("Side","Top","Bottom"))
            factors = controls.pml_thickness_factor
            factor = factors isa Tuple ? factors[direction] : factors
            push!(pml_controls,"Pml$(name)ThicknessFactor = {$(_pro_number(factor)), Name \"Boundary/$(name)/01Relative thickness\"}",
                "Pml$(name)Layers = {$(controls.pml_layers[direction]), Min 1, Step 1, Name \"Boundary/$(name)/02Minimum intervals\"}",
                "Pml$(name)Grading = {$(_pro_number(controls.pml_grading[direction])), Min 0, Name \"Boundary/$(name)/03Grading exponent\"}")
        end
        all_controls = [pml_controls;collect(ui_controls)]
        for (i,line) in enumerate(all_controls)
            println(io, "  ", line, i < length(all_controls) ? "," : "")
        end
        println(io, "];")
        println(io, "// ONELAB loops the complete Gmsh -> GetDP action over the hidden index.")
        println(io, "// Keep the manual selection intact, including after the loop resets its index.")
        println(io, "If(RunFrequencyScan && !StrCmp(OnelabAction, \"compute\"))")
        println(io, "  FrequencyIndex = ScanFrequencyIndex;")
        println(io, "EndIf")
    end
end

function _write_export_bundle(root, stem, model, formulation, controls)
    mkpath(joinpath(root,"geometry"))
    mkpath(joinpath(root,"formulations"))
    for (key,path) in pairs(_getdp_assets(joinpath(root,"formulations")))
        write(path,getproperty(FEM_GETDP_SOURCES,key))
    end
    _write_export_data(joinpath(root,stem*"_data.pro"),model,formulation,controls)
    lock(FEM_SESSION_LOCK) do
        session = _start_gmsh(0)
        saved = Dict(option=>gmsh.option.get_number(option) for option in FEM_EXPORT_MESH_GLOBALS)
        try
            physical = _build_physical_geometry!(model,"export-$(basename(root))")
            _write_physical_geometry(joinpath(root,"geometry","physical.geo"),model,physical)
            gmsh.model.remove()
        finally
            for (option,value) in saved
                gmsh.option.set_number(option,value)
            end
            _finish_gmsh(session)
        end
    end
    write(joinpath(root,stem*".geo"), """
    // Invalidate the displayed mesh identity before parsing any case inputs.
    MeshPublished = DefineString["No mesh", Name "Mesh/Current mesh/00Status", ReadOnly 1];
    MeshPublishedCase = DefineNumber[0, Name "Mesh/Current mesh/01Case index", ReadOnly 1];
    MeshPublishedFrequency = DefineNumber[0, Name "Mesh/Current mesh/02Frequency [Hz]", ReadOnly 1];
    Solver.AutoCheck = 1;
    Include "$(stem)_data.pro";
    Include "formulations/parameters.pro";
    Include "geometry/physical.geo";
    Include "formulations/geometry.geo";
    Include "formulations/mesh.geo";
    // Remove only this project's derived views when inputs are checked/rerun.
    If(!StrCmp(OnelabAction, "check") || !StrCmp(OnelabAction, "compute"))
      Include "views.geo";
    EndIf
    // ONELAB Check/Run and explicit CLI BuildMesh use the same native mesh operations.
    Solver.AutoMesh = -1;
    If(!StrCmp(OnelabAction, "check") || !StrCmp(OnelabAction, "compute") || Exists(BuildMesh))
      Mesh 2;
      Save StrCat(CurrentDirectory, "$(stem).msh");
      MeshPublishedCase = DefineNumber[FrequencyIndex, Name "Mesh/Current mesh/01Case index", ReadOnly 1];
      MeshPublishedFrequency = DefineNumber[FrequencyHz, Name "Mesh/Current mesh/02Frequency [Hz]", ReadOnly 1];
      MeshPublished = DefineString["Generated", Name "Mesh/Current mesh/00Status", ReadOnly 1];
    EndIf
    """)
    write(joinpath(root,stem*".pro"), """
    // Open this file in Gmsh/ONELAB. All numerical operations use native GetDP.
    // Geometry: $(stem).geo; coefficients/connections: $(stem)_data.pro.
    // Regions, constraints and equations: formulations/helmholtz.pro.
    ProjectDirectory = CurrentDirectory;
    ModelDataPath = StrCat[ProjectDirectory,"$(stem)_data.pro"];
    Include ModelDataPath;
    Include "formulations/onelab.pro";
    """)
    # Commands use relative, shell-quoted filenames, including spaces and apostrophes.
    command_stem = replace("./" * stem, "'" => "'\\''")
    readme = replace(read(joinpath(@__DIR__,"onelab_export","README.md"),String),
        "{{MODEL_NAME}}" => stem, "{{COMMAND_STEM}}" => command_stem)
    write(joinpath(root,"README.md"),readme)
    cp(joinpath(@__DIR__,"onelab_export","views.geo"),joinpath(root,"views.geo"))
    files = sort!([relpath(joinpath(dir,name),root) for (dir,_,names) in walkdir(root) for name in names])
    write(joinpath(root,".onelab-export-files"),join(files,"\n")*"\n")
end

"""
    export_data(:onelab, problem, formulation; file_name, mesh_options=(;), solver_options=(;), overwrite=false)

Export a detached Gmsh/ONELAB project with readable native geometry and GetDP
model files. Returns the absolute entry `.pro` path. Load `Gmsh` before calling.

`file_name` selects the entry file in a dedicated bundle directory.
`mesh_options` accepts `domain_size_factor`,
`pml_thickness_factor`, `pml_layers`, `pml_grading`, `pml_reflection`,
`mesh_size_factor`,
`exterior_mesh_size_factor`, `interface_refinement_factor`, `volume_quadrature`,
`physical_volume_quadrature`, `pml_element_family`, and `pml_quadrature`, with their existing FEM
meanings and defaults. `domain_size_factor` multiplies the physical earth sizing
length `min(1/real(q0_e),2pi/abs(q0_e))` \\[m\\] after the native reference-earth ceiling
and prescribed-Γ adjustment.
`solver_options` accepts `linear_solver`, `gmres_iterations_max`,
`gmres_relative_tolerance`, `gmres_absolute_tolerance`,
`mumps_ordering`, `petsc_prealloc`, `mumps_error_analysis`,
`mumps_refinement_max`, `mumps_backward_error_tolerance`, and
`mumps_forward_error_tolerance`
with the same defaults as FEM computation options. These are exported as editable
ONELAB controls; ordering availability depends on the native solver build.
Direct mode prints native MUMPS estimates; GMRES prints convergence reasons
and estimated and true residuals with MUMPS refinement and error analysis disabled.
Detached execution leaves their interpretation to the user. Managed Julia
execution also emits target-based warnings.
Physical CAD is fixed at export. Native frequency, Γ, material and mesh controls
recompute exterior geometry, mesh constraints and solver coefficients together.
The detached runtime uses only native Gmsh/GetDP with ONELAB; it requires neither
Julia nor Python. Export does not execute GetDP or open a GUI. `overwrite=true`
replaces recorded bundle files and removes obsolete owned files, preserving
unrelated files in the destination.
"""
function ImportExport.export_data(::Val{:onelab}, problem::LineParametersProblem,
        formulation::LineCableModelsFEM; file_name::AbstractString,
        mesh_options::NamedTuple=(;), solver_options::NamedTuple=(;), overwrite::Bool=false)
    opts = computation_options(LineCableModelsFEM,ComputationOptions(;mesh_options...,solver_options...))
    unknown = setdiff(keys(mesh_options),FEM_EXPORT_MESH_OPTIONS)
    isempty(unknown) || throw(ArgumentError("Unsupported export mesh options: $(join(unknown, ", "))"))
    unknown_solver = setdiff(keys(solver_options),FEM_EXPORT_SOLVER_OPTIONS)
    isempty(unknown_solver) || throw(ArgumentError("Unsupported export solver options: $(join(unknown_solver, ", "))"))
    entry = abspath(file_name)
    stem, extension = splitext(basename(entry))
    extension == ".pro" && !isempty(stem) || throw(ArgumentError("file_name must end in .pro"))
    occursin(r"[\r\n\"]",entry) && throw(ArgumentError("export path cannot contain quotes or newlines"))
    root = dirname(entry)
    existing = isdir(root) ? readdir(root) : String[]
    marker = joinpath(root,".onelab-export-files")
    !isempty(existing) && (!overwrite || !isfile(marker)) && throw(ArgumentError("destination is not empty; overwrite requires an existing exported bundle"))
    model = _resolved_fem_model(_preflight_fem_problem(problem),formulation)
    mkpath(dirname(root))
    stage = mktempdir(dirname(root);prefix=".onelab-export-")
    try
        _write_export_bundle(stage,stem,model,formulation,opts.data)
        owned = isfile(marker) ? Set(readlines(marker)) : Set{String}()
        for file in owned
            normalized = normpath(file)
            (isabspath(file) || first(splitpath(normalized)) == ".." || normalized == ".") &&
                throw(ArgumentError("invalid path in export ownership manifest: $file"))
            parent = dirname(joinpath(root,normalized))
            while parent != root
                islink(parent) && throw(ArgumentError("export ownership path traverses a symbolic link: $file"))
                parent = dirname(parent)
            end
        end
        files = [readlines(joinpath(stage,".onelab-export-files")); ".onelab-export-files"]

        for file in files
            target = joinpath(root,file)
            ispath(target) && !(file in owned || file==".onelab-export-files") && throw(ArgumentError("refusing to replace unrelated destination file $target"))
        end
        mkpath(root)
        for file in setdiff(owned,Set(files))
            target = joinpath(root,file)
            (isfile(target) || islink(target)) && rm(target)
        end
        for file in files
            mkpath(dirname(joinpath(root,file)))
            cp(joinpath(stage,file),joinpath(root,file);force=true)
        end
    finally
        rm(stage;recursive=true,force=true)
    end
    return entry
end

function ImportExport.export_data(format::Val{:onelab}, system::DataModel.LineCableSystem,
        formulation::LineCableModelsFEM; earth_props, frequencies,
        temperature=20.0, kwargs...)
    problem = LineParametersProblem(system;earth_props,frequencies,temperature)
    return ImportExport.export_data(format,problem,formulation;kwargs...)
end
