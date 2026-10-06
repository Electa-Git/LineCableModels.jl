const FEM_MESH_CONTROL_FIELDS = (
    :domain_size_factor,:pml_reflection,:mesh_size_factor,
    (spec.field for spec in FEM_NATIVE_OVERRIDES if spec.mesh)...)
const FEM_EXPORT_MESH_GLOBALS = (
    "Mesh.MshFileVersion", "Mesh.Binary", "Mesh.SaveAll", "Mesh.MeshSizeMin", "Mesh.MeshSizeMax",
    "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeExtendFromBoundary")

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
            println(io,"  ",name,"(FEMCase) = DefineNumber[",name,"(FEMCase), Name Sprintf[\"Inputs/Cases/%04g/",replace(label,r" \[[^]]*/[^]]*\]"=>""),"\",FEMCase+1], Label \"",label[3:end],"\"];")
        end
        println(io,"EndFor")
        for (name,label) in (("AirSigma","01Sigma [S/m]"),("AirEpsilon","02Epsilon [F/m]"),("AirMu","03Mu [H/m]"))
            println(io,name," = DefineNumber[",name,", Name \"Inputs/Air/",replace(label,r" \[[^]]*/[^]]*\]"=>""),"\", Label \"",label[3:end],"\"];")
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
            println(io,name,"_Mu = DefineNumber[",name,"_Mu, Name \"Materials/",name,"/01Mu\", Label \"Mu [H/m]\"];")
            println(io,"For FEMCase In {0:FrequencyCount-1}")
            for (coefficient,label) in (("Sigma","02Sigma [S/m]"),("Epsilon","03Epsilon [F/m]"))
                println(io,"  ",name,"_",coefficient,"(FEMCase) = DefineNumber[",name,"_",coefficient,"(FEMCase), Name Sprintf[\"Materials/",name,"/Cases/%04g/",replace(label,r" \[[^]]*/[^]]*\]"=>""),"\",FEMCase+1], Label \"",label[3:end],"\"];")
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
            "RunAction = {1, Choices{0=\"Mesh only\",1=\"Mesh and solve\"}, Name \"Inputs/04Run action\"}",
            "MeshSizeFactor = {$(_pro_number(controls.mesh_size_factor)), Min 0., Name \"Physics/MeshSizeFactor\"}",
            "DomainSizeFactor = {$(_pro_number(controls.domain_size_factor)), Min 0., Name \"Physics/DomainSizeFactor\"}",
            "PmlReflection = {$(_pro_number(controls.pml_reflection)), Min 0., Max 1., Name \"Physics/PmlReflection\"}",
            "LinearSolver = {$(Int(controls.linear_solver === :gmres)), Choices{0=\"Direct MUMPS\",1=\"GMRES with LU preconditioning\"}, Name \"Solver/LinearSolver\"}",
            "GmresIterationsMax = {$(controls.gmres_iterations_max), Min 1, Step 1, Name \"Solver/GmresIterationsMax\"}",
            "GmresRelativeTolerance = {$(_pro_number(controls.gmres_relative_tolerance)), Name \"Solver/GmresRelativeTolerance\"}",
            "GmresAbsoluteTolerance = {$(_pro_number(controls.gmres_absolute_tolerance)), Min 0, Name \"Solver/GmresAbsoluteTolerance\"}",
            "PlotFieldMaps = {$(Int(controls.plot_field_maps)), Choices{0,1}, Name \"Outputs/PlotFieldMaps\"}",
            "GmshVerbosity = {$(controls.gmsh_verbosity), Min -1, Max 5, Step 1, Name \"Diagnostics/GmshVerbosity\"}",
            "GetDPVerbosity = {$(controls.getdp_verbosity), Min -1, Max 5, Step 1, Name \"Diagnostics/GetDPVerbosity\"}",
            "GetDPThreads = {$(controls.solver_threads), Min 1, Step 1, Name \"Resources/GetDPThreads\"}",
            "OutputTotal = {$(Int(controls.output_basis === Val(:total))), Choices{0=\"per metre\",1=\"total length\"}, Name \"Outputs/OutputTotal\"}")
        expert_controls=String[]
        for spec in FEM_NATIVE_OVERRIDES
            value=getproperty(controls,spec.field)
            values=length(spec.names)==1 ? (value,) : value
            for (name,v) in zip(spec.names,values)
                limits=spec.closed === :choices ? "Choices{$(join(spec.range,','))}, " :
                    "Min $(spec.closed ? first(spec.range) : nextfloat(Float64(first(spec.range)))), " * (isfinite(last(spec.range)) ? "Max $(last(spec.range)), " : "")
                step=spec.kind === :integer ? "Step 1, " : ""
                push!(expert_controls,"$name = {$(_pro_number(v)), $limits$(step)Name \"Expert/$name\"}")
            end
        end
        all_controls=[collect(ui_controls);expert_controls]
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
    If(GmshVerbosity >= 0) General.Verbosity = GmshVerbosity; EndIf
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
$(TYPEDSIGNATURES)

Export a detached ONELAB bundle using the same native sources and option validator
as managed FEM computation. `options=(;)` accepts the physics controls, native
`overrides`, solver controls, `solver_threads`, `plot_field_maps`, `output_basis`,
and native verbosity. Managed workers, retention, resume, callbacks and logging
controls are rejected. The main ONELAB groups contain physics and solver controls;
all expert controls use their native names in the `Expert` group.

`file_name` must end in `.pro`. `overwrite=true` replaces only files owned by an
existing exported bundle. Returns the entry file path.
"""
function ImportExport.export_data(::Val{:onelab}, problem::LineParametersProblem,
        formulation::LineCableModelsFEM; file_name::AbstractString,
        options::NamedTuple=(;), overwrite::Bool=false)
    managed_keys=(:frequency_workers,:mesh_policy,:mesh_path,:keep_run_directory,
        :resume_run_directory,:on_result,:log_file,:trace,:timing,:verbosity,:getdp_executable)
    rejected=intersect(keys(options),managed_keys)
    isempty(rejected) || throw(ArgumentError("ONELAB export does not accept managed-run options: $(Tuple(rejected))"))
    opts=computation_options(LineCableModelsFEM,ComputationOptions(options))
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
        # Omitted verbosity inherits the existing detached native stage defaults.
        # Explicit values have already passed the shared computation validator.
        controls=merge(opts.data,(gmsh_verbosity=get(options,:gmsh_verbosity,-1),
            getdp_verbosity=get(options,:getdp_verbosity,-1)))
        _write_export_bundle(stage,stem,model,formulation,controls)
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
