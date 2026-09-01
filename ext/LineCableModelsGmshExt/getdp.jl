const FEM_FIELD_QUANTITIES = (
    "az", "b", "bm", "e", "ez", "em", "jz", "jm", "rhoj2"
)

const FEM_FIELD_SPECS = (
    (quantity = "az", label = "Az [T m]", domain = "Domain_Mag"),
    (quantity = "b", label = "B [T]", domain = "Domain_Mag"),
    (quantity = "bm", label = "|B| [T]", domain = "Domain_Mag"),
    (quantity = "e", label = "E [V/m]", domain = "Domain_Mag"),
    (quantity = "ez", label = "Ez [V/m]", domain = "Domain_Mag"),
    (quantity = "em", label = "|E| [V/m]", domain = "Domain_Mag"),
    (quantity = "jz", label = "Jz [A/m2]", domain = "DomainC"),
    (quantity = "jm", label = "|J| [A/m2]", domain = "DomainC"),
    (quantity = "rhoj2", label = "S [W/m3]", domain = "DomainC")
)

function _pro_string(value::AbstractString)
    escaped = replace(String(value), '\\' => "\\\\", '"' => "\\\"")
    return "\"$escaped\""
end

function _pro_number(value::Real)
    return @sprintf("%.17g", Float64(value))
end

_pro_array(values) = "{" * join(_pro_number.(values), ", ") * "}"

function _write_problem_snapshot(path::String, problem::LineParametersProblem)
    open(path, "w") do io
        JSON3.pretty(io, ImportExport.serialize_value(problem))
        write(io, '\n')
    end
    return path
end

function _write_model_data(path::String, model::FEMResolvedModel)
    earth = model.problem.earth_props.layers[2]
    material_tags = getproperty.(model.material_plans, :physical_tag)
    material_conductors = Int[material.kind === :conductor
                              for material in model.material_plans]
    material_sigma = [isinf(material.rho) ? 0.0 : inv(material.rho)
                      for material in model.material_plans]
    material_epsilon = [material.eps_r * 8.8541878128e-12
                        for material in model.material_plans]
    material_mu = [material.mu_r * 4π * 1e-7 for material in model.material_plans]
    material_functions_path = joinpath(dirname(path), "material_functions.pro")
    open(material_functions_path, "w") do io
        println(io, "// Generated immutable material regions and functions")
        println(io, "Group {")
        for index in eachindex(material_tags)
            println(io, "  MaterialRegion", index, " = Region[{",
                material_tags[index], "}];")
        end
        conductor_names = ["MaterialRegion$index"
                           for index in eachindex(material_tags)
                           if material_conductors[index] == 1]
        passive_names = ["MaterialRegion$index"
                         for index in eachindex(material_tags)
                         if material_conductors[index] == 0]
        println(io, "  ConductorMaterialRegions = Region[{",
            join(conductor_names, ", "), "}];")
        println(io, "  PassiveMaterialRegions = Region[{",
            join(passive_names, ", "), "}];")
        println(io, "  DomainCWithI = Region[{Terminals}];")
        println(io,
            "  DomainC = Region[{ConductorMaterialRegions, Earth, EarthInf}];")
        println(io,
            "  DomainCC = Region[{Air, AirInf, PassiveMaterialRegions}];")
        println(io, "  Domain_Mag = Region[{DomainC, DomainCC}];")
        println(io, "}")
        println(io, "Function {")
        for index in eachindex(material_tags)
            println(io, "  nu[MaterialRegion", index, "] = ",
                _pro_number(inv(material_mu[index])), ";")
            println(io, "  sigma_dc[MaterialRegion", index, "] = ",
                _pro_number(material_sigma[index]), ";")
            println(io, "  epsilon[MaterialRegion", index, "] = ",
                _pro_number(material_epsilon[index]), ";")
            println(io, "  mu[MaterialRegion", index, "] = ",
                _pro_number(material_mu[index]), ";")
            println(io, "  tan_delta[MaterialRegion", index, "] = ",
                _pro_number(model.material_plans[index].tan_delta), ";")
        end
        println(io, "}")
    end
    field_map_dispatch_path = joinpath(dirname(path), "field_map_dispatch.pro")
    field_map_operations_path = joinpath(dirname(path), "field_map_operations.pro")
    maps_directory = joinpath(dirname(dirname(path)), "maps")
    open(field_map_dispatch_path, "w") do io
        println(io, "// Generated immutable field-map dispatch")
        println(io, "Macro FEMWriteMaps")
        for frequency in eachindex(model.problem.frequencies)
            for basis in eachindex(model.terminal_ids)
                operation = @sprintf("FEMFieldMaps_f%04d_b%04d", frequency, basis)
                println(io, "Test[\$FEMFrequencyIndex == ", frequency,
                    " && \$FEMBasisTerminal == ", basis, "]{")
                println(io, "  PostOperation[", operation, "];")
                println(io, "}{}")
            end
        end
        println(io, "Return")
    end
    open(field_map_operations_path, "w") do io
        println(io, "// Generated immutable field-map filenames")
        println(io, "PostOperation {")
        for frequency in eachindex(model.problem.frequencies)
            for basis in eachindex(model.terminal_ids)
                operation = @sprintf("FEMFieldMaps_f%04d_b%04d", frequency, basis)
                println(io, "  { Name ", operation,
                    "; NameOfPostProcessing FEMFields;")
                println(io, "    LastTimeStepOnly 1;")
                println(io, "    Operation {")
                for spec in FEM_FIELD_SPECS
                    filename = @sprintf("%s_f%04d_b%04d.pos", spec.quantity, frequency,
                        basis)
                    label = @sprintf("%s; f=%.17g Hz; basis=%s",
                        spec.label,
                        model.problem.frequencies[frequency],
                        model.terminal_ids[basis])
                    println(io, "      Print[", spec.quantity,
                        ", OnElementsOf ", spec.domain, ", Name ",
                        _pro_string(label), ", File ",
                        _pro_string(joinpath(maps_directory, filename)), "];")
                end
                println(io, "    }")
                println(io, "  }")
            end
        end
        println(io, "}")
    end
    open(path, "w") do io
        println(io, "// Generated immutable LineCableModels FEM model data")
        println(io, "NumTerminals = ", length(model.terminal_ids), ";")
        println(io, "NumCables = ", length(model.problem.system.designs), ";")
        println(io, "NumMaterialRegions = ", length(model.material_plans), ";")
        println(io, "FrequencyCount = ", length(model.problem.frequencies), ";")
        println(io, "AIR_EM = ", model.tags.air, ";")
        println(io, "EARTH_EM = ", model.tags.earth, ";")
        println(io, "AIR_INF = ", model.tags.air_infinite, ";")
        println(io, "EARTH_INF = ", model.tags.earth_infinite, ";")
        println(io, "DOMAIN_INF = ", model.tags.infinite_domain, ";")
        println(io, "OUTBND_EM = ", model.tags.outer_boundary, ";")
        println(io, "OUTBND_ELE_INS = ", model.tags.outer_air_boundary, ";")
        println(io, "OUTBND_ELE_REF = ", model.tags.outer_earth_boundary, ";")
        println(io, "INTERFACE_AIR_SOIL = ", model.tags.interface, ";")
        println(io, "INNER_INF_BND = ", model.tags.inner_shell_boundary, ";")
        println(io, "TERMINAL = ", model.tags.terminal_base + 1, ";")
        println(io, "TERMINAL_CONTOUR = ", model.tags.terminal_contour_base + 1, ";")
        println(io, "CABLE_CONTOUR = ", model.tags.cable_contour_base + 1, ";")
        println(io, "MaterialRegionTags() = ", _pro_array(material_tags), ";")
        println(io, "MaterialIsConductor() = ", _pro_array(material_conductors), ";")
        println(io, "MaterialSigma() = ", _pro_array(material_sigma), ";")
        println(io, "MaterialEpsilon() = ", _pro_array(material_epsilon), ";")
        println(io, "MaterialMu() = ", _pro_array(material_mu), ";")
        println(io, "MaterialTanDelta() = ",
            _pro_array(
                getproperty.(model.material_plans, :tan_delta)
            ), ";")
        println(io, "MaterialFunctionsPath = ",
            _pro_string(material_functions_path), ";")
        println(io, "FieldMapDispatchPath = ",
            _pro_string(field_map_dispatch_path), ";")
        println(io, "FieldMapOperationsPath = ",
            _pro_string(field_map_operations_path), ";")
        println(io, "sigma_earth = ", _pro_number(inv(earth.rho)), ";")
        println(io, "eps_earth = ", _pro_number(earth.eps_r * 8.8541878128e-12), ";")
        println(io, "mu_earth = ", _pro_number(earth.mu_r * 4π * 1e-7), ";")
        println(io, "eps0 = 8.8541878128e-12;")
        println(io, "mu0 = 1.2566370614359173e-6;")
        println(io, "UnitSource = 1.0;")
        println(io, "GammaQuasiTEMRe = 0.0;")
        println(io, "GammaQuasiTEMIm = 1.0e-12;")
        println(io, "If(!Exists(Val_Rint))")
        println(io, "  Val_Rint = ", _pro_number(model.domain_radius), ";")
        println(io, "EndIf")
        println(io, "If(!Exists(Val_Rext))")
        println(io, "  Val_Rext = ", _pro_number(model.shell_outer_radius), ";")
        println(io, "EndIf")
        println(io, "Xcenter = ", _pro_number(model.centre[1]), ";")
        println(io, "Ycenter = ", _pro_number(model.centre[2]), ";")
        println(io, "Zcenter = 0.0;")
    end
    return path
end

function _getdp_assets()
    root = joinpath(@__DIR__, "getdp")
    return (
        model = joinpath(root, "model.pro"),
        jacobian = joinpath(root, "jacobian_integration.pro"),
        quasi_tem = joinpath(root, "quasi_tem.pro")
    )
end

function _resolve_getdp(formulation::LineCableModelsFEM, run::FEMRun)
    explicit = formulation.execution.getdp_executable
    executable = explicit === nothing ? Sys.which("getdp") : abspath(explicit)
    executable === nothing && _fem_error(
        :getdp,
        "GetDP",
        :getdp_executable,
        "GetDP was not found; pass getdp_executable or add getdp to PATH";
        run_directory = run.path
    )
    isfile(executable) || _fem_error(
        :getdp,
        "GetDP",
        :getdp_executable,
        "GetDP executable does not exist: $executable";
        run_directory = run.path
    )
    output = try
        buffer = IOBuffer()
        Base.run(pipeline(
            Cmd([executable, "-info"]), stdout = buffer, stderr = buffer
        ))
        String(take!(buffer))
    catch exception
        _fem_error(
            :getdp,
            "GetDP",
            :getdp_executable,
            "failed to execute $executable -info: $(sprint(showerror, exception))";
            run_directory = run.path
        )
    end
    occursin("getdp", lowercase(output)) &&
    occursin(r"\d+\.\d+", output) || _fem_error(
        :getdp,
        "GetDP",
        :getdp_executable,
        "executable identity check did not report GetDP: $executable";
        run_directory = run.path
    )
    return executable
end

function _command_argument(value::AbstractString)
    return "\"" * replace(String(value), '\\' => "\\\\", '"' => "\\\"") * "\""
end

function _getdp_command(
        executable::String,
        model_path::String,
        mesh_path::String,
        verbosity::Int,
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        mesh_plan::FEMMeshPlan,
        frequency_index::Int,
        basis_terminal::Int
)
    frequency_index in eachindex(model.problem.frequencies) || throw(
        ArgumentError("frequency index $frequency_index is out of range")
    )
    basis_terminal in eachindex(model.terminal_ids) || throw(
        ArgumentError("basis terminal $basis_terminal is out of range")
    )
    job_name = @sprintf(
        "getdp-f%04d-b%04d", frequency_index, basis_terminal
    )
    arguments = [
        executable,
        model_path,
        "-solve",
        "LineCableModelsFEMScan",
        "-msh",
        mesh_path,
        "-name",
        joinpath(run.path, "input", job_name),
        "-v",
        string(verbosity),
        "-setstring",
        "RunDirectory",
        run.path,
        "-setstring",
        "RawOutputStem",
        job_name,
        "-setnumber",
        "FrequencyCount",
        string(length(model.problem.frequencies)),
        "-setnumber",
        "FrequencyIndex",
        string(frequency_index),
        "-setnumber",
        "FrequencyHz",
        _pro_number(mesh_plan.frequency),
        "-setnumber",
        "Val_Rint",
        _pro_number(mesh_plan.domain_radius),
        "-setnumber",
        "Val_Rext",
        _pro_number(mesh_plan.shell_outer_radius),
        "-setnumber",
        "BasisTerminal",
        string(basis_terminal),
        "-setnumber",
        "PlotFieldMaps",
        formulation.execution.plot_field_maps ? "1" : "0"
    ]
    return join(_command_argument.(arguments), " ")
end

function _job_raw_paths(run::FEMRun, job_name::String)
    root = joinpath(run.path, "raw", "jobs")
    return (
        Z = joinpath(root, "$job_name-Z.tsv"),
        P = joinpath(root, "$job_name-P.tsv")
    )
end

function _job_raw_data_row_count(path::String)
    isfile(path) || return 0
    return count(!isempty, strip.(readlines(path)))
end

function _read_job_raw_rows(
        path::String,
        terminal_count::Int,
        run::FEMRun
)
    isfile(path) || _fem_error(
        :getdp,
        basename(path),
        :raw_output,
        "GetDP job output is missing";
        run_directory = run.path
    )
    rows = readlines(path)
    length(rows) == terminal_count || _fem_error(
        :getdp,
        basename(path),
        :raw_output,
        "expected $terminal_count physical rows, found $(length(rows))";
        run_directory = run.path
    )
    for (line_index, row) in pairs(rows)
        isempty(strip(row)) && _fem_error(
            :getdp,
            basename(path),
            :raw_output,
            "unexpected blank row at line $line_index";
            run_directory = run.path
        )
        length(split(row, '\t'; keepempty = true)) == 6 || _fem_error(
            :getdp,
            basename(path),
            :raw_output,
            "line $line_index does not contain six tab-separated fields";
            run_directory = run.path
        )
    end
    return rows
end

function _append_job_raw!(
        aggregate_path::String,
        job_path::String,
        terminal_count::Int,
        run::FEMRun
)
    rows = _read_job_raw_rows(job_path, terminal_count, run)
    open(aggregate_path, "a") do io
        foreach(row -> println(io, row), rows)
    end
    return length(rows)
end

function _raw_data_row_count(path::String)
    isfile(path) || return 0
    lines = filter(!isempty, strip.(readlines(path)))
    return max(length(lines) - 1, 0)
end

function _write_scan_completion!(run::FEMRun, model::FEMResolvedModel)
    frequency_count = length(model.problem.frequencies)
    terminal_count = length(model.terminal_ids)
    expected_rows = frequency_count * terminal_count * terminal_count
    path = joinpath(run.path, "raw", "scan_complete.tsv")
    temporary = tempname(dirname(path))
    open(temporary, "w") do io
        println(io, join(FEM_COMPLETE_HEADER, '\t'))
        println(io, join((
            frequency_count,
            terminal_count,
            expected_rows,
            expected_rows,
            expected_rows,
            1
        ), '\t'))
    end
    mv(temporary, path; force = true)
    return path
end

function _persist_gmsh_log(path::String, records)
    mkpath(dirname(path))
    open(path, "w") do io
        for record in records
            println(io, record)
        end
    end
    return path
end

function _log_tail(records; count::Int = 20)
    isempty(records) && return "(no Gmsh log records)"
    first_index = max(firstindex(records), lastindex(records) - count + 1)
    return join(@view(records[first_index:lastindex(records)]), '\n')
end

function _run_getdp_unlocked!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        mesh_paths::AbstractVector{<:AbstractString}
)
    length(mesh_paths) == length(model.problem.frequencies) || throw(
        DimensionMismatch("one FEM mesh path is required per frequency")
    )
    executable = _resolve_getdp(formulation, run)
    assets = _getdp_assets()
    all(isfile, values(assets)) || _fem_error(
        :getdp,
        "LineCableModelsGmshExt",
        :assets,
        "one or more bundled GetDP formulation files are missing";
        run_directory = run.path
    )
    log_path = joinpath(run.path, "logs", "getdp.log")
    records = String[]
    terminal_count = length(model.terminal_ids)
    completed_jobs = 0
    z_path = joinpath(run.path, "raw", "Z.tsv")
    p_path = joinpath(run.path, "raw", "P.tsv")
    for frequency_index in eachindex(model.problem.frequencies)
        mesh_plan = model.mesh_plans[frequency_index]
        mesh_path = String(mesh_paths[frequency_index])
        for basis_terminal in eachindex(model.terminal_ids)
            job_name = @sprintf(
                "getdp-f%04d-b%04d", frequency_index, basis_terminal
            )
            job_raw = _job_raw_paths(run, job_name)
            foreach(path -> rm(path; force = true), job_raw)
            command = _getdp_command(
                executable,
                assets.model,
                mesh_path,
                formulation.execution.getdp_verbosity,
                run,
                model,
                formulation,
                mesh_plan,
                frequency_index,
                basis_terminal
            )
            client_name = @sprintf(
                "LineCableModelsFEMGetDP-f%04d-b%04d",
                frequency_index,
                basis_terminal
            )
            job_records = String[]
            deferred_exception = nothing
            run.getdp_invocations += 1
            gmsh.logger.start()
            try
                gmsh.onelab.run(client_name, command)
            catch exception
                try
                    append!(job_records, gmsh.logger.get())
                catch
                end
                z_rows = _job_raw_data_row_count(job_raw.Z)
                p_rows = _job_raw_data_row_count(job_raw.P)
                if z_rows == terminal_count && p_rows == terminal_count
                    deferred_exception = exception
                else
                    diagnostic_records = vcat(records, job_records)
                    _persist_gmsh_log(log_path, diagnostic_records)
                    _fem_error(
                        :getdp,
                        "GetDP",
                        :client,
                        "GetDP client failed for frequency index " *
                        "$frequency_index, basis terminal $basis_terminal: " *
                        "$(sprint(showerror, exception))\n" *
                        "raw row counts: Z=$z_rows, P=$p_rows, " *
                        "expected=$terminal_count\n" *
                        "Gmsh log tail:\n$(_log_tail(diagnostic_records))";
                        run_directory = run.path
                    )
                end
            finally
                try
                    isempty(job_records) && append!(job_records, gmsh.logger.get())
                catch
                end
                try
                    gmsh.logger.stop()
                catch
                end
                isempty(job_records) || append!(records, job_records)
                _persist_gmsh_log(log_path, records)
            end
            _read_job_raw_rows(job_raw.Z, terminal_count, run)
            _read_job_raw_rows(job_raw.P, terminal_count, run)
            _append_job_raw!(z_path, job_raw.Z, terminal_count, run)
            _append_job_raw!(p_path, job_raw.P, terminal_count, run)
            completed_jobs += 1
            expected_job_rows = completed_jobs * terminal_count
            z_rows = _raw_data_row_count(z_path)
            p_rows = _raw_data_row_count(p_path)
            z_rows == expected_job_rows && p_rows == expected_job_rows || _fem_error(
                :getdp,
                "GetDP",
                :raw_output,
                "GetDP job for frequency index $frequency_index, basis " *
                "terminal $basis_terminal returned without one complete " *
                "response column; raw row counts are Z=$z_rows, P=$p_rows, " *
                "expected=$expected_job_rows";
                run_directory = run.path
            )
            deferred_exception === nothing || @debug(
                "Gmsh reported a client error after GetDP wrote a complete " *
                "response column; deferring success to strict result validation",
                frequency_index,
                basis_terminal,
                exception = sprint(showerror, deferred_exception)
            )
        end
    end
    _write_scan_completion!(run, model)
    return nothing
end

function _getdp_socket_name(run::FEMRun)
    Sys.iswindows() && return "127.0.0.1:0"
    digest = bytes2hex(sha256(codeunits(run.path)))[1:20]
    return ".linecablemodels-gmsh-$digest.sock"
end

function _run_getdp!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        mesh_paths::AbstractVector{<:AbstractString}
)
    previous_socket = gmsh.option.get_string("Solver.SocketName")
    socket_name = _getdp_socket_name(run)
    gmsh.option.set_string("Solver.SocketName", socket_name)
    @info "Assigned an isolated Gmsh/GetDP ONELAB socket" socket_name
    try
        return _run_getdp_unlocked!(run, model, formulation, mesh_paths)
    finally
        gmsh.option.set_string("Solver.SocketName", previous_socket)
    end
end

function _run_getdp!(
        run::FEMRun,
        model::FEMResolvedModel,
        formulation::LineCableModelsFEM,
        mesh_path::String
)
    length(model.problem.frequencies) == 1 || throw(DimensionMismatch(
        "a single FEM mesh path is valid only for a one-frequency problem"
    ))
    return _run_getdp!(run, model, formulation, [mesh_path])
end
