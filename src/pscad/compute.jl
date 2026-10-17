function computation_options(
        ::Union{Type{PSCADFormulation}, Type{PSCADTransientFormulation}},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    allowed = (:output_stem, :remote, :verbosity, :output_basis, :on_result,
        :resume_run_directory, :solver_identity, :work_root, :timing)
    unknown = filter(key -> key ∉ allowed, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown PSCAD computation options: $(sort!(collect(unknown)))",
    ))
    haskey(options, :remote) || throw(ArgumentError(
        "PSCAD computation options must define remote::RemoteConfig",
    ))
    options.remote isa RemoteConfig || throw(ArgumentError(
        "PSCAD computation option remote must be a RemoteConfig",
    ))
    normalized = merge(
        (
            output_stem = "lcm",
            work_root = options.remote.local_root,
            verbosity = (default = 0,),
            output_basis = :pul,
            on_result = nothing,
            resume_run_directory = nothing,
            solver_identity = nothing,
            timing = false
        ),
        options
    )
    normalized.output_stem isa AbstractString || throw(ArgumentError(
        "PSCAD output_stem must be a string",
    ))
    output_stem = String(normalized.output_stem)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_]{0,19}$", output_stem) ||
        throw(ArgumentError(
            "PSCAD output_stem must contain 1–20 ASCII letters, digits, or underscores",
        ))
    levels = verbosity(normalized.verbosity)
    normalized.timing isa Bool || throw(ArgumentError("timing must be Bool"))
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    resume = normalized.resume_run_directory
    (resume === nothing || resume === :latest || resume isa AbstractString) ||
        throw(ArgumentError("PSCAD resume_run_directory must be nothing, :latest, or a completed run path"))
    resume isa AbstractString && isempty(resume) &&
        throw(ArgumentError(
            "PSCAD resume_run_directory cannot be empty"))
    identity = normalized.solver_identity
    (identity === nothing || identity isa AbstractDict{String, String}) ||
        throw(ArgumentError("PSCAD solver_identity must be the record returned by identify(remote)"))
    identity = identity === nothing ? nothing :
               Dict{String, String}(validate(identity, options.remote))
    options.remote.transport === :local && !Sys.iswindows() && throw(ArgumentError(
        "PSCAD local transport requires a Windows caller"))
    work_root = abspath(normalized.work_root)
    first(splitpath(relpath(work_root, options.remote.local_root))) == ".." &&
        throw(ArgumentError("PSCAD work_root must be inside remote.local_root"))
    return ComputationOptions(;
        work_root,
        output_stem,
        remote = options.remote,
        verbosity = levels,
        timing = normalized.timing,
        output_basis = Val(basis_value),
        on_result = normalized.on_result,
        resume_run_directory = resume isa AbstractString ? abspath(resume) : resume,
        solver_identity = identity
    )
end

function _pscad_size(problem::LineParametersProblem)
    assignments = problem.system.connection_order
    isempty(assignments) && throw(ArgumentError(
        "PSCAD computation requires at least one explicit terminal",
    ))
    any(iszero, assignments) && throw(ArgumentError(
        "PSCAD computation does not permit conductor elimination",
    ))
    all(>(0), assignments) || throw(ArgumentError(
        "PSCAD computation phase assignments must identify active phases",
    ))
    length(unique(assignments)) == length(assignments) || throw(ArgumentError(
        "PSCAD computation does not permit bundled terminals",
    ))
    sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
        "PSCAD computation active-phase assignments must be contiguous from 1",
    ))
    return (length(assignments), length(assignments), length(problem.frequencies))
end

function _pscad_basis(parameters, ::LineParametersProblem, ::Val{:pul})
    parameters
end

function _pscad_basis(
        parameters::LineParameters,
        problem::LineParametersProblem,
        ::Val{:total}
)
    return LineParameters(
        PhaseDomain,
        Z(parameters) .* problem.system.line_length,
        Y(parameters) .* problem.system.line_length,
        frequencies(parameters);
        basis = :total, details = details(parameters)
    )
end

function _pscad_inputs(problem::LineParametersProblem, formulation::PSCADFormulation, blueprints)
    _pscad_deterministic(eltype(problem), typeof(formulation.options.data.base_frequency))
    setting = pscad_setting(formulation, problem, blueprints)
    frequency = formulation.options.data.base_frequency
    temperature = problem.temperature
    # PSCAD writes one cable for each coaxial group and numbers its conductors cable by
    # cable, in the order of the blueprints.
    groups = [(blueprint, conductors) for blueprint in blueprints
              for conductors in blueprint.assembly_ranges]
    native_order = [(cable = design, terminal = conductor.terminal)
                    for (design, blueprint) in enumerate(blueprints)
                    for conductor in blueprint.conductors]
    native_cables = [cable for (cable, (_, conductors)) in enumerate(groups)
                     for _ in conductors]
    allunique(native_order) && Set(native_order) == Set(problem.system.terminal_order) ||
        throw(ArgumentError("PSCAD equivalent components must represent every terminal exactly once"))
    requested_order = problem.system.terminal_order[sortperm(problem.system.connection_order)]
    permutation = [only(findall(==(terminal), native_order)) for terminal in requested_order]
    coordinates = terminal_labels((cable = native_cables[index],
        terminal = native_order[index].terminal) for index in permutation)
    document = _pscad_project(problem.system, problem.earth_props, frequency, blueprints,
        formulation, temperature; native_settings = setting[(:ground, :frequency)])
    project = sprint(print, document)
    # The loss tangent G/(ωC) of each group, as the project writer evaluates it.
    dielectric_losses = map(enumerate(groups)) do (cable, (blueprint, conductors))
        map(enumerate(conductors)) do (layer, index)
            requested = Float64(_pscad_part_parameters(blueprint, index, layer, frequency,
                formulation, temperature).loss)
            exported = parse(Float64, _pscad_value(requested; maximum = 10))
            (; cable, layer, requested, exported, capped = requested > 10)
        end
    end
    return (; setting, project, native_order, permutation, coordinates, dielectric_losses)
end

function _stage_pscad_project(inputs, work_root::AbstractString)
    root = mktempdir(mkpath(abspath(work_root)); prefix = "run-", cleanup = false)
    @debug "Exporting PSCAD computation project"
    staged = joinpath(root, "generated.pscx")
    write(staged, inputs.project)
    return (; root, staged)
end

const PSCAD_COMPLETED_OUTPUTS = ("result_zm.out", "result_zp.out", "result_ym.out",
    "result_yp.out", "solver.toml", "native-settings.toml", "timing.txt")

_pscad_digest(record::AbstractDict) =
    bytes2hex(sha256(sprint(io -> TOML.print(io, record; sorted = true))))

function _compute_pscad(problem::LineParametersProblem, formulation::PSCADFormulation,
        execution_options, inputs)
    setting = inputs.setting
    input = Dict{String, Any}(
        "schema_version" => 4,
        "project_sha256" => bytes2hex(sha256(inputs.project)),
        "frequencies" => Float64.(problem.frequencies),
        "matrix_size" => collect(_pscad_size(problem)),
        "native_settings" => Dict(string(component) => Dict(string(field) => Dict(
            "value" => control.value, "readback" => control.readback)
            for (field, control) in pairs(getproperty(setting, component)))
            for component in (:ground, :frequency, :configuration)),
        "toolkit" => Dict(name => bytes2hex(sha256(source)) for (name, source) in PSCAD_REMOTE_SOURCES))
    @debug "Computing PSCAD line parameters" system=problem.system.system_id
    return _compute_pscad(problem, formulation, execution_options, inputs, input,
        PSCAD_COMPLETED_OUTPUTS)
end

function _compute_pscad(problem, formulation, execution_options, inputs, input, completed_outputs)
    config = execution_options.data.remote
    # Expected identity constrains execution, not numerical compatibility. The
    # completion record independently authenticates the complete request.
    signature = _pscad_digest(input)
    expected = execution_options.data.solver_identity
    expected === nothing || (input["expected_solver"] = expected)
    resume = execution_options.data.resume_run_directory
    parent = execution_options.data.work_root
    candidates = resume === nothing ? String[] : resume === :latest ?
        (isdir(parent) ? sort!(filter(path -> isdir(path) && isfile(joinpath(path, "complete.toml")),
            readdir(parent; join = true)); by = path -> mtime(joinpath(path, "complete.toml")), rev = true) : String[]) : [resume]
    source_root = nothing
    station_identity = nothing
    actual = nothing
    for candidate in candidates
        record_path = joinpath(candidate, "complete.toml")
        isfile(record_path) || throw(ArgumentError("PSCAD run has no completion record: $candidate"))
        record = TOML.parsefile(record_path)
        if get(record, "schema_version", nothing) != 4
            resume === :latest && continue
            throw(ArgumentError("PSCAD completed run requires a fresh version-4 computation: $candidate"))
        end
        if get(record, "input_sha256", nothing) != signature
            resume === :latest && continue
            throw(ArgumentError("PSCAD completed run has different numerical inputs or solver implementation: $candidate"))
        end
        stored = TOML.parsefile(joinpath(candidate, "computation.toml"))
        _pscad_digest(stored) == get(record, "request_sha256", nothing) || throw(ArgumentError(
            "PSCAD completed-run input integrity check failed: $candidate"))
        numerical = copy(stored)
        pop!(numerical, "expected_solver", nothing)
        _pscad_digest(numerical) == signature || throw(ArgumentError(
            "PSCAD completed-run numerical input integrity check failed: $candidate"))
        bytes2hex(open(sha256, joinpath(candidate, "generated.pscx"))) == stored["project_sha256"] ||
            throw(ArgumentError("PSCAD completed-run exported project changed: $candidate"))
        for (name, digest) in stored["toolkit"]
            bytes2hex(open(sha256, joinpath(candidate, "toolkit", name))) == digest ||
                throw(ArgumentError("PSCAD completed-run solver source changed: $candidate/toolkit/$name"))
        end
        for name in completed_outputs
            path = joinpath(candidate, "outputs", name)
            isfile(path) && bytes2hex(open(sha256, path)) == get(record["outputs"], name, nothing) ||
                throw(ArgumentError("PSCAD completed-run output integrity check failed: $path"))
        end
        stored_solver_identity = Dict{String, String}(validate(
            TOML.parsefile(joinpath(candidate, "outputs", "solver.toml")), config))
        stored_solver_identity == get(record, "observed_solver", nothing) || throw(ArgumentError(
            "PSCAD completed run has no matching solver attestation: $candidate"))
        station_identity === nothing && (station_identity = identify(config))
        expected === nothing || expected == station_identity || throw(ArgumentError(
            "PSCAD station does not match the expected solver identity"))
        if stored_solver_identity != station_identity
            resume === :latest && continue
            throw(ArgumentError("PSCAD completed run belongs to a different solver installation: $candidate"))
        end
        source_root, actual = candidate, stored_solver_identity
        break
    end
    reused = source_root !== nothing
    if reused
        @debug "PSCAD reuses a verified completed run" source_run=source_root
        output = joinpath(source_root, "outputs")
        execution = (exit_code = 0,
            stdout_path = joinpath(output, "stdout.txt"), stderr_path = joinpath(output, "stderr.txt"),
            console_path = joinpath(output, "pscad-console.txt"), output_dir = output)
    else
        root, staged = _stage_pscad_project(inputs, parent)
        source_root = root
        open(joinpath(root, "computation.toml"), "w") do io
            TOML.print(io, input; sorted = true)
        end
        execution = run_remote_pscad(config, staged, joinpath(root, "outputs"),
            formulation, input["frequencies"]; output_stem = execution_options.data.output_stem,
            verbosity = verbosity(execution_options, :PSCAD))
        actual = Dict{String, String}(validate(
            TOML.parsefile(joinpath(execution.output_dir, "solver.toml")), config))
        expected === nothing || actual == expected || throw(ArgumentError(
            "PSCAD result does not match the expected solver identity: $root"))
    end
    native_readback = TOML.parsefile(joinpath(execution.output_dir, "native-settings.toml"))
    validate(native_readback, formulation, input, execution)
    result = try
        read_pscad_result(execution.output_dir, problem, formulation;
            output_stem = execution_options.data.output_stem)
    catch error
        verbosity(execution_options, :progress) > 0 && @info "PSCAD result validation failed" _group=:progress exception=error
        throw(ErrorException("PSCAD result validation failed: $(sprint(showerror, error))" *
            "\nLast PSCAD diagnostics:\n$(_diagnostic_tail(execution.console_path))" *
            "\nFull PSCAD diagnostics: $(execution.console_path)"))
    end
    source_elapsed = parse(Float64, strip(read(joinpath(execution.output_dir, "timing.txt"), String)))
    isfinite(source_elapsed) && source_elapsed >= 0 || throw(ArgumentError("invalid PSCAD execution timing"))
    if !reused
        completion = Dict("schema_version" => 4, "input_sha256" => signature,
            "request_sha256" => _pscad_digest(input), "observed_solver" => actual,
            "outputs" => Dict(name => bytes2hex(open(sha256, joinpath(execution.output_dir, name)))
                for name in completed_outputs))
        temporary = tempname(source_root)
        try
            open(temporary, "w") do io
                TOML.print(io, completion; sorted = true)
            end
            mv(temporary, joinpath(source_root, "complete.toml"); force = false)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
    files = [(path = relpath(joinpath(directory, name), source_root),
        source = joinpath(directory, name), sha256 = bytes2hex(open(sha256, joinpath(directory, name))))
        for (directory, _, names) in walkdir(source_root) for name in sort(names)]
    # Reused and fresh native records have the same execution-fact schema.
    # The public remote call's duration is transient. The scan projects it once.
    execution = (; (key => value for (key, value) in pairs(execution)
        if key ∉ (:elapsed_seconds, :elapsed_scope))...)
    execution = merge(execution, (backend = :pscad, pscad_version = config.pscad_version,
        reused, source_run = source_root, input_sha256 = signature, solver_identity = actual))
    return (; result, native_readback, files, execution, source_elapsed)
end

function _pscad_result(problem, formulation, options, inputs, native; batch_reuse = false)
    permutation = inputs.permutation
    parameters = LineParameters(PhaseDomain,
        Z(native.result)[permutation, permutation, :], Y(native.result)[permutation, permutation, :],
        copy(frequencies(native.result)); details = details(native.result))
    parameters = _pscad_basis(parameters, problem, options.data.output_basis)
    all(isfinite, Z(parameters)) || throw(DomainError(Z(parameters),
        "completed PSCAD series impedance must contain only finite entries"))
    all(isfinite, Y(parameters)) || throw(DomainError(Y(parameters),
        "completed PSCAD shunt admittance must contain only finite entries"))
    execution = native.execution
    batch_reuse && (execution=merge(execution, (reused=true,)))
    retained = (files = deepcopy(native.files), coordinates = copy(inputs.coordinates),
        native_terminal_order = copy(inputs.native_order), requested_frequencies = copy(problem.frequencies),
        formulations = computation_details(formulation).data, native_setting = inputs.setting,
        base_frequency = formulation.options.data.base_frequency, loss_tangent_limit = 10.0,
        aerial_shunt_conductance = 1e-38, native_readback = deepcopy(native.native_readback),
        native_frequencies = copy(details(native.result).data.native_frequencies),
        dielectric_losses = deepcopy(inputs.dielectric_losses), exported_project = inputs.project,
        execution = deepcopy(execution))
    if options.data.timing
        timing = execution.reused ? (;) : (compile_call_seconds=native.source_elapsed,)
        retained = merge(retained, (; timing))
    end
    return LineParameters(parameters.domain, parameters.Z, parameters.Y, parameters.f, Engine.completion_details(retained))
end

function compute(problem::LineParametersProblem, formulation::PSCADFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions(),
        modal=nothing, modal_options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
    modal===nothing || return compute(problem,formulation,
        LineCableModels.ModalAnalysisFormulation(modal);options,modal_options)
    isempty(ComputationOptions(modal_options).data) ||
        throw(ArgumentError("modal_options require a modal formulation"))
    return first(compute(problem, [formulation]; options))
end

function compute(problem::LineParametersProblem, formulations::AbstractVector{<:PSCADFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions(),
        modal=nothing, modal_options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
    isempty(formulations) && throw(ArgumentError("PSCAD formulation collections cannot be empty"))
    options = ComputationOptions(options)
    modal===nothing || return compute(problem,formulations,
        LineCableModels.ModalAnalysisFormulation(modal);options,modal_options)
    isempty(ComputationOptions(modal_options).data) ||
        throw(ArgumentError("modal_options require a modal formulation"))
    execution = computation_options(PSCADFormulation, options)
    _pscad_deterministic(eltype(problem))
    validate(problem.frequencies, PSCADFormulation)
    _pscad_size(problem)
    blueprints = _pscad_blueprints(problem.system)
    inputs = [_pscad_inputs(problem, formulation, blueprints) for formulation in formulations]
    logger = LineCableModels.VerbosityLogger(Logging.current_logger(), execution.data.verbosity)
    return Logging.with_logger(logger) do
        _compute_pscad(problem, formulations, execution, inputs)
    end
end

function _compute_pscad(problem::LineParametersProblem,
        formulations::AbstractVector{<:PSCADFormulation}, execution::ComputationOptions, inputs)
    progress = verbosity(execution, :progress) > 0
    started = progress ? time_ns() : UInt64(0)
    last_log = started
    previous = started
    average_seconds = 0.0
    progress && @info "PSCAD computation started" _group=:progress total=length(formulations)
    physical_inputs = Engine.completed_inputs(problem)
    source_id = gridpoint_id().source_id
    keys = [(project = value.project, setting = value.setting[(:ground, :frequency, :configuration)])
            for value in inputs]
    # Include native execution, readback, and final result construction, but
    # exclude batch input construction, measurement attachment, and callbacks.
    scan_started = execution.data.timing ? time_ns() : UInt64(0)
    native = _compute_pscad(problem, first(formulations), execution, first(inputs))
    execution = ComputationOptions(merge(execution.data, (solver_identity = native.execution.solver_identity,)))
    first_result = Engine.retain_gridpoint(_pscad_result(problem, first(formulations), execution,
        first(inputs), native), gridpoint_id(; source_id);
        fields = merge(Engine.completed_formulation(first(formulations)), (inputs = physical_inputs,)))
    if execution.data.timing && !isempty(first_result.details.data.timing)
        wall_seconds = (time_ns() - scan_started) * 1e-9
        first_result = Engine.retain_gridpoint(first_result, first_result.details.data.gridpoint;
            fields=(timing=merge((; wall_seconds), first_result.details.data.timing),))
    end
    values = Vector{typeof(first_result)}(undef, length(formulations))
    values[1] = first_result
    execution.data.on_result === nothing || execution.data.on_result(problem, 1, first_result)
    if progress
        now = time_ns()
        average_seconds = (now - previous) * 1e-9
        previous = now
        if now - last_log >= 5_000_000_000
            @info "PSCAD progress" _group=:progress completed=1 total=length(formulations) elapsed_seconds=(now-started)*1e-9 eta_hours=(length(formulations)-1)*average_seconds/3600
            last_log = now
        end
    end
    completed = Dict(first(keys) => native)
    for index in 2:length(formulations)
        shared = haskey(completed, keys[index])
        scan_started = execution.data.timing ? time_ns() : UInt64(0)
        native = if shared
            @debug "PSCAD reuses identical exported inputs" formulation=index
            completed[keys[index]]
        else
            _compute_pscad(problem, formulations[index], execution, inputs[index])
        end
        value = _pscad_result(problem, formulations[index], execution, inputs[index], native; batch_reuse = shared)
        value = Engine.retain_gridpoint(value,
            gridpoint_id(; source_id, formulation_index = index);
            fields = merge(Engine.completed_formulation(formulations[index]), (inputs = physical_inputs,)))
        if execution.data.timing && !isempty(value.details.data.timing)
            wall_seconds = (time_ns() - scan_started) * 1e-9
            value = Engine.retain_gridpoint(value, value.details.data.gridpoint;
                fields=(timing=merge((; wall_seconds), value.details.data.timing),))
        end
        typeof(value) === eltype(values) || throw(ArgumentError("PSCAD formulations produced inconsistent result types"))
        values[index] = value
        completed[keys[index]] = native
        execution.data.on_result === nothing || execution.data.on_result(problem, index, values[index])
        if progress
            now = time_ns()
            interval = (now - previous) * 1e-9
            average_seconds = 0.2 * interval + 0.8 * average_seconds
            previous = now
            if now - last_log >= 5_000_000_000
                @info "PSCAD progress" _group=:progress completed=index total=length(formulations) elapsed_seconds=(now-started)*1e-9 eta_hours=(length(formulations)-index)*average_seconds/3600
                last_log = now
            end
        end
    end
    progress && @info "PSCAD computation completed successfully" _group=:progress completed=length(values) total=length(formulations) elapsed_seconds=(time_ns()-started)*1e-9
    return values
end

function validate(readback::AbstractDict, ::PSCADFormulation, input::AbstractDict, execution)
    expected = Dict(component => Dict(field => control["readback"]
        for (field, control) in controls) for (component, controls) in input["native_settings"])
    readback == expected || throw(ArgumentError(
        "PSCAD result has no matching native-settings readback: $(execution.output_dir)"))
    return readback
end

function validate(readback::AbstractDict, ::PSCADTransientFormulation, input::AbstractDict, execution)
    expected = Dict(component => Dict(field => control["readback"]
        for (field, control) in controls) for (component, controls) in input["native_settings"])
    Set(keys(readback)) == Set((keys(expected)..., "project")) ||
        throw(ArgumentError("PSCAD transient result has incomplete settings readback"))
    all(readback[name] == controls for (name, controls) in expected) ||
        throw(ArgumentError("PSCAD result has no matching native-settings readback: $(execution.output_dir)"))
    project = readback["project"]
    Set(keys(project)) == Set(keys(input["project_settings"])) ||
        throw(ArgumentError("PSCAD transient project settings readback is incomplete"))
    for field in ("time_duration", "time_step", "sample_step", "PlotType")
        parse(Float64, string(project[field])) == input["project_settings"][field] ||
            throw(ArgumentError("PSCAD transient project setting $field differs from the request"))
    end
    project["output_filename"] == input["project_settings"]["output_filename"] ||
        throw(ArgumentError("PSCAD transient output filename differs from the request"))
    return readback
end

function read_pscad_result(directory::AbstractString, problem::LineParametersProblem,
        ::PSCADFormulation; output_stem)
    dimensions = _pscad_size(problem)
    result = read_pscad_result(directory, problem.frequencies, dimensions)
    all(isfinite, Z(result)) || throw(DomainError(Z(result), "PSCAD series impedance must be finite"))
    all(isfinite, Y(result)) || throw(DomainError(Y(result), "PSCAD shunt admittance must be finite"))
    return result
end

function read_pscad_result(directory::AbstractString,
        problem::TransientAnalysisProblem{<:LineParametersProblem},
        ::PSCADTransientFormulation; output_stem)
    execution = TOML.parsefile(joinpath(directory, "execution.toml"))
    execution == Dict("mode" => "transient", "stem" => output_stem) ||
        throw(ArgumentError("PSCAD transient execution record differs from the request"))
    record = read_pscad_transient(directory; stem = output_stem)
    return read_pscad_transient(record, problem)
end

function _pscad_inputs(problem::TransientAnalysisProblem{<:LineParametersProblem},
        formulation::PSCADTransientFormulation, work_root::AbstractString)
    return mktempdir(mkpath(work_root)) do directory
        path = export_data(:pscad, problem, formulation; file_name = joinpath(directory, "transient.pscx"))
        document = readxml(path)
        frequency = only(findall("//User[@defn='master:Line_FrePhase_Options' or @defn='master:Line_FreDep_Options']", document))
        numeric_fields = ("FS", "FE", "Numf", "FT", "YMaxP", "AMaxP", "YMaxE", "AMaxE",
            "ZMaxP", "ZMaxE", "MaxRPtol", "W1", "W2", "W3")
        frequency_controls = Dict(node["name"] => Dict("value" => parse(Float64, node["value"]),
            "readback" => parse(Float64, node["value"]))
            for node in findall("./paramlist/param", frequency) if node["name"] in numeric_fields)
        ground = only(findall("//User[@defn='master:Line_Ground']", document))
        ground_controls = Dict(node["name"] => Dict("value" => parse(Float64, node["value"]),
            "readback" => parse(Float64, node["value"]))
            for node in findall("./paramlist/param", ground)
            if node["name"] in ("GRRES", "GPERM", "GRP", "K0", "K1", "alpha"))
        native_settings = Dict("ground" => ground_controls, "frequency" => frequency_controls,
            "configuration" => Dict("Freq" => Dict("value" => Float64(formulation.frequency),
                "readback" => Float64(formulation.frequency))))
        project_settings = Dict{String, Any}(node["name"] => parse(Float64, node["value"])
            for node in findall("/project/paramlist/param", document)
            if node["name"] in ("time_duration", "time_step", "sample_step", "PlotType"))
        project_settings["output_filename"] = formulation.output_stem * ".out"
        channels = [only(findall("./paramlist/param[@name='Name']", node))["value"]
            for node in findall("//User[@defn='master:pgb']", document)]
        return (; project = read(path, String), native_settings, project_settings, channels)
    end
end

function _compute_pscad(problem::TransientAnalysisProblem{<:LineParametersProblem},
        formulation::PSCADTransientFormulation, execution_options, inputs)
    input = Dict{String, Any}("schema_version" => 4, "mode" => "transient",
        "project_sha256" => bytes2hex(sha256(inputs.project)),
        "frequencies" => Float64.(problem.line.frequencies),
        "native_settings" => inputs.native_settings, "project_settings" => inputs.project_settings,
        "toolkit" => Dict(name => bytes2hex(sha256(source)) for (name, source) in PSCAD_REMOTE_SOURCES))
    stem = execution_options.data.output_stem
    completed_outputs = ["solver.toml", "native-settings.toml", "timing.txt", "execution.toml",
        stem * ".inf", [stem * "_" * lpad(string(index), 2, '0') * ".out"
            for index in 1:cld(length(inputs.channels), 10)]...]
    @debug "Computing PSCAD transient response" system=problem.line.system.system_id
    return _compute_pscad(problem, formulation, execution_options, inputs, input, completed_outputs)
end

function compute(problem::TransientAnalysisProblem{<:LineParametersProblem},
        formulation::PSCADTransientFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    execution = computation_options(PSCADTransientFormulation,
        ComputationOptions(merge((output_stem = formulation.output_stem,), ComputationOptions(options).data)))
    formulation = PSCADTransientFormulation(; formulation.model, formulation.fit,
        formulation.frequency, formulation.sample, output_stem = execution.data.output_stem)
    validate(problem, PSCADTransientFormulation)
    inputs = _pscad_inputs(problem, formulation, execution.data.work_root)
    logger = LineCableModels.VerbosityLogger(Logging.current_logger(), execution.data.verbosity)
    return Logging.with_logger(logger) do
        started = time_ns()
        native = _compute_pscad(problem, formulation, execution, inputs)
        response = native.result
        retained = merge(details(response).data, (formulations = computation_details(formulation).data,
            window = problem.window, gridpoint = gridpoint_id(), coordinates = copy(response.coordinates),
            files = native.files, exported_project = inputs.project, execution = native.execution,
            native_readback = native.native_readback))
        if execution.data.timing
            timing = native.execution.reused ? (;) :
                (wall_seconds = (time_ns() - started) * 1e-9, run_call_seconds = native.source_elapsed)
            retained = merge(retained, (; timing))
        end
        result = Engine.TransientAnalysis.TransientResponse(response.t, response.voltage,
            response.current, response.coordinates, ComputationDetails(retained))
        execution.data.on_result === nothing || execution.data.on_result(problem, 1, result)
        return result
    end
end

function computation_details(formulation::PSCADTransientFormulation)::ComputationDetails
    return ComputationDetails(; backend = :pscad, model = formulation.model, fit = formulation.fit,
        frequency = formulation.frequency, sample = formulation.sample, output_stem = formulation.output_stem)
end

function computation_details(::Type{<:PSCADTransientFormulation},
        result::Engine.TransientAnalysis.TransientResponse)::ComputationDetails
    return details(result)
end
