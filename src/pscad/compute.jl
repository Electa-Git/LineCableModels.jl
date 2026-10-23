# A normalized record is returned as it is.
computation_options(::Type{<:PSCADFormulation},
    record::ComputationOptions{<:NamedTuple{(:work_root, :remote,
        :resume_run_directory, :solver_identity),
        <:Tuple{String, RemoteConfig, Any, Any}}}) = record

function computation_options(
        ::Type{<:PSCADFormulation},
        record::ComputationOptions
)::ComputationOptions
    Base.get_extension(LineCableModels, :LineCableModelsPSCADExt) === nothing &&
        throw(ArgumentError("PSCAD native runs require `using PythonCall`"))
    options = record.data
    allowed = (:remote, :resume_run_directory, :solver_identity, :work_root)
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
            work_root = options.remote.local_root,
            resume_run_directory = nothing,
            solver_identity = nothing
        ),
        options
    )
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
    work_root = abspath(normalized.work_root)
    first(splitpath(relpath(work_root, options.remote.local_root))) == ".." &&
        throw(ArgumentError("PSCAD work_root must be inside remote.local_root"))
    return ComputationOptions(;
        work_root,
        remote = options.remote,
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
    terminals = ["cable:$(native_cables[index]):$(native_order[index].terminal)"
                 for index in permutation]
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
    return (; setting, project, native_order, permutation, terminals, dielectric_losses)
end

# The run directory contains the project as `<namespace>.pscx`, since PSCAD names a loaded
# project after its file. The run's line takes the name of the directory, and the detailed
# files that PSCAD writes into the shared directory carry that name.
function _stage_pscad_project(inputs, work_root::AbstractString, system::LineCableSystem)
    root = mktempdir(mkpath(abspath(work_root)); prefix = "run-", cleanup = false)
    @debug "Exporting PSCAD computation project"
    staged = joinpath(root, _pscad_namespace(system) * ".pscx")
    write(staged, inputs.project)
    line_name = "lcm_" * chop(basename(root); head = 4, tail = 0)
    return (; root, staged, line_name)
end

const PSCAD_COMPLETED_OUTPUTS = ("result_zm.out", "result_zp.out", "result_ym.out",
    "result_yp.out", "solver.toml", "native-settings.toml", "timing.txt")

_pscad_digest(record::AbstractDict) =
    bytes2hex(sha256(sprint(io -> TOML.print(io, record; sorted = true))))

function _compute_pscad(problem::LineParametersProblem, formulation::PSCADFormulation,
        execution_options, inputs)
    config = execution_options.data.remote
    setting = inputs.setting
    input = Dict{String, Any}(
        "schema_version" => 5,
        "project_sha256" => bytes2hex(sha256(inputs.project)),
        "frequencies" => Float64.(problem.frequencies),
        "matrix_size" => collect(_pscad_size(problem)),
        "native_settings" => Dict(string(component) => Dict(string(field) => Dict(
            "value" => control.value, "readback" => control.readback)
            for (field, control) in pairs(getproperty(setting, component)))
            for component in (:ground, :frequency, :configuration)))
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
        if get(record, "schema_version", nothing) != 5
            resume === :latest && continue
            throw(ArgumentError("PSCAD completed run requires a fresh version-5 computation: $candidate"))
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
        for name in PSCAD_COMPLETED_OUTPUTS
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
        execution = (console_path = joinpath(output, "pscad-console.txt"), output_dir = output)
    else
        root, staged, line_name = _stage_pscad_project(inputs, parent, problem.system)
        source_root = root
        open(joinpath(root, "computation.toml"), "w") do io
            TOML.print(io, input; sorted = true)
        end
        @debug "Computing PSCAD line parameters" system=problem.system.system_id
        execution = run_remote_pscad(config, staged, joinpath(root, "outputs"),
            input["native_settings"], problem.frequencies; line_name, expected_identity = expected)
        actual = Dict{String, String}(validate(
            TOML.parsefile(joinpath(execution.output_dir, "solver.toml")), config))
        expected === nothing || actual == expected || throw(ArgumentError(
            "PSCAD result does not match the expected solver identity: $root"))
    end
    native_readback = TOML.parsefile(joinpath(execution.output_dir, "native-settings.toml"))
    expected_readback = Dict(component => Dict(field => control["readback"]
        for (field, control) in controls) for (component, controls) in input["native_settings"])
    native_readback == expected_readback || throw(ArgumentError(
        "PSCAD result has no matching native-settings readback: $(execution.output_dir)"))
    parameters = try
        read_pscad_result(execution.output_dir, problem.frequencies, _pscad_size(problem),
            problem.system.line_length)
    catch error
        throw(ErrorException("PSCAD result validation failed: $(sprint(showerror, error))" *
            "\nLast PSCAD diagnostics:\n$(_diagnostic_tail(execution.console_path))" *
            "\nFull PSCAD diagnostics: $(execution.console_path)"))
    end
    source_elapsed = parse(Float64, strip(read(joinpath(execution.output_dir, "timing.txt"), String)))
    isfinite(source_elapsed) && source_elapsed >= 0 || throw(ArgumentError("invalid PSCAD execution timing"))
    if !reused
        completion = Dict("schema_version" => 5, "input_sha256" => signature,
            "request_sha256" => _pscad_digest(input), "observed_solver" => actual,
            "outputs" => Dict(name => bytes2hex(open(sha256, joinpath(execution.output_dir, name)))
                for name in PSCAD_COMPLETED_OUTPUTS))
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
    execution = merge(execution, (backend = :pscad, pscad_version = PSCAD_VERSION,
        reused, source_run = source_root, input_sha256 = signature, solver_identity = actual,
        compile_call_seconds = source_elapsed))
    return (; parameters, native_readback, files, execution)
end

function _pscad_result(problem, formulation, inputs, native; batch_reuse = false)
    permutation = inputs.permutation
    parameters = LineParameters(PhaseDomain,
        Z(native.parameters)[permutation, permutation, :], Y(native.parameters)[permutation, permutation, :],
        copy(frequencies(native.parameters)); terminals = copy(inputs.terminals),
        line_length = problem.system.line_length, details = details(native.parameters))
    all(isfinite, Z(parameters)) || throw(DomainError(Z(parameters),
        "completed PSCAD series impedance must contain only finite entries"))
    all(isfinite, Y(parameters)) || throw(DomainError(Y(parameters),
        "completed PSCAD shunt admittance must contain only finite entries"))
    execution = native.execution
    batch_reuse && (execution=merge(execution, (reused=true,)))
    retained = (files = deepcopy(native.files),
        native_terminal_order = copy(inputs.native_order), requested_frequencies = copy(problem.frequencies),
        native_setting = inputs.setting,
        base_frequency = formulation.options.data.base_frequency, loss_tangent_limit = 10.0,
        aerial_shunt_conductance = 1e-38, native_readback = deepcopy(native.native_readback),
        native_frequencies = copy(details(native.parameters).data.native_frequencies),
        dielectric_losses = deepcopy(inputs.dielectric_losses), exported_project = inputs.project,
        execution = deepcopy(execution))
    return LineParameters(parameters.domain, parameters.Z, parameters.Y, parameters.f,
        parameters.terminals, parameters.line_length, ComputationDetails(retained))
end

"""
$(TYPEDSIGNATURES)

Resolve the parts of a line-parameter problem that its PSCAD formulations share. Once per
problem, it checks the problem and builds the blueprints and the completed inputs. It
exports the inputs of each formulation. The formulations share one table of native runs by
the content keys of their exported inputs. A point whose run is not in the table launches
it, and a later point with the same native input reuses it. The formulations also share the
solver identity of the problem. It is the `solver_identity` option when given, otherwise the
identity of the first run. Every later run must match it.
"""
function resolve(problem::LineParametersProblem,
        formulations::Union{Tuple{PSCADFormulation, Vararg{PSCADFormulation}},
            AbstractVector{<:PSCADFormulation}};
        options)
    execution = computation_options(PSCADFormulation, ComputationOptions(options))
    _pscad_deterministic(eltype(problem))
    validate(problem.frequencies, PSCADFormulation)
    _pscad_size(problem)
    blueprints = _pscad_blueprints(problem.system)
    physical_inputs = Engine.completed_inputs(problem)
    identity = Ref{Union{Nothing, Dict{String, String}}}(execution.data.solver_identity)
    runs = Dict{Any, Any}()
    return map(formulations) do formulation
        inputs = _pscad_inputs(problem, formulation, blueprints)
        key = (project = inputs.project,
            setting = inputs.setting[(:ground, :frequency, :configuration)])
        Resolved(problem, (; inputs, key, runs, identity, physical_inputs))
    end
end

"""
$(TYPEDSIGNATURES)

Compute the line parameters of one point with PSCAD: the native run of its exported inputs,
or the run that an earlier point of the problem launched for the same native input, checked
against the solver identity of the problem, and the result completed for its own
formulation.
"""
function compute(resolved::Resolved{<:LineParametersProblem}, formulation::PSCADFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    (; problem, input) = resolved
    normalized = computation_options(PSCADFormulation, ComputationOptions(options))
    execution = ComputationOptions(merge(normalized.data, (solver_identity = input.identity[],)))
    reused = haskey(input.runs, input.key)
    native = reused ? input.runs[input.key] :
             (input.runs[input.key] = _compute_pscad(problem, formulation, execution, input.inputs))
    input.identity[] === nothing && (input.identity[] = native.execution.solver_identity)
    value = _pscad_result(problem, formulation, input.inputs, native; batch_reuse = reused)
    return Engine.retain_gridpoint(value, gridpoint_id();
        fields = (formulation = formulation, choices = (), inputs = input.physical_inputs))
end

