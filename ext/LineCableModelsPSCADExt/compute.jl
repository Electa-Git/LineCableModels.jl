function computation_options(
        ::Type{PSCADFormulation},
        options::NamedTuple
)::ComputationOptions
    allowed = (:output_stem, :remote, :verbosity, :output_basis, :on_result,
        :resume_run_directory, :solver_identity, :work_root)
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
            solver_identity = nothing
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
    verbosity_values = normalized.verbosity
    verbosity_values isa NamedTuple || throw(ArgumentError(
        "verbosity must be a named tuple",
    ))
    haskey(verbosity_values, :default) || throw(ArgumentError(
        "verbosity must define a default level",
    ))
    all(value -> value isa Integer && value in 0:2, values(verbosity_values)) ||
        throw(ArgumentError("verbosity levels must be integers from 0 to 2"))
    basis_value = normalized.output_basis
    basis_value in (:pul, :total) || throw(ArgumentError(
        "output_basis must be :pul or :total; got $(repr(basis_value))",
    ))
    levels = NamedTuple{keys(verbosity_values)}(Int.(values(verbosity_values)))
    resume = normalized.resume_run_directory
    (resume === nothing || resume === :latest || resume isa AbstractString) ||
        throw(ArgumentError("PSCAD resume_run_directory must be nothing, :latest, or a completed run path"))
    resume isa AbstractString && isempty(resume) &&
        throw(ArgumentError(
            "PSCAD resume_run_directory cannot be empty"))
    identity = normalized.solver_identity
    (identity === nothing || identity isa AbstractDict{String, String}) ||
        throw(ArgumentError("PSCAD solver_identity must be the record returned by identify(remote)"))
    work_root = abspath(normalized.work_root)
    first(splitpath(relpath(work_root, options.remote.local_root))) == ".." &&
        throw(ArgumentError("PSCAD work_root must be inside remote.local_root"))
    return (
        work_root,
        output_stem,
        remote = options.remote,
        verbosity = levels,
        output_basis = Val(basis_value),
        on_result = normalized.on_result,
        resume_run_directory = resume isa AbstractString ? abspath(resume) : resume,
        solver_identity = identity
    )
end

function _formulation_label(formulation::PSCADFormulation)
    return join(
        (
            string(NamedTuple(formulation).requested.earth_impedance),
            "PSCAD native earth admittance",
            description(formulation.methods.insulation_admittance)
        ),
        '/')
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
        "PSCAD computation phase assignments must be positive",
    ))
    length(unique(assignments)) == length(assignments) || throw(ArgumentError(
        "PSCAD computation does not permit bundled terminals",
    ))
    sort(assignments) == collect(1:length(assignments)) || throw(ArgumentError(
        "PSCAD computation phase assignments must be contiguous from 1",
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

function _stage_pscad_project(
        problem::LineParametersProblem, formulation::PSCADFormulation,
        setting, work_root::AbstractString)
    parent = abspath(work_root)
    mkpath(parent)
    root = mktempdir(parent; prefix = "run-", cleanup = false)
    @info "Exporting PSCAD computation project" system = problem.system.system_id
    project = export_data(
        :pscad,
        problem.system,
        problem.earth_props;
        formulation = formulation,
        base_freq = formulation.options.base_frequency,
        temperature = problem.temperature,
        native_settings = setting[(:ground, :frequency)],
        file_name = joinpath(root, "generated.pscx")
    )
    project isa AbstractString && isfile(project) || throw(ArgumentError(
        "the PSCAD exporter did not create the requested project",
    ))
    staged = joinpath(root, "generated.pscx")
    project == staged || cp(project, staged; force = true)
    dielectric_losses = map(enumerate(problem.system.designs)) do (cable, design)
        map(enumerate(_pscad_components(design, formulation.options.base_frequency,
            formulation, problem.temperature))) do (layer, component)
            dielectric = component.dielectric
            requested = iszero(dielectric.shunt_capacitance) ? 0.0 :
                        dielectric.shunt_conductance /
                        (2pi * formulation.options.base_frequency *
                         dielectric.shunt_capacitance)
            exported = parse(Float64, _pscad_value(requested; maximum = 10))
            (; cable, layer, requested, exported, capped = requested > 10)
        end
    end
    return (; root, staged, dielectric_losses)
end

function _compute_pscad(problem::LineParametersProblem, formulation::PSCADFormulation,
        execution_options, prepared, setting)
    config = execution_options.remote
    root, staged = prepared.root, prepared.staged
    started = time_ns()
    input = Dict{String, Any}(
        "schema_version"=>3,
        "project_sha256"=>bytes2hex(open(sha256, staged)),
        "frequencies"=>Float64.(problem.frequencies),
        "matrix_size"=>collect(_pscad_size(problem)),
        "native_settings"=>Dict(string(component)=>Dict(string(field)=>Dict(
                                                            "value"=>control.value, "readback"=>control.readback)
                                for (field, control) in
                                    pairs(getproperty(setting, component)))
        for component in (:ground, :frequency, :configuration)),
        "solver"=>execution_options.solver_identity,
        "toolkit"=>Dict(name=>bytes2hex(sha256(source))
        for (name, source) in PSCAD_REMOTE_SOURCES))
    signature = bytes2hex(sha256(sprint(io -> TOML.print(io, input; sorted = true))))
    resume = execution_options.resume_run_directory
    candidates = resume === nothing ? String[] :
                 resume === :latest ?
                 sort!(
        filter(path -> isdir(path) && isfile(joinpath(path, "complete.toml")),
            readdir(dirname(root); join = true));
        by = path -> mtime(joinpath(path, "complete.toml")), rev = true) : [resume]
    source_root = nothing
    completion = nothing
    for candidate in candidates
        record_path = joinpath(candidate, "complete.toml")
        isfile(record_path) ||
            throw(ArgumentError("PSCAD run has no completion record: $candidate"))
        record = TOML.parsefile(record_path)
        if get(record, "input_sha256", nothing) != signature
            resume === :latest && continue
            throw(ArgumentError("PSCAD completed run has different numerical inputs or solver implementation: $candidate"))
        end
        get(record, "schema_version", nothing) == 3 || throw(ArgumentError(
            "unsupported PSCAD completion record: $record_path"))
        stored = TOML.parsefile(joinpath(candidate, "computation.toml"))
        bytes2hex(sha256(sprint(io -> TOML.print(io, stored; sorted = true)))) ==
        signature || throw(ArgumentError(
            "PSCAD completed-run input integrity check failed: $candidate"))
        bytes2hex(open(sha256, joinpath(candidate, "generated.pscx"))) ==
        stored["project_sha256"] ||
            throw(ArgumentError("PSCAD completed-run exported project changed: $candidate"))
        for (name, digest) in stored["toolkit"]
            bytes2hex(open(sha256, joinpath(candidate, "toolkit", name))) == digest ||
                throw(ArgumentError("PSCAD completed-run solver source changed: $candidate/toolkit/$name"))
        end
        for name in ("result_zm.out", "result_zp.out", "result_ym.out", "result_yp.out",
            "solver.toml", "native-settings.toml", "timing.txt")
            path = joinpath(candidate, "outputs", name)
            isfile(path) &&
            bytes2hex(open(sha256, path)) == get(record["outputs"], name, nothing) ||
                throw(ArgumentError("PSCAD completed-run output integrity check failed: $path"))
        end
        TOML.parsefile(joinpath(candidate, "outputs", "solver.toml")) == input["solver"] ||
            throw(ArgumentError("PSCAD completed run has no matching solver attestation: $candidate"))
        source_root, completion = candidate, record
        break
    end
    reused = source_root !== nothing
    if reused
        @info "PSCAD reuses a verified completed run" source_run=source_root
        output = joinpath(source_root, "outputs")
        execution = (elapsed_seconds = 0.0,
            elapsed_scope = "completed-run reuse; no solver execution", exit_code = 0,
            stdout_path = joinpath(output, "stdout.txt"), stderr_path = joinpath(output, "stderr.txt"),
            console_path = joinpath(output, "pscad-console.txt"), output_dir = output)
    else
        open(joinpath(root, "computation.toml"), "w") do io
            TOML.print(io, input; sorted = true)
        end
        @info "Computing PSCAD line parameters" system = problem.system.system_id
        execution = run_remote_pscad(config, staged, joinpath(root, "outputs"),
            formulation, problem.frequencies;
            output_stem = execution_options.output_stem,
            verbosity = verbosity(execution_options, :PSCAD))
        TOML.parsefile(joinpath(execution.output_dir, "solver.toml")) == input["solver"] ||
            throw(ArgumentError("PSCAD result has no matching solver attestation: $root"))
    end
    native_readback = TOML.parsefile(joinpath(execution.output_dir, "native-settings.toml"))
    expected_readback = Dict(component => Dict(field => control["readback"]
                             for (field, control) in controls)
    for (component, controls) in input["native_settings"])
    native_readback == expected_readback || throw(ArgumentError(
        "PSCAD result has no matching native-settings readback: $(execution.output_dir)"))
    parameters = try
        read_pscad_result(
            execution.output_dir,
            problem.frequencies,
            _pscad_size(problem)
        )
    catch error
        console_path=execution.console_path
        throw(ErrorException(
            "PSCAD result validation failed: $(sprint(showerror, error))" *
            "\nLast PSCAD diagnostics:\n$(_diagnostic_tail(console_path))" *
            "\nFull PSCAD diagnostics: $console_path",
        ))
    end
    source_elapsed = parse(Float64, strip(read(joinpath(execution.output_dir, "timing.txt"), String)))
    isfinite(source_elapsed) && source_elapsed >= 0 ||
        throw(ArgumentError("invalid PSCAD execution timing"))
    if !reused
        # Publish completion only after all four matrices have parsed and the
        # remote implementation has been checked. Interrupted runs stay intact.
        completion = Dict("schema_version"=>3, "input_sha256"=>signature,
            "outputs"=>Dict(name=>bytes2hex(open(sha256, joinpath(execution.output_dir, name)))
            for name in
                ("result_zm.out", "result_zp.out", "result_ym.out", "result_yp.out",
                "solver.toml", "native-settings.toml", "timing.txt")))
        temporary = tempname(root)
        try
            open(temporary, "w") do io
                TOML.print(io, completion; sorted = true)
            end
            mv(temporary, joinpath(root, "complete.toml"); force = false)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
    parameters = _pscad_basis(parameters, problem, execution_options.output_basis)
    evidence_root=reused ? source_root : root
    files=[(path = relpath(joinpath(directory, name), evidence_root),
               source = joinpath(directory, name), sha256 = bytes2hex(open(sha256, joinpath(directory, name))))
           for (directory, _, names) in walkdir(evidence_root) for name in sort(names)]
    names=["cable:$(terminal.cable):$(terminal.terminal)" for terminal in problem.system.terminal_order]
    coordinates=names[sortperm(problem.system.connection_order)]
    retained = (files, coordinates, requested_frequencies = copy(problem.frequencies),
        formulations = computation_details(formulation), native_setting = setting,
        base_frequency = formulation.options.base_frequency, loss_tangent_limit = 10.0, aerial_shunt_conductance = 1e-38,
        native_readback, native_frequencies = parameters.details.native_frequencies,
        dielectric_losses = prepared.dielectric_losses,
        exported_project = read(staged, String),
        execution = merge(execution,
            (backend = :pscad, pscad_version = config.pscad_version,
                reused, source_run = reused ? source_root : root,
                source_elapsed_seconds = source_elapsed,
                source_elapsed_scope = PSCAD_TIMING_SCOPE,
                wall_seconds = (time_ns() - started) * 1.0e-9,
                input_sha256 = signature, solver_identity = execution_options.solver_identity)))
    return LineParameters(
        parameters.domain, parameters.Z, parameters.Y, parameters.f, retained)
end

function compute(problem::LineParametersProblem, formulation::PSCADFormulation;
        options::NamedTuple = (;))
    return first(compute(problem, [formulation]; options))
end

function compute(
        problem::LineParametersProblem, formulations::AbstractVector{<:PSCADFormulation};
        options::NamedTuple = (;))
    isempty(formulations) &&
        throw(ArgumentError("PSCAD formulation collections cannot be empty"))
    settings = [pscad_setting(value, problem) for value in formulations]
    _validate_frequencies(problem.frequencies)
    _pscad_size(problem)
    execution = computation_options(PSCADFormulation, options)
    observed = identify(execution.remote)
    execution.solver_identity === nothing || execution.solver_identity == observed ||
        throw(ArgumentError("PSCAD solver installation changed during the campaign; start a new campaign"))
    execution = merge(execution, (solver_identity = observed,))
    projects = [_stage_pscad_project(problem, value, setting, execution.work_root)
                for (value, setting) in zip(formulations, settings)]
    # Same problem, frequency vector and execution settings throughout this
    # batch; reuse only byte-identical exported inputs and native solver choices.
    keys = [(project = read(project.staged, String),
                setting = setting[(:ground, :frequency)])
            for (project, setting) in zip(projects, settings)]
    first_result = _compute_pscad(
        problem, first(formulations), execution, first(projects), first(settings))
    values = Vector{typeof(first_result)}(undef, length(formulations))
    values[1] = first_result
    execution.on_result === nothing || execution.on_result(problem, 1, first_result)
    completed = Dict(first(keys)=>1)
    for index in 2:length(formulations)
        previous = get(completed, keys[index], nothing)
        if previous === nothing
            values[index] = _compute_pscad(
                problem, formulations[index], execution, projects[index], settings[index])
        else
            source = values[previous]
            @info "PSCAD reuses identical exported inputs" formulation=index source_formulation=previous
            retained = merge(deepcopy(source.details),
                (formulations = computation_details(formulations[index]),
                    native_setting = settings[index],
                    execution = merge(source.details.execution,
                        (reused = true, elapsed_seconds = 0.0,
                            elapsed_scope = "identical-input reuse; no solver execution", wall_seconds = 0.0))))
            values[index] = LineParameters(source.domain,
                SeriesImpedance(copy(source.Z.values); basis = basis(source)),
                ShuntAdmittance(copy(source.Y.values); basis = basis(source)), copy(source.f), retained)
        end
        completed[keys[index]] = index
        execution.on_result === nothing ||
            execution.on_result(problem, index, values[index])
    end
    return values
end
