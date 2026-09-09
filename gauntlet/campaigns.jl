function record_calculation(result::AbstractCoreResult, model)
    domain(result) === PhaseDomain || throw(ArgumentError("only phase-domain calculations can be retained"))
    return (kind=:gauntlet_calculation, frequencies=copy(frequencies(result)), basis=basis(result),
        domain=:PhaseDomain, Z=copy(observe(result, Z)), Y=copy(observe(result, Y)),
        comparison_unsupported=get(details(result), :comparison_unsupported, (;)),
        data_sha256=semantic_sha256(result, (; port_order=model.port_order)), computation_details=_selection_value(details(result)))
end

function record_calculation(result::LineCableModels.AbstractUncertaintyResult, model)
    moments = NamedTuple(extract_moments(result, model.port_order))
    propagation = formulation_record(result.formulation)
    sampling = result isa LineCableModels.MonteCarloResult ?
               (root_seed = result.root_seed, point_seeds = copy(result.point_seeds),
        trial_counts = copy(result.trial_counts), distribution = result.formulation.distribution) :
               nothing
    return (kind = :gauntlet_moments, moments, frequencies = copy(moments.frequencies),
        basis = moments.basis, domain = moments.domain, data_sha256 = semantic_sha256(MomentResult(moments), (; port_order=model.port_order)),
        parameter_manifest = parameter_manifest(model), applied_variation = variation_record(model.variation),
        correlation = correlation_record(model), propagation, sampling,
        computation_details = _selection_value(LineCableModels.details(result)))
end

function record_calculation(result::ParametricResult, model)
    points=[record_calculation(value, model) for value in result]
    axes=(
        problems = [(terminal_order = copy(problem.system.terminal_order),
                        connection_order = copy(problem.system.connection_order),
                        declaration = ImportExport.serialize_value(problem))
                    for problem in NamedTuple(result).axes.problems],
        formulations = [_selection_value(formulation)
                        for formulation in NamedTuple(result).axes.formulations])
    return (kind = :gauntlet_result_space, points, axes,
        frequencies = getproperty.(points, :frequencies), basis = getproperty.(points, :basis),
        domain = getproperty.(points, :domain), data_sha256 = semantic_sha256(result, (; port_order=model.port_order, axes)))
end

function _execute(calculation::BenchmarkCalculation; directory = nothing, model = nothing,
        implementation = nothing)
    started = time_ns()
    keywords=isempty(calculation.options) ? (;) : (; options = calculation.options)
    directory === nothing && return (
        result = compute(calculation.problem,
            calculation.formulation; keywords...),
        elapsed_seconds = (time_ns()-started)*1e-9, reused = false)
    mkpath(directory)
    implementation === nothing && (implementation = implementation_record())
    source_identity = [(; value.path, value.sha256) for value in implementation]
    declaration=calculation_record(calculation)
    signature=semantic_sha256((
        numerical = Base.structdiff(declaration, (; id = declaration.id)),
        implementation = source_identity))
    path = joinpath(directory, "calculation.jld2")
    marker = joinpath(directory, "complete.toml")
    if isfile(marker)
        record = TOML.parsefile(marker)
        record["signature"] == signature || throw(ArgumentError(
            "calculation inputs or implementation changed: $directory"))
        bytes2hex(open(sha256, path)) == record["sha256"] || throw(ArgumentError(
            "calculation payload integrity check failed: $path"))
        read_calculation(path)
        result = JLD2.load(path, "result")
        return (; result, elapsed_seconds = 0.0, reused = true)
    end
    temporary = tempname(directory)
    try
        result = compute(calculation.problem, calculation.formulation; keywords...)
        elapsed_seconds = (time_ns()-started)*1e-9
        current = Dict(value.path=>value.sha256 for value in implementation_record())
        for value in implementation
            if !haskey(current, value.path)
                dependency_path=joinpath(REPOSITORY_ROOT, value.path)
                isfile(dependency_path) || throw(ArgumentError("execution source disappeared: $dependency_path"))
                current[value.path]=bytes2hex(open(sha256, dependency_path))
            end
        end
        current == Dict(value.path=>value.sha256 for value in source_identity) ||
            throw(ArgumentError("runtime sources changed during calculation"))
        calculation_record(calculation) == declaration ||
            throw(ArgumentError("calculation inputs changed during execution"))
        retained_files=NamedTuple[]
        evidence=result isa ParametricResult ?
                 [(point = index, file) for (index, value) in enumerate(result)
                  for file in get(details(value), :files, ())] :
                 [(point = 0, file) for file in get(details(result), :files, ())]
        for (point, file) in evidence
            isabspath(file.path) &&
                throw(ArgumentError("backend evidence paths must be relative"))
            first(splitpath(normpath(file.path))) == ".." &&
                throw(ArgumentError("backend evidence path escapes its bundle"))
            bytes2hex(open(sha256, file.source)) == file.sha256 ||
                throw(ArgumentError("backend evidence changed: $(file.source)"))
            target=joinpath(directory, "files", string(point), file.path)
            mkpath(dirname(target));
            cp(file.source, target; force = true, follow_symlinks = true)
            bytes2hex(open(sha256, target)) == file.sha256 ||
                throw(ArgumentError("copied evidence integrity check failed"))
            push!(retained_files,
                (point, path = relpath(target, directory),
                    sha256 = file.sha256, original_source = file.source))
        end
        payload = record_calculation(result, model)
        JLD2.jldsave(temporary; schema_version = 2, status = :complete, result, payload...,
            declaration = calculation, retained_files, computation_signature = signature,
            problem = ImportExport.serialize_value(
                calculation.problem isa LineParametersProblem ? calculation.problem :
                model.nominal_problem),
            case_id = string(model.id), backend = string(nameof(typeof(calculation.formulation))),
            selection = string(calculation.id), formulation = calculation_record(calculation),
            port_order = copy(model.port_order), implementation = source_identity, source_evidence = implementation,
            elapsed_at_completion_seconds = elapsed_seconds)
        mv(temporary, path; force = true)
        digest = bytes2hex(open(sha256, path))
        write(path * ".sha256", digest * "  calculation.jld2\n")
        _write_toml(marker, Dict("schema"=>2, "signature"=>signature, "sha256"=>digest))
        return (; result, elapsed_seconds, reused = false)
    catch error
        _write_toml(joinpath(directory, "failure.toml"),
            Dict("message"=>sprint(showerror, error), "signature"=>signature, "time"=>string(now(UTC))))
        rethrow()
    finally
        isfile(temporary) && rm(temporary)
    end
end

"Execute declared benchmarks without changing their problems, formulations or options."
function run_campaign(
        directory::AbstractString, definitions::AbstractVector{<:BenchmarkDefinition};
        on_error::Symbol = :continue, resume::Bool = false, execution_sources = ())
    on_error in (:continue, :fail) ||
        throw(ArgumentError("on_error must be :continue or :fail"))
    isempty(definitions) && throw(ArgumentError("campaign needs benchmark definitions"))
    allunique(getproperty.(definitions, :id)) ||
        throw(ArgumentError("campaign benchmark IDs must be unique"))
    foreach(validate, definitions)
    root = abspath(directory)
    mkpath(root)
    isfile(joinpath(root, "bundle.toml")) &&
        throw(ArgumentError("locked bundles are immutable; create a new campaign"))
    lock = open(joinpath(root, "execution.lock"), "a+")
    acquired = Sys.iswindows() ?
               ccall(:_locking, Cint, (Cint, Cint, Clong), fd(lock), 2, 1) == 0 :
               ccall(:flock, Cint, (Cint, Cint), fd(lock), 6) == 0
    acquired || (close(lock); throw(ArgumentError("another process owns campaign $root")))
    try
        source = implementation_record()
        for entry in execution_sources
            path=abspath(entry.path)
            relative=relpath(path, REPOSITORY_ROOT)
            any(value -> value.path==relative, source) && continue
            push!(source, (path = relative, sha256 = bytes2hex(open(sha256, path)),
                source = read(path)))
        end
        declaration_path = joinpath(root, "declarations.jld2")
        if isfile(declaration_path)
            resume || throw(ArgumentError("campaign already exists; use resume_campaign"))
            expected=strip(read(declaration_path*".sha256", String))
            bytes2hex(open(sha256, declaration_path)) == expected ||
                throw(ArgumentError("campaign declaration integrity check failed"))
            saved = JLD2.load(declaration_path)
            [(; x.path, x.sha256) for x in saved["implementation"]] ==
            [(; x.path, x.sha256) for x in source] ||
                throw(ArgumentError("campaign runtime sources changed"))
        else
            declaration_sources=map(sort!(unique(vcat(
                [value.source_file for value in definitions],
                [value.model.source_file for value in definitions],
                [normpath(joinpath(dirname(value.model.source_file), asset))
                 for value in definitions for asset in value.model.definition.assets])))) do path
                (path, sha256 = bytes2hex(open(sha256, path)), source = read(path))
            end
            temporary=tempname(root)
            execution_evidence=[(path = abspath(entry.path),
                                    module_name = entry.module_name,
                                    source = read(entry.path), sha256 = bytes2hex(open(sha256, entry.path)))
                                for entry in execution_sources]
            JLD2.jldsave(temporary; definitions, implementation = source,
                declaration_sources, execution_evidence, on_error)
            mv(temporary, declaration_path)
            write(declaration_path*".sha256", bytes2hex(open(sha256, declaration_path)))
            _write_toml(joinpath(root, "campaign.toml"),
                Dict("schema"=>2,
                    "benchmarks"=>string.(getproperty.(definitions, :id)),
                    "packages"=>[Dict("name"=>id.name, "uuid"=>string(id.uuid))
                                 for id in keys(Base.loaded_modules)
                                 if id.uuid !== nothing &&
                        haskey(Pkg.dependencies(), id.uuid)],
                    "execution_sources"=>[Dict("path"=>abspath(entry.path),
                                              "module"=>string(entry.module_name),
                                              "sha256"=>bytes2hex(open(sha256, entry.path)))
                                          for entry in execution_sources],
                    "case_sources"=>[value.model.source_file for value in definitions],
                    "case_source_sha256"=>[bytes2hex(open(sha256, value.model.source_file))
                                           for value in definitions], "created"=>string(now(UTC))))
        end
        outcomes = NamedTuple[]
        for definition in definitions
            output = joinpath(root, string(definition.id))
            state = joinpath(output, "state.toml")
            _write_toml(state, Dict("state"=>"running", "pid"=>getpid()))
            try
                value = run_benchmark(definition; directory = output, implementation = source, mode = :live)
                _write_toml(state, Dict("state"=>"complete"))
                push!(outcomes, (id = definition.id, state = :complete, result = value))
            catch error
                status = error isa InterruptException ? "interrupted" : "failed"
                _write_toml(state, Dict("state"=>status, "message"=>sprint(showerror, error)))
                error isa InterruptException && rethrow()
                on_error === :fail && rethrow()
                push!(outcomes, (id = definition.id, state = :failed,
                    message = sprint(showerror, error)))
            end
        end
        return outcomes
    finally
        close(lock)
    end
end

function resume_campaign(directory::AbstractString)
    path=joinpath(directory, "declarations.jld2")
    bytes2hex(open(sha256, path)) == strip(read(path*".sha256", String)) ||
        throw(ArgumentError("campaign declaration integrity check failed"))
    manifest=TOML.parsefile(joinpath(directory, "campaign.toml"))
    for package in get(manifest, "packages", [])
        Base.require(Base.PkgId(Base.UUID(package["uuid"]), package["name"]))
    end
    for source in get(manifest, "execution_sources", [])
        source_path=source["path"]
        isfile(source_path) && bytes2hex(open(sha256, source_path)) == source["sha256"] ||
            throw(ArgumentError("execution source changed or is unavailable: $source_path"))
        owner=source["module"] == "Main" ? Main : @__MODULE__
        Base.include(owner, source_path)
    end
    for (source, digest) in zip(manifest["case_sources"], manifest["case_source_sha256"])
        isfile(source) && bytes2hex(open(sha256, source)) == digest ||
            throw(ArgumentError("case builder source changed or is unavailable: $source"))
        Base.include(@__MODULE__, source)
    end
    declaration = JLD2.load(path)
    return Base.invokelatest(run_campaign, directory, declaration["definitions"];
        on_error = declaration["on_error"], resume = true,
        execution_sources = [(path = source["path"],
                                 module_name = Symbol(source["module"]))
                             for source in get(manifest, "execution_sources", [])])
end

function campaign_status(directory::AbstractString)
    root=abspath(directory)
    locked=isfile(joinpath(root, "bundle.toml"))
    lock=locked ? nothing : open(joinpath(root, "execution.lock"), "a+")
    available=locked || (Sys.iswindows() ?
               ccall(:_locking, Cint, (Cint, Cint, Clong), fd(lock), 2, 1) == 0 :
               ccall(:flock, Cint, (Cint, Cint), fd(lock), 6) == 0)
    lock === nothing || close(lock)
    return map(filter(path -> isfile(joinpath(path, "state.toml")), readdir(root; join = true))) do path
        record=TOML.parsefile(joinpath(path, "state.toml"))
        state=record["state"] == "running" && available ? :interrupted :
              Symbol(record["state"])
        (id = basename(path), state, message = get(record, "message", ""))
    end
end

export run_campaign, resume_campaign, campaign_status
