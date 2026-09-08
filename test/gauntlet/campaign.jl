# Manual execution uses the same calculation records and compute dispatch as
# benchmark definitions; this file owns persistence, not numerical algorithms.
module CampaignCatalogue
include("fem_catalogue.jl")
end

function parse_selections(args)
    choices = Pair{Symbol, Any}[]
    for index in findall(==("--select"), args)
        index < length(args) || throw(ArgumentError("--select requires SLOT=ID,ID"))
        parts = split(args[index + 1], '='; limit=2)
        length(parts) == 2 || throw(ArgumentError("--select requires SLOT=ID,ID"))
        name = Symbol(strip(first(parts)))
        any(pair -> first(pair) === name, choices) && throw(ArgumentError(
            "duplicate --select for $name; list its identifiers in one comma-separated value"))
        identifiers = strip.(split(last(parts), ','))
        all(value -> occursin(r"^[A-Za-z][A-Za-z0-9_]*$", value), identifiers) ||
            throw(ArgumentError("--select $name needs nonempty formula identifiers"))
        push!(choices, name => Grid(Symbol.(identifiers)))
    end
    return (; choices...)
end

function campaign_selections(model, backend::Symbol, catalogue::Bool;
        choices::NamedTuple=(;), combine::Symbol=:product, dielectric::Symbol=:default)
    combine in (:product, :zip) || throw(ArgumentError("campaign combine must be product or zip"))
    if !isempty(choices)
        catalogue && throw(ArgumentError(
            "explicit --select axes require --formulas default; the catalogue is an independent earth-formula sweep"))
        slots = keys(Formulation().definitions)
        unknown = setdiff(keys(choices), slots)
        isempty(unknown) || throw(ArgumentError(
            "unknown formulation slots $(join(unknown, ", ")); choose from $(join(slots, ", "))"))
        defaults = merge(NamedTuple{slots}(map(_ -> :default, slots)),
            (insulation_admittance=dielectric, semicon_admittance=dielectric))
        # The public constructor owns product/zip, singleton broadcasting and
        # formula resolution. Campaign persistence only records its selections.
        space = Formulation(; merge(defaults, choices)..., combine)
        values = space isa Gridspace ? collect(space) : [space]
        selections = map(values) do value
            all(selection -> selection isa Symbol, value.definitions) || throw(ArgumentError(
                "manual campaign selections accept formula identifiers; route overrides belong to the Julia formulation API"))
            id = join((string(name, "_", lowercase(string(getproperty(value.definitions, name))))
                for name in keys(choices)), "__")
            return merge((; id), value.definitions)
        end
        return (; selections, skipped=NamedTuple[])
    end
    baseline = (id="default", earth_impedance=:default, earth_admittance=:default)
    selections = [baseline]
    skipped = NamedTuple[]
    catalogue || return (; selections, skipped)
    if backend === :pscad
        heights = getproperty.(model.problem.system.positions, :y)
        placement = all(>(0), heights) ? Val(:overhead) :
                    all(<(0), heights) ? Val(:underground) : Val(:mixed)
        identifiers = PSCADBenchmarks.formulas(placement)
        append!(selections, [(id=lowercase(string(id)), earth_impedance=id,
            earth_admittance=:default) for id in identifiers])
    else
        workspace = backend === :coaxial ? CampaignCatalogue.prepare_case(model) : nothing
        for record in CampaignCatalogue.catalogue()
            selected = CampaignCatalogue.variant(record)
            reason = backend === :coaxial ? CampaignCatalogue.case_skip_reason(model, selected, workspace) : nothing
            if reason !== nothing
                push!(skipped, (id=string(selected.id), reason))
                continue
            end
            push!(selections, (id=string(selected.id),
                earth_impedance=record.kind === :earth_impedance ? record.identifier : :default,
                earth_admittance=record.kind === :earth_admittance ? record.identifier : :default))
        end
    end
    return (; selections, skipped)
end

function campaign_formulation(backend::Symbol, selection, dielectric::Symbol)
    options = (reduce_bundle=false, kron_reduction=false,
        ideal_transposition=false, temperature_correction=true)
    requested = (; (Symbol(name) => Symbol(value) for (name, value) in pairs(selection)
        if name != "id")...)
    keywords = merge((insulation_admittance=dielectric, semicon_admittance=dielectric),
        requested, (; options))
    backend === :coaxial && return Formulation(; keywords...)
    backend === :pscad && return Formulation(:pscad; keywords...)
    if backend === :fem
        return Formulation(:LineCableModelsFEM; keywords...,
            fem_options=(gmsh_verbosity=0, getdp_verbosity=0,
                keep_run_directory=true))
    end
    throw(ArgumentError("unsupported campaign backend :$backend"))
end

function campaign_implementation(formulation::LineParametersFormulation)
    return implementation_record(formulation)
end

function campaign_implementation(formulation::LineCableModelsFEM)
    paths = [joinpath("ext", "LineCableModelsGmshExt", file) for file in (
        "LineCableModelsGmshExt.jl", "model.jl", "formulations.jl", "geometry.jl",
        "mesh.jl", "onelab.jl", "getdp.jl", "results.jl", "compute.jl",
        "getdp/model.pro", "getdp/quasi_tem.pro", "getdp/jacobian_integration.pro")]
    append!(paths, ["src/engine/formulations.jl", "src/engine/matrixops.jl",
        "src/engine/reduction.jl", "src/materials/material.jl"])
    for (family, selected) in (("insulationadmittance", formulation.methods.insulation_admittance),
            ("semiconadmittance", formulation.methods.semicon_admittance))
        push!(paths, "src/engine/$family/interface.jl")
        append!(paths, _formula_paths(family, selected))
    end
    physical = Formulation(; formulation.definitions..., options=formulation.options)
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    extension === nothing && error("load Gmsh before recording a FEM implementation")
    resolved = extension._getdp_selection(formulation)
    executable_identity = extension._getdp_identity(resolved.path)
    execution = merge(_selection_value(formulation.execution),
        (getdp_executable=nothing,))
    selection = merge(formulation_record(physical),
        (backend=:fem, execution,
            executable=executable_identity, gmsh_version=Gmsh.gmsh.GMSH_API_VERSION,
            gmsh_library=String(Gmsh.gmsh.lib), julia_version=string(VERSION)))
    return (selection, selection_sha256=semantic_sha256(selection),
        blobs=git_blob_record.(sort!(unique(paths))))
end

function campaign_implementation(formulation::PSCADBenchmarks.PSCADFormulation)
    paths = ["test/gauntlet/pscad/$file" for file in (
        "formulations.jl", "remote.jl", "outputs.jl", "remote/files.jl",
        "remote/runner.jl", "remote/supervisor.ps1", "remote/identity.py",
        "remote/Project.toml", "remote/Manifest.toml")]
    append!(paths, ["src/importexport/pscad/$file" for file in ("pscad.jl", "project.jl")])
    append!(paths, FLATTEN_IMPLEMENTATION_PATHS)
    append!(paths, ["src/engine/blueprint.jl",
        "src/engine/admittance.jl", "src/materials/material.jl"])
    for (family, selected) in (("insulationadmittance", formulation.methods.insulation_admittance),
            ("semiconadmittance", formulation.methods.semicon_admittance))
        push!(paths, "src/engine/$family/interface.jl")
        append!(paths, _formula_paths(family, selected))
    end
    selection = formulation_record(formulation)
    return (selection, selection_sha256=semantic_sha256(selection),
        blobs=git_blob_record.(sort!(unique(paths))))
end

_owned_formulation_record(value::PSCADBenchmarks.PSCADFormulation) = formulation_record(value)

function campaign_implementation(formulation::Union{LineCableModels.LinearError, LineCableModels.MonteCarlo})
    inner = campaign_implementation(formulation.inner)
    paths = ["src/grid.jl", "src/gridspace.jl", "src/uq/formulations.jl",
        "src/uq/results.jl", "src/grammar/uncertainty.jl",
        "src/engine/lineparameters/lineparameters.jl",
        "ext/LineCableModelsMeasurementsExt.jl", "test/gauntlet/case_loader.jl",
        "test/gauntlet/comparisons/uq_moments.jl"]
    append!(paths, formulation isa LineCableModels.LinearError ?
        ["src/uq/linearerror.jl", "src/parametricbuilder/traversal.jl"] :
        ["src/uq/montecarlo/compute.jl", "src/uq/statistics.jl"])
    selection = (physical=inner.selection, propagation=_owned_formulation_record(formulation),
        julia_version=string(VERSION), measurements_version=string(Base.pkgversion(Measurements)))
    return (selection, selection_sha256=semantic_sha256(selection),
        blobs=unique(vcat(inner.blobs, git_blob_record.(paths))))
end

function campaign_propagation(kind::Symbol, backend::Symbol, inner, settings)
    kind === :deterministic && return inner
    if kind === :linear_error
        backend === :coaxial || throw(ArgumentError(
            "$backend does not propagate Measurements derivatives; linear_error requires the coaxial backend"))
        return LineCableModels.LinearError(inner; options=(retain_details=true,))
    elseif kind === :monte_carlo
        backend === :fem && throw(ArgumentError("Gauntlet Monte Carlo on FEM is explicitly disabled"))
        return LineCableModels.MonteCarlo(inner; trials=settings["trials"],
            seed=parse(UInt64, settings["seed"]), distribution=:normal,
            return_samples=false, return_histograms=false, options=(retain_details=true,))
    end
    throw(ArgumentError("propagation must be deterministic, linear_error or monte_carlo"))
end

function campaign_models(ids, uncertainty)
    models = Dict{Tuple{Symbol, Bool}, LoadedCase}()
    for id in ids
        model = reference_case(id)
        models[(id, false)] = model
        uncertainty === nothing && continue
        selected = reference_grid(model.nominal_problem.frequencies)
        models[(id, true)] = load_case(id; variation=compose_variations(
            ExactOverrides(frequencies=selected), uncertainty))
        length(models[(id, true)].problem) == 1 || throw(ArgumentError(
            "a UQ campaign case must describe one uncertain design point"))
    end
    return models
end

function campaign_input(model, propagation::Symbol)
    nominal = numerical_input_sha256(model.nominal_problem)
    propagation === :deterministic && return nominal
    return semantic_sha256((nominal, parameters=parameter_manifest(model),
        variation=variation_record(model.variation), correlation=correlation_record(model)))
end

function record_calculation(result::LineParameters, model)
    data = (Z=vec(result.Z.values), Y=vec(result.Y.values), frequencies=result.f,
        port_order=model.port_order, basis=LineCableModels.basis(result))
    return (kind=:gauntlet_calculation, frequencies=copy(result.f), basis=data.basis,
        domain=:PhaseDomain, Z=copy(result.Z.values), Y=copy(result.Y.values),
        data_sha256=semantic_sha256(data), computation_details=LineCableModels.details(result))
end

function record_calculation(result::LineCableModels.AbstractUncertaintyResult, model)
    moments = NamedTuple(extract_moments(result, model.port_order))
    propagation = _owned_formulation_record(result.formulation)
    sampling = result isa LineCableModels.MonteCarloResult ?
        (root_seed=result.root_seed, point_seeds=copy(result.point_seeds),
            trial_counts=copy(result.trial_counts), distribution=result.formulation.distribution) : nothing
    return (kind=:gauntlet_moments, moments, frequencies=copy(moments.frequencies),
        basis=moments.basis, domain=moments.domain, data_sha256=semantic_sha256(moments),
        parameter_manifest=parameter_manifest(model), applied_variation=variation_record(model.variation),
        correlation=correlation_record(model), propagation, sampling,
        computation_details=LineCableModels.details(result))
end

function write_campaign_state(path, document)
    mkpath(dirname(path))
    temporary = tempname(dirname(path))
    try
        open(temporary, "w") do io
            TOML.print(io, document; sorted=true)
        end
        mv(temporary, path; force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

function run_campaign(directory::AbstractString, ids;
        backends=(:coaxial, :fem, :pscad), catalogue=true, dielectric=:default,
        choices::NamedTuple=(;), combine::Symbol=:product,
        propagation=(:deterministic,), uncertainty=nothing, trials=UQ_MONTE_CARLO_TRIALS, seed=nothing)
    haskey(ENV, "CI") && throw(ArgumentError("Gauntlet campaigns are manual, not CI simulations"))
    isempty(ids) && throw(ArgumentError("a campaign needs at least one case"))
    isempty(backends) && throw(ArgumentError("a campaign needs at least one backend"))
    all(in((:coaxial, :fem, :pscad)), backends) || throw(ArgumentError(
        "campaign backends must be coaxial, fem or pscad"))
    length(unique(ids)) == length(ids) || throw(ArgumentError("duplicate campaign cases"))
    length(unique(backends)) == length(backends) || throw(ArgumentError("duplicate campaign backends"))
    dielectric in (:default, :Ametani2004) || throw(ArgumentError(
        "campaign dielectric selection must be :default or :Ametani2004"))
    combine in (:product, :zip) || throw(ArgumentError("campaign combine must be product or zip"))
    explicit = isempty(choices) ? nothing : campaign_selections(nothing, :coaxial,
        catalogue; choices, combine, dielectric)
    isempty(propagation) && throw(ArgumentError("a campaign needs at least one propagation method"))
    allunique(propagation) || throw(ArgumentError("duplicate propagation methods"))
    all(in((:deterministic, :linear_error, :monte_carlo)), propagation) || throw(ArgumentError(
        "propagation must be deterministic, linear_error or monte_carlo"))
    uses_uncertainty = any(!=(:deterministic), propagation)
    if uses_uncertainty
        uncertainty isa RelativeStandardUncertainty || throw(ArgumentError(
            "UQ campaigns require explicit RelativeStandardUncertainty(percent; tags)"))
    else
        uncertainty === nothing || throw(ArgumentError("uncertainty requires a UQ propagation method"))
    end
    :monte_carlo in propagation && seed === nothing && throw(ArgumentError(
        "Monte Carlo campaigns require an explicit seed for reproducible resumption"))
    :monte_carlo in propagation || seed === nothing || throw(ArgumentError(
        "seed only applies to Monte Carlo propagation"))
    sampling = :monte_carlo in propagation ? Dict("trials"=>trials, "seed"=>string(UInt64(seed))) : false
    for backend in backends, method in propagation
        campaign_propagation(method, backend, Formulation(), sampling)
    end
    root = abspath(directory)
    ispath(root) && throw(ArgumentError("campaign directory already exists; use resume: $root"))
    models = campaign_models(ids, uncertainty)
    jobs = Dict{String, Any}[]
    for id in ids, backend in backends, method in propagation
        model = models[(id, method !== :deterministic)]
        selected = explicit === nothing ? campaign_selections(models[(id, false)], backend, catalogue) : explicit
        selections = [Dict(string(name)=>string(item) for (name, item) in pairs(value))
            for value in selected.selections]
        job_id = "$(id)_$backend" * (method === :deterministic ? "" : "_$method")
        push!(jobs, Dict("id"=>job_id, "case"=>string(id), "backend"=>string(backend),
            "propagation"=>string(method),
            "description"=>model.definition.description,
            "input_sha256"=>campaign_input(model, method),
            "selections"=>selections,
            "skipped"=>[Dict("id"=>value.id, "reason"=>value.reason) for value in selected.skipped]))
    end
    repository = repository_provenance()
    plan = Dict("schema_version"=>1, "created_at_utc"=>string(now(UTC)),
        "repository_commit"=>repository.commit, "repository_dirty"=>repository.dirty,
        "dielectric"=>string(dielectric), "combine"=>string(combine),
        "monte_carlo"=>sampling, "jobs"=>jobs)
    if uncertainty !== nothing
        plan["uncertainty"] = Dict("percent"=>uncertainty.percent, "tags"=>string.(collect(uncertainty.tags)))
    end
    mkpath(root)
    write_campaign_state(joinpath(root, "campaign.toml"), plan)
    return execute_campaign(root, plan, models)
end

function campaign_plan(directory)
    root = abspath(directory)
    path = joinpath(root, "campaign.toml")
    isfile(path) || throw(ArgumentError("campaign manifest is missing: $path"))
    plan = TOML.parsefile(path)
    get(plan, "schema_version", nothing) == 1 || throw(ArgumentError(
        "unsupported campaign schema in $path"))
    for job in plan["jobs"]
        occursin(_CASE_IDENTIFIER, job["id"]) || throw(ArgumentError("invalid campaign job ID"))
        Symbol(job["case"]) in keys(case_index()) || throw(ArgumentError("unknown campaign case $(job["case"])"))
        job["backend"] in ("coaxial", "fem", "pscad") || throw(ArgumentError("unknown campaign backend"))
        method = Symbol(get(job, "propagation", "deterministic"))
        campaign_propagation(method, Symbol(job["backend"]), Formulation(), plan["monte_carlo"])
        method === :deterministic || haskey(plan, "uncertainty") || throw(ArgumentError(
            "UQ campaign manifest is missing its uncertainty declaration"))
    end
    return root, plan
end

function resume_campaign(directory)
    haskey(ENV, "CI") && throw(ArgumentError("Gauntlet campaigns are manual, not CI simulations"))
    root, plan = campaign_plan(directory)
    ids = unique(Symbol(job["case"]) for job in plan["jobs"])
    uncertainty = haskey(plan, "uncertainty") ? RelativeStandardUncertainty(
        plan["uncertainty"]["percent"]; tags=Symbol.(plan["uncertainty"]["tags"])) : nothing
    models = campaign_models(ids, uncertainty)
    return execute_campaign(root, plan, models)
end

function campaign_status(directory)
    root, plan = campaign_plan(directory)
    lock_path = joinpath(root, "execution.lock")
    active = isfile(lock_path) && open(lock_path, "r") do io
        # Observe a live kernel lock, not merely the presence of its file.
        ccall(:flock, Cint, (Cint, Cint), Base.fd(io), 6) != 0
    end
    rows = NamedTuple[]
    for job in plan["jobs"]
        path = joinpath(root, job["id"], "state.toml")
        state = isfile(path) ? TOML.parsefile(path) : Dict("state"=>"pending", "completed"=>0)
        observed = state["state"] == "running" && !active ? "interrupted" : state["state"]
        push!(rows, (id=job["id"], state=observed, completed=state["completed"],
            requested=length(job["selections"]), skipped=length(job["skipped"]),
            message=get(state, "message", "")))
    end
    return rows
end

function execute_campaign(root, plan, models)
    haskey(ENV, "CI") && throw(ArgumentError("Gauntlet campaigns are manual, not CI simulations"))
    Sys.isunix() || throw(ArgumentError("the manual campaign CLI currently requires a Unix host for process locking"))
    lock = open(joinpath(root, "execution.lock"), "a+")
    if ccall(:flock, Cint, (Cint, Cint), Base.fd(lock), 6) != 0
        close(lock)
        throw(ArgumentError("another process is executing this campaign: $root"))
    end
    failures = 0
    try
        for job in plan["jobs"]
            propagation = Symbol(get(job, "propagation", "deterministic"))
            model = models[(Symbol(job["case"]), propagation !== :deterministic)]
            directory = joinpath(root, job["id"])
            state_path = joinpath(directory, "state.toml")
            completed = 0
            started = time_ns()
            attempt_path = joinpath(directory, "attempts", basename(tempname()) * ".toml")
            try
                campaign_input(model, propagation) == job["input_sha256"] || throw(ArgumentError(
                    "case numerical inputs changed; start a new campaign directory"))
                backend = Symbol(job["backend"])
                selections = job["selections"]
                formulations = [campaign_formulation(backend, selected, Symbol(plan["dielectric"]))
                    for selected in selections]
                # Capture implementation evidence before execution, never after a
                # remote solve that may outlive working-tree changes.
                calculations = [campaign_propagation(propagation, backend, inner, plan["monte_carlo"])
                    for inner in formulations]
                implementations = campaign_implementation.(calculations)
                options = backend === :fem ? (trace=true, resume_run_directory=:latest) : (;)
                if backend === :pscad
                    remote = PSCADBenchmarks._load_config()
                    solver_identity = PSCADBenchmarks.identify(remote)
                    options = (; remote, solver_identity, resume_run_directory=:latest)
                    implementations = [merge(value, (; solver_identity)) for value in implementations]
                end
                repository = repository_provenance()
                signatures = [semantic_sha256((input=job["input_sha256"], implementation=value))
                    for value in implementations]
                paths = [joinpath(directory, lpad(string(index), 4, '0') * ".jld2")
                    for index in eachindex(selections)]
                pending = Int[]
                for index in eachindex(paths)
                    if isfile(paths[index])
                        document = JLD2.load(paths[index])
                        document["computation_signature"] == signatures[index] || throw(ArgumentError(
                            "stored calculation $(paths[index]) has different numerical inputs or implementation; start a new campaign"))
                        document["status"] === :complete || throw(ArgumentError("incomplete stored calculation: $(paths[index])"))
                        data = get(document, "kind", nothing) === :gauntlet_moments ?
                            document["moments"] : (Z=vec(document["Z"]), Y=vec(document["Y"]),
                                frequencies=document["frequencies"], port_order=document["port_order"],
                                basis=document["basis"])
                        semantic_sha256(data) == document["data_sha256"] || throw(ArgumentError(
                            "stored numerical payload failed its integrity check: $(paths[index])"))
                        digest = bytes2hex(sha256(read(paths[index])))
                        checksum_path = paths[index] * ".sha256"
                        if isfile(checksum_path)
                            first(split(read(checksum_path, String))) == digest || throw(ArgumentError(
                                "stored artifact checksum differs: $(paths[index])"))
                        else
                            # Complete a sidecar interrupted after the atomic JLD2
                            # write, only after checking the internal payload hash.
                            write(checksum_path, digest * "  " * basename(paths[index]) * "\n")
                        end
                        completed += 1
                    else
                        push!(pending, index)
                    end
                end
                if isempty(pending)
                    write_campaign_state(state_path, Dict("state"=>"complete", "completed"=>completed))
                    println("REUSE\t", job["id"], "\t", completed)
                    continue
                end
                write_campaign_state(state_path, Dict("state"=>"running", "completed"=>completed,
                    "pid"=>getpid(), "started_at_utc"=>string(now(UTC))))
                println("BEGIN\t", job["id"], "\tpending=", length(pending))
                flush(stdout)
                chosen = formulations[pending]
                # One immutable input declaration per design point, shared by
                # its formulation callbacks. Replay must not need today's case file.
                problem_definition = LineCableModels.ImportExport.serialize_value(model.nominal_problem)
                # Gridspace owns the formulation axis; scalar problem normalization
                # and backend batch dispatch own lowering and numerical reuse.
                target = backend === :coaxial ? LineParametersFormulation :
                    backend === :fem ? LineCableModelsFEM : PSCADBenchmarks.PSCADFormulation
                space = Gridspace{target}(identity, (Grid(chosen),))
                started = time_ns()
                function save_result(problem, offset, result)
                    index = pending[offset]
                    temporary = tempname(directory)
                    try
                        payload = record_calculation(result, model)
                        JLD2.jldsave(temporary; schema_version=1, payload...,
                            status=:complete, numerical_reference_approval=:unreviewed,
                            case_id=string(model.id), input_sha256=job["input_sha256"],
                            problem=problem_definition,
                            formulation=(definitions=formulations[index].definitions,
                                options=formulations[index].options),
                            computation_signature=signatures[index], implementation=implementations[index],
                            backend, selection=selections[index],
                            repository_commit=repository.commit, repository_dirty=repository.dirty,
                            campaign_commit=plan["repository_commit"],
                            port_order=copy(model.port_order),
                            elapsed_at_completion_seconds=(time_ns() - started) * 1.0e-9,
                            batch_selection_count=length(pending),
                            recorded_at_utc=string(now(UTC)))
                        mv(temporary, paths[index]; force=false)
                        write(paths[index] * ".sha256", bytes2hex(sha256(read(paths[index]))) *
                            "  " * basename(paths[index]) * "\n")
                    finally
                        isfile(temporary) && rm(temporary)
                    end
                    completed += 1
                    write_campaign_state(state_path, Dict("state"=>"running", "completed"=>completed))
                    println("SAVED\t", job["id"], "\t", selections[index]["id"])
                    flush(stdout)
                    return nothing
                end
                if propagation === :deterministic
                    calculation = benchmark_calculation(Symbol(job["id"]),
                        backend === :coaxial ? :engine : :external,
                        LineCableModels.ParametricProblem(model.problem, merge(options, (on_result=save_result,))),
                        LineCableModels.Combinatorial(space; options=(retain_details=true,)))
                    _compute_owned(calculation)
                else
                    # UQ owns its uncertain Gridspace and sampling/derivatives.
                    # Checkpoint only whole moment results, never individual trials.
                    for (offset, inner) in enumerate(space)
                        calculation = benchmark_calculation(Symbol(job["id"]), :uq,
                            LineCableModels.ParametricProblem(model.problem, options),
                            campaign_propagation(propagation, backend, inner, plan["monte_carlo"]))
                        save_result(model.problem, offset, _compute_owned(calculation))
                    end
                end
                completed == length(selections) || throw(ArgumentError(
                    "backend did not report every completed formulation for checkpointing"))
                state = Dict("state"=>"complete", "completed"=>completed,
                    "batch_execution_seconds"=>(time_ns() - started) * 1.0e-9,
                    "recorded_at_utc"=>string(now(UTC)))
                write_campaign_state(attempt_path, state)
                write_campaign_state(state_path, state)
                println("COMPLETE\t", job["id"], "\t", completed)
            catch exception
                failures += 1
                message = sprint(showerror, exception)
                state = Dict("state"=>exception isa InterruptException ? "interrupted" : "failed",
                    "completed"=>completed, "message"=>message,
                    "elapsed_seconds"=>(time_ns() - started) * 1.0e-9,
                    "recorded_at_utc"=>string(now(UTC)))
                write_campaign_state(attempt_path, state)
                write_campaign_state(state_path, state)
                exception isa InterruptException && rethrow()
                println(stderr, "FAIL\t", job["id"], "\t", message)
            end
            flush(stdout)
            flush(stderr)
        end
    finally
        close(lock)
    end
    return failures == 0
end
