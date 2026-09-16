function record_calculation(result::AbstractCoreResult, model)
    domain(result) === PhaseDomain || throw(ArgumentError("only phase-domain calculations can be retained"))
    return (kind=:gauntlet_calculation, frequencies=copy(frequencies(result)), basis=basis(result),
        domain=:PhaseDomain, Z=copy(observe(result, Z)), Y=copy(observe(result, Y)),
        comparison_unsupported=get(details(result).data, :comparison_unsupported, (;)),
        data_sha256=semantic_sha256(result, (; port_order=get(details(result).data,:coordinates,model.port_order))), computation_details=_selection_value(details(result).data))
end

function record_calculation(result::LineCableModels.AbstractUncertaintyResult, model)
    scientific_result=ImportExport.serialize_value(result)
    core=first(result)
    propagation = formulation_record(NamedTuple(result).formulation)
    sampling = result isa LineCableModels.MonteCarloResult ?
               (root_seed = LineCableModels.root_seed(result), point_seeds = [LineCableModels.point_seed(result,i) for i in eachindex(result)],
        trial_counts = [LineCableModels.trial_count(result,i) for i in eachindex(result)], distribution = LineCableModels.sampling_distribution(result)) :
               nothing
    return (kind = :gauntlet_uncertainty, scientific_result, frequencies = copy(frequencies(core)),
        basis = basis(core), domain = nameof(domain(core)), data_sha256 = semantic_sha256((scientific=scientific_result,port_order=get(details(core).data,:coordinates,model.port_order))),
        parameter_manifest = parameter_manifest(model), applied_variation = variation_record(model.variation),
        correlation = correlation_record(model), propagation, sampling,
        computation_details = _selection_value(LineCableModels.details(result).data))
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
        domain = getproperty.(points, :domain), data_sha256 = semantic_sha256(result, (; port_order=get(details(first(result)).data,:coordinates,model.port_order), axes)))
end

_numerical_record(record) = Base.structdiff(record, (; id=record.id))

# Native execution objects can contain deeply nested geometry and closures.
# Julia owns their serialization; JLD2 retains the portable numerical fields and
# evidence separately. Older typed JLD2 checkpoints remain readable.
function _execution_bytes(value)
    io=IOBuffer()
    Serialization.serialize(io,value)
    return take!(io)
end

function _read_execution(file,key)
    bytes_key=key*"_bytes"
    haskey(file,bytes_key) || return file[key]
    try
        return Serialization.deserialize(IOBuffer(file[bytes_key]))
    catch exception
        if key == "definitions" && exception isa TypeError
            throw(ArgumentError(
                "saved executable declarations do not match the current model layout; " *
                "rebuild them with `lcm gauntlet run --definition FILE --directory DIR` " *
                "without --resume. Previous attempts and results are preserved. " *
                "Deserialization error: " * sprint(showerror, exception)))
        end
        rethrow()
    end
end

_read_execution(path::AbstractString,key)=jldopen(file->_read_execution(file,key),path,"r")

function _campaign_request(definition)
    return (reference=_numerical_record(calculation_record(definition.reference)),
        candidate=_numerical_record(calculation_record(definition.candidate)),
        comparison=definition.comparison_settings,tolerances=definition.tolerances)
end

function _attempt_request(attempt)
    declaration=joinpath(attempt,"declarations.jld2")
    isfile(declaration) || return nothing
    request=jldopen(declaration,"r") do file
        haskey(file,"request") ? file["request"] : nothing
    end
    request === nothing || return request
    # Legacy completed attempts already retain portable calculation/report records.
    # Reading these fields does not restore a builder or load numerical arrays.
    operands=map((:reference,:candidate)) do role
        path=joinpath(attempt,string(role),"calculation.jld2")
        isfile(path) ? jldopen(file->file["calculation"],path,"r") : nothing
    end
    any(isnothing,operands) && return nothing
    analyses=joinpath(attempt,"analyses")
    isdir(analyses) || return nothing
    snapshots=[joinpath(folder,name) for (folder,_,names) in walkdir(analyses)
        for name in names if name=="snapshot.jld2"]
    isempty(snapshots) && return nothing
    snapshot=last(sort!(snapshots;by=mtime))
    comparison,tolerances=jldopen(snapshot,"r") do file
        (file["comparison_settings"],file["tolerances"])
    end
    return (reference=_numerical_record(operands[1]),
        candidate=_numerical_record(operands[2]),comparison,tolerances)
end

function _calculation_matches(directory,record)
    isfile(joinpath(directory,"complete.toml")) || return false
    path=joinpath(directory,"calculation.jld2")
    isfile(path) || throw(ArgumentError("completed calculation is missing: $path"))
    saved=jldopen(path,"r") do file
        haskey(file,"calculation") ? file["calculation"] : nothing
    end
    return saved !== nothing && _numerical_record(saved)==record
end

function _campaign_plan(root,definitions;resume=false,force=false)
    resume && force && throw(ArgumentError("force and resume are mutually exclusive"))
    return map(definitions) do definition
        benchmark_root=joinpath(root,string(definition.id))
        path=joinpath(benchmark_root,"state.toml")
        state=isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
        previous=get(state,"current",nothing)
        sources=unique([joinpath(benchmark_root,state[key]) for key in ("attempt","current")
            if haskey(state,key)])
        request=_campaign_request(definition)
        matched=force ? nothing : findfirst(path->_attempt_request(path)==request,sources)
        attempt=resume && haskey(state,"attempt") ? joinpath(benchmark_root,state["attempt"]) :
            matched === nothing ? nothing : sources[matched]
        saved_request=attempt === nothing ? nothing : _attempt_request(attempt)
        if resume && saved_request !== nothing
            all(getproperty(saved_request,role)==getproperty(request,role)
                for role in (:reference,:candidate)) ||
                throw(ArgumentError("campaign calculation inputs changed: $(definition.id)"))
        end
        if attempt !== nothing
            retained=jldopen(joinpath(attempt,"declarations.jld2"),"r") do file
                haskey(file,"reuse_attempts") ? file["reuse_attempts"] : nothing
            end
            # In particular, resuming a forced attempt cannot resurrect the
            # older calculations it was explicitly intended to replace.
            if retained !== nothing
                sources=unique([attempt;[joinpath(benchmark_root,path) for path in retained]])
            elseif resume
                sources=[attempt]
            end
        end
        completed=attempt !== nothing && previous !== nothing &&
            attempt==joinpath(benchmark_root,previous) &&
            saved_request==request
        actions=map((:reference,:candidate)) do role
            record=getproperty(request,role)
            !force && any(source->_calculation_matches(joinpath(source,string(role)),record),sources) ?
                :reuse : :compute
        end
        return (;id=definition.id,action=completed ? :skip : attempt === nothing ? :new : :continue,
            reference=actions[1],candidate=actions[2],definition,request,attempt,previous,
            sources=force ? String[] : sources)
    end
end

# A completed skip checks numerical/report payloads, not the multi-gigabyte
# native evidence tree. Full evidence validation belongs to the artifact readers.
function _skip_campaign_benchmark(root,id,attempt)
    state=TOML.parsefile(joinpath(root,string(id),"state.toml"))
    receiver=LineCableModels.progress_receiver()
    receiver === nothing || LineCableModels.report_progress(receiver,
        (kind=:benchmark,benchmark=id,attempt=relpath(attempt,joinpath(root,string(id))),
            state=:running,stage=:validating))
    declaration=joinpath(attempt,"declarations.jld2")
    bytes2hex(open(sha256,declaration))==strip(read(declaration*".sha256",String)) ||
        throw(ArgumentError("campaign declaration integrity check failed: $declaration"))
    digests=map((:reference,:candidate)) do role
        directory=joinpath(attempt,string(role))
        path=joinpath(directory,"calculation.jld2")
        marker=TOML.parsefile(joinpath(directory,"complete.toml"))
        digest=bytes2hex(open(sha256,path))
        digest==marker["sha256"] && first(split(read(path*".sha256",String)))==digest ||
            throw(ArgumentError("calculation checksum mismatch: $path"))
        digest
    end
    snapshots=[joinpath(folder,name) for (folder,_,names) in walkdir(joinpath(attempt,"analyses"))
        for name in names if name=="snapshot.jld2"]
    isempty(snapshots) && throw(ArgumentError("completed benchmark has no report: $attempt"))
    for path in snapshots
        bytes2hex(open(sha256,path))==first(split(read(path*".sha256",String))) ||
            throw(ArgumentError("benchmark checksum mismatch: $path"))
        operands=jldopen(file->file["calculations"],path,"r")
        all(operand.sha256==digest for (operand,digest) in zip(operands,digests)) ||
            throw(ArgumentError("benchmark operand checksum mismatch: $path"))
    end
    receiver=LineCableModels.progress_receiver()
    receiver === nothing || LineCableModels.report_progress(receiver,
        (kind=:benchmark,benchmark=id,state=:skipped,stage=:skipped,
            attempt=relpath(attempt,joinpath(root,string(id))),
            reused=true))
    return (id,state=:complete,result=nothing,identity=state["identity"],skipped=true)
end

function _retain_checkpoint(source,destination)
    # Completed files are immutable: hard links retain independent directory
    # ownership without duplicating native evidence. Writes use new files/rename.
    for (folder,_,names) in walkdir(source), name in names
        name in ("failure.toml","complete.toml") && continue
        original=joinpath(folder,name)
        target=joinpath(destination,relpath(original,source))
        mkpath(dirname(target))
        temporary=tempname(dirname(target))
        try
            try
                hardlink(original,temporary)
            catch
                cp(original,temporary;follow_symlinks=true)
            end
            mv(temporary,target;force=true)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
    # Completion records are published last, including independent sweep points.
    for (folder,_,names) in reverse(collect(walkdir(source)))
        "complete.toml" in names || continue
        target=joinpath(destination,relpath(folder,source),"complete.toml")
        mkpath(dirname(target))
        temporary=tempname(dirname(target))
        try
            cp(joinpath(folder,"complete.toml"),temporary)
            mv(temporary,target;force=true)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
end

function _saved_execution(directory,declaration;destination=nothing)
    path=joinpath(directory,"calculation.jld2")
    record=TOML.parsefile(joinpath(directory,"complete.toml"))
    saved=read_calculation(path;sha256_expected=record["sha256"])
    saved.metadata.calculation !== nothing &&
        _numerical_record(saved.metadata.calculation)==_numerical_record(declaration) ||
        throw(ArgumentError("calculation inputs changed: $directory"))
    result=_read_execution(path,"result")
    destination === nothing || _retain_checkpoint(directory,destination)
    original=saved.metadata.session === nothing ?
        (id="legacy",repository=saved.metadata.repository,active_project=saved.metadata.active_project) :
        saved.metadata.session
    return (;result,reused=true,session=original,timing=saved.metadata.timing)
end

function _execute(calculation::BenchmarkCalculation; directory = nothing, model = nothing,
        implementation = (), session = nothing, recover_solvers::Bool = false,
        reuse_directories=String[])
    started = time_ns()
    if directory === nothing
        began=time_ns()
        result=_compute_calculation(calculation)
        seconds=(time_ns()-began)*1e-9
        _accepted_scan_report(LineCableModels.progress_receiver(),result)
        return (;result,elapsed_seconds=seconds,reused=false,session,
            timing=(schema=1,scope=:compute_call_wall,seconds,
                callback_policy=:declared,diagnostic_policy=:declared,
                source_timings=_source_timings(result)))
    end
    validate(Base.write,directory)
    mkpath(directory)
    session === nothing && (session=execution_record())
    implementation === nothing && (implementation=())
    source_identity = [(; value.path, value.sha256) for value in implementation]
    declaration=calculation_record(calculation)
    signature=semantic_sha256(_numerical_record(declaration))
    path = joinpath(directory, "calculation.jld2")
    marker = joinpath(directory, "complete.toml")
    if isfile(marker)
        saved=_saved_execution(directory,declaration)
        _accepted_scan_report(LineCableModels.progress_receiver(),saved.result;reused=true)
        return (;saved...,elapsed_seconds=(time_ns()-started)*1e-9)
    end
    inputs=joinpath(directory,"inputs.toml")
    if isfile(inputs)
        TOML.parsefile(inputs)["signature"] == signature ||
            throw(ArgumentError("calculation inputs changed: $directory"))
    end
    for source in reuse_directories
        abspath(source)==abspath(directory) && continue
        _calculation_matches(source,_numerical_record(declaration)) || continue
        saved=_saved_execution(source,declaration;destination=directory)
        _accepted_scan_report(LineCableModels.progress_receiver(),saved.result;reused=true)
        return (;saved...,elapsed_seconds=(time_ns()-started)*1e-9)
    end
    if !isfile(inputs)
        _write_toml(inputs,Dict("signature"=>signature))
    end
    temporary = tempname(directory)
    try
        point_sessions=NamedTuple[]
        point_timings=NamedTuple[]
        reused_points=0
        jobs_reused=nothing
        compute_seconds=0.0
        receiver=LineCableModels.progress_receiver()
        receiver === nothing || LineCableModels.report_progress(receiver,(stage=:computing,
            backend=_execution_backend(calculation.formulation)))
        result = if calculation.problem isa LineParametersProblem && calculation.formulation isa Gridspace
            # A formulation sweep is a sequence of independently recoverable
            # scalar calculations. Preserve the public problem/formulation axes.
            formulations=collect(calculation.formulation)
            point_sources=[joinpath(parent,name)
                for source in reuse_directories
                for parent in (joinpath(source,"points"),) if isdir(parent)
                for name in readdir(parent) if isdir(joinpath(parent,name))]
            jobs_reused=0
            values=LineCableModels.with_scan_progress(;total=length(formulations)) do scan_receiver
              map(eachindex(formulations)) do index
                scan_receiver === nothing || LineCableModels.report_progress(scan_receiver,(kind=:scan,stage=:computing,
                    completed=index-1,total=length(formulations),formulation=index))
                options=calculation.options
                if haskey(options.data,:on_result) && options.data.on_result !== nothing
                    callback=options.data.on_result
                    options=Grammar.ComputationOptions(merge(options.data,(on_result=(problem,_,result)->callback(problem,index,result),)))
                end
                point=BenchmarkCalculation(calculation.id,calculation.problem,formulations[index];
                    options)
                value=LineCableModels.with_progress_scope(formulation=index) do
                    _execute(point;directory=joinpath(directory,"points",string(index)),
                        model,implementation,session,recover_solvers,
                        reuse_directories=point_sources)
                end
                push!(point_sessions,value.session)
                push!(point_timings,value.timing)
                reused_points += value.reused
                jobs_reused += value.reused || _result_reused(value.result; partial=false)
                compute_seconds += value.reused ? 0.0 : value.timing.seconds
                scan_receiver === nothing || LineCableModels.report_progress(scan_receiver,
                    (kind=:scan,stage=:computing,completed=index,total=length(formulations),reused=jobs_reused))
                value.result
              end
            end
            ParametricResult(LineCableModels.Combinatorial(calculation.formulation),values,
                (problems=[calculation.problem],formulations), Grammar.ComputationDetails())
        else
            options=calculation.options
            if recover_solvers && calculation.formulation isa Union{Engine.LineCableModelsFEM,PSCAD.PSCADFormulation}
                options=Grammar.ComputationOptions(merge((resume_run_directory=:latest,),options.data))
            end
            began=time_ns()
            value=_compute_calculation(calculation;options)
            compute_seconds=(time_ns()-began)*1e-9
            value
        end
        elapsed_seconds = (time_ns()-started)*1e-9
        timing=(schema=1,scope=:compute_call_wall,seconds=compute_seconds,
            callback_policy=:declared,diagnostic_policy=:declared,
            source_timings=_source_timings(result),points=point_timings,reused_points)
        _accepted_scan_report(receiver,result)
        receiver === nothing || LineCableModels.report_progress(receiver,(stage=:saving,))
        calculation_record(calculation) == declaration ||
            throw(ArgumentError("calculation inputs changed during execution"))
        retained_files=NamedTuple[]
        evidence=result isa ParametricResult ?
                 [(point = index, file) for (index, value) in enumerate(result)
                  for file in get(details(value).data, :files, ())] :
                 [(point = 0, file) for file in get(details(result).data, :files, ())]
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
        output_coordinates=get(details(result isa ParametricResult ? first(result) : result).data,:coordinates,model.port_order)
        JLD2.jldsave(temporary; schema_version = 3, status = :complete,
            result_bytes=_execution_bytes(result),payload...,
            retained_files, computation_signature = signature,
            problem = ImportExport.serialize_value(
                calculation.problem isa LineParametersProblem ? calculation.problem :
                model.nominal_problem),
            case_id = string(model.id), backend = string(nameof(typeof(calculation.formulation))),
            selection = string(calculation.id), formulation = declaration.formulation, calculation=declaration,
            repository=session.repository, active_project=session.active_project, session, point_sessions,
            port_order = copy(output_coordinates), implementation = source_identity, source_evidence = implementation,
            elapsed_at_completion_seconds = elapsed_seconds, timing)
        mv(temporary, path; force = true)
        digest = bytes2hex(open(sha256, path))
        write(path * ".sha256", digest * "  calculation.jld2\n")
        _write_toml(marker, Dict("schema"=>2, "signature"=>signature, "sha256"=>digest,
            "jobs"=>_calculation_jobs(calculation)))
        execution_wall=(time_ns()-started)*1e-9
        _write_toml(joinpath(directory,"timing.toml"),Dict(
            "schema"=>1,"execution_wall_seconds"=>execution_wall,
            "compute_call_seconds"=>compute_seconds,
            "scope"=>"execution through required result persistence; excludes this timing record"))
        return (; result, elapsed_seconds=execution_wall, reused = false, session,timing,jobs_reused)
    catch error
        _write_toml(joinpath(directory, "failure.toml"),
            Dict("message"=>sprint(showerror, error), "signature"=>signature, "time"=>string(now(UTC))))
        rethrow()
    finally
        isfile(temporary) && rm(temporary)
    end
end

"""
    run_campaign(directory, definitions; on_error=:continue, resume=false,
        force=false, benchmark=nothing, dry_run=false, recover_solvers=false,
        progress=:auto)

Reconcile selected definitions with retained work. Skip identical completed
benchmarks, continue identical drafts, and reuse matching operands and formulation
points when a changed declaration requires a new attempt. Save the entire selected
queue before entering a solver. Historical attempts remain readable.
Source edits alone do not invalidate completed calculations.
Each invocation records a new execution session; reused operands keep their original
execution records.

`benchmark` selects explicit IDs. `dry_run=true` returns per-operand scheduling
decisions without writing campaign files or checking full artifact integrity.
`force=true` requests fresh calculations for the selection. It is incompatible
with `resume=true`, which requires unchanged saved numerical declarations.
Skipped outcomes have `state=:complete`, `skipped=true` and `result=nothing`;
use `read_benchmark` when numerical results are needed.

Set `recover_solvers=true` to ask FEM and PSCAD to recover compatible native run
directories using their own input and output validation.

`progress=:auto` publishes snapshots and prints a command for a separate watcher
plus a final summary. `:plain` adds throttled single-line status at outer boundaries;
`:off` disables optional observation, estimation and publication. Controlled compute
samples suspend observation before the timing boundary. Operational wall time includes
monitoring; compute-call and native timing records retain their separate scopes.
"""
function run_campaign(directory::AbstractString, definitions::AbstractVector{<:BenchmarkDefinition};
        on_error::Symbol=:continue,resume::Bool=false,execution_sources=(),
        session=execution_record(),recover_solvers::Bool=false,progress::Symbol=:auto,
        force::Bool=false,benchmark=nothing,dry_run::Bool=false)
    ids=_campaign_selection(string.(getproperty.(definitions,:id)),benchmark)
    definitions=filter(definition->string(definition.id) in ids,definitions)
    isempty(definitions) && throw(ArgumentError("campaign needs benchmark definitions"))
    on_error in (:continue,:fail) || throw(ArgumentError("on_error must be :continue or :fail"))
    foreach(validate,definitions)
    plan=_campaign_plan(abspath(directory),definitions;resume,force)
    dry_run && return [(;item.id,item.action,item.reference,item.candidate) for item in plan]
    # A rejected immutable destination must not receive even a UI snapshot.
    validate(Base.write,directory)
    started=time_ns()
    return _with_campaign_progress(directory,getproperty.(definitions,:id),session.id;
            progress,selected=benchmark !== nothing) do tracker
        _record_campaign_wall(directory,session;progress,started) do
            _progress_declarations!(tracker,definitions,plan)
            _run_campaign(directory,definitions;on_error,resume,execution_sources,session,recover_solvers,plan)
        end
    end
end

# This measurement belongs to the execution owner, including when monitoring is
# off. An outer resume invocation replaces intermediate inner-run observations
# with its complete wall duration. It never includes downtime between invocations.
function _record_campaign_wall(f,directory,session;progress,started=time_ns())
    try
        return f()
    finally
        root=joinpath(abspath(directory),"sessions")
        if isfile(joinpath(root,session.id*".jld2"))
            _write_toml(joinpath(root,session.id*".timing.toml"),Dict(
                "schema"=>1,"session"=>session.id,"scope"=>"invocation_wall",
                "seconds"=>(time_ns()-started)*1e-9,"progress"=>string(progress),
                "includes"=>"preparation, calculations, reports, persistence and monitoring; excludes this final timing write"))
        end
    end
end

function _run_campaign(directory, definitions;
        on_error,resume,execution_sources,session,recover_solvers,plan)
    on_error in (:continue,:fail) || throw(ArgumentError("on_error must be :continue or :fail"))
    isempty(definitions) && throw(ArgumentError("campaign needs benchmark definitions"))
    allunique(getproperty.(definitions,:id)) || throw(ArgumentError("campaign benchmark IDs must be unique"))
    foreach(validate,definitions)
    root=abspath(directory)
    validate(Base.write,root)
    mkpath(root)
    manifest_path=joinpath(root,"campaign.toml")
    manifest_lock=open(joinpath(root,"execution.lock"),"a+")
    acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(manifest_lock),2,1)==0 :
        ccall(:flock,Cint,(Cint,Cint),fd(manifest_lock),6)==0
    acquired || (close(manifest_lock);throw(ArgumentError("another process is updating campaign $root")))
    try
        manifest=isfile(manifest_path) ? TOML.parsefile(manifest_path) :
            Dict("schema"=>3,"benchmarks"=>String[],"created"=>string(now(UTC)))
        manifest["schema"] == 3 || throw(ArgumentError("historical campaigns remain readable; run new drafts in a new staging directory"))
        manifest["benchmarks"]=sort!(unique(vcat(manifest["benchmarks"],string.(getproperty.(definitions,:id)))))
        tracker=_CAMPAIGN_TRACKER[]
        tracker === nothing || (tracker.selected |= length(tracker.rows)<length(manifest["benchmarks"]))
        _write_toml(manifest_path,manifest)
    finally
        close(manifest_lock)
    end
    sessions=joinpath(root,"sessions")
    mkpath(sessions)
    session_path=joinpath(sessions,session.id*".jld2")
    isfile(session_path) || JLD2.jldsave(session_path;session)
    prepared=NamedTuple[]
    # Prepare the entire queue before entering any solver. In particular, an
    # interruption in the first benchmark must not erase later declarations.
    for item in plan
        definition=item.definition
        benchmark_root=joinpath(root,string(definition.id))
        mkpath(benchmark_root)
        lease=open(joinpath(benchmark_root,"execution.lock"),"a+")
        acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
            ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
        acquired || (close(lease);throw(ArgumentError("another process owns benchmark $(definition.id)")))
        state_path=joinpath(benchmark_root,"state.toml")
        state=isfile(state_path) ? TOML.parsefile(state_path) : Dict{String,Any}()
        previous=get(state,"current",nothing)
        attempt=nothing
        try
            if item.attempt !== nothing
                attempt=item.attempt
                declaration=joinpath(attempt,"declarations.jld2")
                bytes2hex(open(sha256,declaration)) == strip(read(declaration*".sha256",String)) ||
                    throw(ArgumentError("campaign declaration integrity check failed"))
                if resume
                    saved=only(_read_execution(declaration,"definitions"))
                    for role in (:reference,:candidate)
                        _numerical_record(calculation_record(getproperty(saved,role))) ==
                            _numerical_record(calculation_record(getproperty(definition,role))) ||
                            throw(ArgumentError("campaign calculation inputs changed: $(definition.id)/$role"))
                    end
                end
                if state["attempt"] != relpath(attempt,benchmark_root)
                    state["attempt"]=relpath(attempt,benchmark_root)
                    state["state"]=item.action===:skip ? "complete" : "pending"
                    _write_toml(state_path,state)
                end
            else
                attempts=joinpath(benchmark_root,"attempts")
                mkpath(attempts)
                attempt=mktempdir(attempts;cleanup=false)
                declaration=joinpath(attempt,"declarations.jld2")
                declaration_paths=unique(vcat([definition.source_file,definition.model.source_file],
                    [normpath(joinpath(dirname(definition.model.source_file),asset)) for asset in definition.model.definition.assets]))
                declaration_sources=[(path,sha256=bytes2hex(open(sha256,path)),source=read(path)) for path in declaration_paths]
                execution_evidence=[(path=abspath(entry.path),module_name=entry.module_name,
                    source=read(entry.path),sha256=bytes2hex(open(sha256,entry.path))) for entry in execution_sources]
                dependencies=Set(entry.uuid for entry in session.packages if entry.extension_of === nothing)
                packages=[(name=id.name,uuid=string(id.uuid)) for id in keys(Base.loaded_modules)
                    if id.uuid !== nothing && string(id.uuid) in dependencies]
                JLD2.jldsave(declaration;definitions_bytes=_execution_bytes([definition]),implementation=(),session,
                    request=item.request,
                    reuse_attempts=[relpath(path,benchmark_root) for path in item.sources],
                    declaration_sources,execution_evidence,packages,on_error,case_sources=[definition.model.source_file],
                    repository=session.repository,active_project=session.active_project)
                write(declaration*".sha256",bytes2hex(open(sha256,declaration)))
                history=state
                state=Dict{String,Any}("schema"=>3,"attempt"=>relpath(attempt,benchmark_root))
                for key in ("fresh_wall_seconds","fresh_timing_key",
                        "fresh_reference_seconds","fresh_candidate_seconds")
                    haskey(history,key) && (state[key]=history[key])
                end
                previous === nothing || (state["current"]=previous)
                state["state"]="pending"
                _write_toml(state_path,state)
            end
            push!(prepared,(;definition,attempt,previous,reuse_directories=item.sources,
                skip=item.action===:skip))
        finally
            close(lease)
        end
    end
    outcomes=NamedTuple[]
    for (;definition,attempt,previous,reuse_directories,skip) in prepared
        benchmark_root=joinpath(root,string(definition.id))
        lease=open(joinpath(benchmark_root,"execution.lock"),"a+")
        acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
            ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
        acquired || (close(lease);throw(ArgumentError("another process owns benchmark $(definition.id)")))
        state_path=joinpath(benchmark_root,"state.toml")
        state=TOML.parsefile(state_path)
        if joinpath(benchmark_root,state["attempt"]) != attempt
            close(lease)
            throw(ArgumentError("benchmark attempt was replaced by another process: $(definition.id)"))
        end
        try
            if skip
                push!(outcomes,_skip_campaign_benchmark(root,definition.id,attempt))
                continue
            end
            completed=resume && get(state,"state","")=="complete"
            benchmark_started=time_ns()
            receiver=LineCableModels.progress_receiver()
            attempt_relative=relpath(attempt,benchmark_root)
            receiver === nothing || LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=definition.id,
                state=:running,stage=:preparing,attempt=attempt_relative))
            state["state"]="running"
            state["pid"]=getpid()
            state["session"]=session.id
            delete!(state,"message")
            _write_toml(state_path,state)
            value=LineCableModels.with_progress_scope(benchmark=definition.id,
                    attempt=attempt_relative) do
                run_benchmark(definition;directory=attempt,session,mode=:live,
                    measure_performance=!completed,recover_solvers,reuse_directories)
            end
            identity=semantic_sha256(read_benchmark,attempt)
            state["state"]="complete"
            state["current"]=relpath(attempt,benchmark_root)
            state["identity"]=identity
            state["wall_seconds"]=(time_ns()-benchmark_started)*1e-9
            state["reused"]=any(values(value.timings.execution)) do execution
                execution.reused || get(execution.compute, :reused_points, 0) > 0
            end
            if !state["reused"]
                state["fresh_wall_seconds"]=state["wall_seconds"]
                state["fresh_reference_seconds"]=value.timings.execution.reference.seconds
                state["fresh_candidate_seconds"]=value.timings.execution.candidate.seconds
            end
            _write_toml(state_path,state)
            push!(outcomes,(id=definition.id,state=:complete,result=value,identity))
            receiver === nothing || LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=definition.id,
                state=:complete,stage=:complete,reused=state["reused"],seconds=state["wall_seconds"],
                finalization_seconds=max(0.0, state["wall_seconds"] -
                    value.timings.execution.reference.seconds -
                    value.timings.execution.candidate.seconds)))
        catch error
            state["state"]=error isa InterruptException ? "interrupted" : "failed"
            state["message"]=sprint(showerror,error;context=:limit=>true)
            _write_toml(state_path,state)
            receiver=LineCableModels.progress_receiver()
            receiver === nothing || LineCableModels.report_progress(receiver,
                (kind=:benchmark,benchmark=definition.id,state=Symbol(state["state"]),
                    stage=Symbol(state["state"])))
            error isa InterruptException && rethrow()
            on_error === :fail && rethrow()
            push!(outcomes,(id=definition.id,state=:failed,message=state["message"]))
        finally
            close(lease)
        end
    end
    return outcomes
end

"""
    resume_campaign(directory; benchmark=nothing, recover_solvers=false, progress=:auto)

Continue the saved work order using checksummed declarations. Restore captured
declaration code for constructors without consulting the original source files.
Complete entries are skipped before restoring executable declarations or loading
results. Skip verification checks payload checksums; use artifact readers for full
backend-evidence validation. `benchmark` optionally selects explicit saved IDs.
Completed numerical operands keep their original execution records; unfinished
work uses the current Julia environment. A legacy queue without saved declarations
must first be supplied to `run_campaign` with `resume=true`.
One progress tracker covers the whole resumed invocation; `progress` has the same
meaning as in [`run_campaign`](@ref).
"""
function resume_campaign(directory::AbstractString;recover_solvers::Bool=false,progress::Symbol=:auto,
        benchmark=nothing)
    root=abspath(directory)
    validate(Base.write,root)
    manifest=TOML.parsefile(joinpath(root,"campaign.toml"))
    manifest["schema"] == 3 || throw(ArgumentError("historical campaigns remain readable; new execution requires a new staging directory"))
    manifest["benchmarks"]=_campaign_selection(manifest["benchmarks"],benchmark)
    outcomes=NamedTuple[]
    session=execution_record()
    started=time_ns()
    return _with_campaign_progress(root,manifest["benchmarks"],session.id;
            progress,selected=benchmark !== nothing) do tracker
        _record_campaign_wall(root,session;progress,started) do
            sessions=joinpath(root,"sessions")
            mkpath(sessions)
            JLD2.jldsave(joinpath(sessions,session.id*".jld2");session)
            _resume_campaign(root,manifest,outcomes,session;recover_solvers,progress)
        end
    end
end

function _resume_campaign(root,manifest,outcomes,session;recover_solvers,progress)
    for id in manifest["benchmarks"]
        state_path=joinpath(root,id,"state.toml")
        isfile(state_path) || throw(ArgumentError(
            "legacy campaign has no saved declaration for $id; supply the original definitions to run_campaign(...; resume=true)"))
        state=TOML.parsefile(state_path)
        if state["state"]=="complete"
            lease=open(joinpath(root,id,"execution.lock"),"a+")
            acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
                ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
            acquired || (close(lease);throw(ArgumentError("another process owns benchmark $id")))
            try
                current=TOML.parsefile(state_path)
                current==state || throw(ArgumentError("benchmark state changed during resume: $id"))
                push!(outcomes,_skip_campaign_benchmark(root,Symbol(id),joinpath(root,id,state["current"])))
            finally
                close(lease)
            end
            continue
        end
        receiver=LineCableModels.progress_receiver()
        receiver === nothing || LineCableModels.report_progress(receiver,
            (kind=:benchmark,benchmark=id,attempt=state["attempt"],state=:running,stage=:restoring))
        declaration=joinpath(root,id,state["attempt"],"declarations.jld2")
        bytes2hex(open(sha256,declaration)) == strip(read(declaration*".sha256",String)) ||
            throw(ArgumentError("campaign declaration integrity check failed"))
        evidence=jldopen(declaration,"r") do file
            (packages=file["packages"],execution=file["execution_evidence"],sources=file["declaration_sources"],case_sources=file["case_sources"])
        end
        for package in evidence.packages
            Base.require(Base.PkgId(Base.UUID(package.uuid),package.name))
        end
        snapshot=joinpath(dirname(declaration),"sources")
        captured_path(path)=joinpath(snapshot,splitpath(abspath(path))[2:end]...)
        for entry in vcat(evidence.execution,evidence.sources)
            bytes2hex(sha256(entry.source)) == entry.sha256 ||
                throw(ArgumentError("captured declaration integrity check failed: $(entry.path)"))
            target=captured_path(entry.path)
            mkpath(dirname(target))
            write(target,entry.source)
        end
        for entry in evidence.execution
            owner=entry.module_name === :Main ? Main : @__MODULE__
            Base.include(owner,captured_path(entry.path))
        end
        for path in evidence.case_sources
            Base.include(@__MODULE__,captured_path(path))
        end
        definitions=_read_execution(declaration,"definitions")
        on_error=JLD2.load(declaration,"on_error")
        append!(outcomes,Base.invokelatest(run_campaign,root,definitions;
            on_error,resume=true,session,recover_solvers,progress))
    end
    return outcomes
end

function _campaign_selection(ids,selection)
    allunique(ids) || throw(ArgumentError("campaign benchmark IDs must be unique"))
    selection === nothing && return ids
    selected=selection isa Union{Symbol,AbstractString} ? [string(selection)] : string.(collect(selection))
    isempty(selected) && throw(ArgumentError("benchmark selection cannot be empty"))
    allunique(selected) || throw(ArgumentError("benchmark selection contains duplicate IDs"))
    unknown=setdiff(selected,ids)
    isempty(unknown) || throw(ArgumentError("unknown benchmark IDs: $(join(unknown,", "))"))
    return filter(id->id in selected,ids)
end

"""
    campaign_status(directory; verify=true)

Report draft states and the identity of the last complete result. `verify=false`
reads only campaign metadata and ownership locks, without loading numerical data
or checking result integrity. Watching uses this lightweight read-only mode.
"""
function campaign_status(directory::AbstractString; verify::Bool=true)
    root=abspath(directory)
    manifest=TOML.parsefile(joinpath(root,"campaign.toml"))
    locked=isfile(joinpath(root,"bundle.toml"))
    return map(manifest["benchmarks"]) do id
        path=joinpath(root,id,"state.toml")
        isfile(path) || return (id,state=:missing,message="No completed attempt",identity=nothing,previous=false)
        record=TOML.parsefile(path)
        state=Symbol(record["state"])
        if state === :running && !locked
            lock_path=joinpath(root,id,"execution.lock")
            if !isfile(lock_path)
                return (id,state=:interrupted,message=get(record,"message",""),
                    identity=get(record,"identity",nothing),previous=haskey(record,"current"))
            end
            lease=open(lock_path,"r")
            acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
                ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
            close(lease)
            acquired && (state=:interrupted)
        end
        identity=state === :complete && verify ? semantic_sha256(read_benchmark,joinpath(root,id)) : get(record,"identity",nothing)
        return (id,state,message=get(record,"message",""),identity,
            previous=state !== :complete && haskey(record,"current"))
    end
end

export run_campaign,resume_campaign,campaign_status

function _calculation_jobs(calculation)
    problem, formulation = calculation.problem, calculation.formulation
    points = problem isa ParametricProblem ? length(problem.space) : 1
    forms = formulation isa Gridspace ? length(formulation) :
        formulation isa Union{LinearError,LineCableModels.Combinatorial} &&
        formulation.inner isa Union{Gridspace,AbstractVector} ? length(formulation.inner) : 1
    return points * forms
end

function _calculation_progress(calculation)
    jobs=_calculation_jobs(calculation)
    form=calculation.formulation
    total=form isa MonteCarlo ? (form.options.data.trials === nothing ? -1 : jobs*form.options.data.trials) : jobs
    batch=calculation.problem isa ParametricProblem && form isa Union{LinearError,LineCableModels.Combinatorial} &&
        form.inner isa Union{Gridspace,AbstractVector} ? length(form.inner) : 1
    return (backend=_execution_backend(form),mode=_execution_mode(form),total,batch,
        warmup=!_external_formulation(form))
end

_accepted_scans(result) = 1
_accepted_scans(result::ParametricResult) = sum(_accepted_scans,result;init=0)
_accepted_scans(result::LineCableModels.AbstractUncertaintyResult) = length(result)
_accepted_scans(result::LineCableModels.MonteCarloResult) = sum(i->LineCableModels.trial_count(result,i),eachindex(result);init=0)
_recovered_scans(result) = _result_reused(result;partial=false) ? _accepted_scans(result) : 0
_recovered_scans(result::ParametricResult) = sum(_recovered_scans,result;init=0)
function _accepted_scan_report(receiver,result;reused=false)
    receiver === nothing && return
    count=_accepted_scans(result)
    recovered=reused ? count : _recovered_scans(result)
    LineCableModels.report_progress(receiver,(kind=:scan_result,completed=count,total=count,
        reused=recovered,partial_recovery=recovered>0 || _result_reused(result)))
end
