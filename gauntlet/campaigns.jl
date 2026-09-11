function record_calculation(result::AbstractCoreResult, model)
    domain(result) === PhaseDomain || throw(ArgumentError("only phase-domain calculations can be retained"))
    return (kind=:gauntlet_calculation, frequencies=copy(frequencies(result)), basis=basis(result),
        domain=:PhaseDomain, Z=copy(observe(result, Z)), Y=copy(observe(result, Y)),
        comparison_unsupported=get(details(result), :comparison_unsupported, (;)),
        data_sha256=semantic_sha256(result, (; port_order=get(details(result),:coordinates,model.port_order))), computation_details=_selection_value(details(result)))
end

function record_calculation(result::LineCableModels.AbstractUncertaintyResult, model)
    moments = NamedTuple(extract_moments(result, model.port_order))
    propagation = formulation_record(result.formulation)
    sampling = result isa LineCableModels.MonteCarloResult ?
               (root_seed = result.root_seed, point_seeds = copy(result.point_seeds),
        trial_counts = copy(result.trial_counts), distribution = result.formulation.distribution) :
               nothing
    return (kind = :gauntlet_moments, moments, frequencies = copy(moments.frequencies),
        basis = moments.basis, domain = moments.domain, data_sha256 = semantic_sha256(MomentResult(moments), (; port_order=get(details(result),:coordinates,model.port_order))),
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
        domain = getproperty.(points, :domain), data_sha256 = semantic_sha256(result, (; port_order=get(details(first(result)),:coordinates,model.port_order), axes)))
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

function _execute(calculation::BenchmarkCalculation; directory = nothing, model = nothing,
        implementation = (), session = nothing, recover_solvers::Bool = false)
    started = time_ns()
    keywords=isempty(calculation.options) ? (;) : (; options = calculation.options)
    directory === nothing && return (
        result = compute(calculation.problem,
            calculation.formulation; keywords...),
        elapsed_seconds = (time_ns()-started)*1e-9, reused = false, session)
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
        record = TOML.parsefile(marker)
        bytes2hex(open(sha256, path)) == record["sha256"] || throw(ArgumentError(
            "calculation payload integrity check failed: $path"))
        saved=read_calculation(path)
        # Older files included the source tree in their signature. Compare their
        # retained numerical declaration directly, without rewriting provenance.
        saved.metadata.calculation !== nothing &&
            _numerical_record(saved.metadata.calculation) == _numerical_record(declaration) ||
            throw(ArgumentError("calculation inputs changed: $directory"))
        result = _read_execution(path, "result")
        original=saved.metadata.session === nothing ?
            (id="legacy", repository=saved.metadata.repository,
                active_project=saved.metadata.active_project) : saved.metadata.session
        return (; result, elapsed_seconds = 0.0, reused = true, session=original)
    end
    inputs=joinpath(directory,"inputs.toml")
    if isfile(inputs)
        TOML.parsefile(inputs)["signature"] == signature ||
            throw(ArgumentError("calculation inputs changed: $directory"))
    else
        _write_toml(inputs,Dict("signature"=>signature))
    end
    temporary = tempname(directory)
    try
        point_sessions=NamedTuple[]
        result = if calculation.problem isa LineParametersProblem && calculation.formulation isa Gridspace
            # A formulation sweep is a sequence of independently recoverable
            # scalar calculations. Preserve the public problem/formulation axes.
            formulations=collect(calculation.formulation)
            values=map(eachindex(formulations)) do index
                options=calculation.options
                if haskey(options,:on_result) && options.on_result !== nothing
                    callback=options.on_result
                    options=merge(options,(on_result=(problem,_,result)->callback(problem,index,result),))
                end
                point=BenchmarkCalculation(calculation.id,calculation.problem,formulations[index];
                    options)
                value=_execute(point;directory=joinpath(directory,"points",string(index)),
                    model,implementation,session,recover_solvers)
                push!(point_sessions,value.session)
                value.result
            end
            ParametricResult(LineCableModels.Combinatorial(calculation.formulation),values,
                (problems=[calculation.problem],formulations), (;))
        else
            if recover_solvers && calculation.formulation isa Union{Engine.LineCableModelsFEM,PSCAD.PSCADFormulation}
                keywords=(options=merge((resume_run_directory=:latest,),calculation.options),)
            end
            compute(calculation.problem, calculation.formulation; keywords...)
        end
        elapsed_seconds = (time_ns()-started)*1e-9
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
        output_coordinates=get(details(result isa ParametricResult ? first(result) : result),:coordinates,model.port_order)
        JLD2.jldsave(temporary; schema_version = 2, status = :complete,
            result_bytes=_execution_bytes(result),payload...,
            retained_files, computation_signature = signature,
            problem = ImportExport.serialize_value(
                calculation.problem isa LineParametersProblem ? calculation.problem :
                model.nominal_problem),
            case_id = string(model.id), backend = string(nameof(typeof(calculation.formulation))),
            selection = string(calculation.id), formulation = declaration.formulation, calculation=declaration,
            repository=session.repository, active_project=session.active_project, session, point_sessions,
            port_order = copy(output_coordinates), implementation = source_identity, source_evidence = implementation,
            elapsed_at_completion_seconds = elapsed_seconds)
        mv(temporary, path; force = true)
        digest = bytes2hex(open(sha256, path))
        write(path * ".sha256", digest * "  calculation.jld2\n")
        _write_toml(marker, Dict("schema"=>2, "signature"=>signature, "sha256"=>digest))
        return (; result, elapsed_seconds, reused = false, session)
    catch error
        _write_toml(joinpath(directory, "failure.toml"),
            Dict("message"=>sprint(showerror, error), "signature"=>signature, "time"=>string(now(UTC))))
        rethrow()
    finally
        isfile(temporary) && rm(temporary)
    end
end

"""
    run_campaign(directory, definitions; on_error=:continue, resume=false)

Save every selected declaration before executing the first benchmark. A fresh run
replaces only those drafts once complete; resume continues recorded attempts with
unchanged numerical inputs. Source edits do not invalidate completed calculations.
Each invocation records a new execution session; reused operands keep their original
provenance. Previous results remain readable while a replacement is incomplete.

Set `recover_solvers=true` to ask FEM and PSCAD to recover compatible native run
directories using their own input and output validation.
"""
function run_campaign(directory::AbstractString, definitions::AbstractVector{<:BenchmarkDefinition};
        on_error::Symbol=:continue,resume::Bool=false,execution_sources=(),
        session=execution_record(),recover_solvers::Bool=false)
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
    for definition in definitions
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
            if resume && haskey(state,"attempt")
                attempt=joinpath(benchmark_root,state["attempt"])
                declaration=joinpath(attempt,"declarations.jld2")
                bytes2hex(open(sha256,declaration)) == strip(read(declaration*".sha256",String)) ||
                    throw(ArgumentError("campaign declaration integrity check failed"))
                saved=only(_read_execution(declaration,"definitions"))
                for role in (:reference,:candidate)
                    _numerical_record(calculation_record(getproperty(saved,role))) ==
                        _numerical_record(calculation_record(getproperty(definition,role))) ||
                        throw(ArgumentError("campaign calculation inputs changed: $(definition.id)/$role"))
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
                    declaration_sources,execution_evidence,packages,on_error,case_sources=[definition.model.source_file],
                    repository=session.repository,active_project=session.active_project)
                write(declaration*".sha256",bytes2hex(open(sha256,declaration)))
                state=Dict{String,Any}("schema"=>3,"attempt"=>relpath(attempt,benchmark_root))
                previous === nothing || (state["current"]=previous)
                state["state"]="pending"
                _write_toml(state_path,state)
            end
            push!(prepared,(;definition,attempt,previous))
        finally
            close(lease)
        end
    end
    outcomes=NamedTuple[]
    for (;definition,attempt,previous) in prepared
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
            completed=resume && get(state,"state","")=="complete"
            state["state"]="running"
            state["pid"]=getpid()
            state["session"]=session.id
            delete!(state,"message")
            _write_toml(state_path,state)
            value=run_benchmark(definition;directory=attempt,session,mode=:live,
                measure_performance=!completed,recover_solvers)
            read_benchmark(attempt)
            identity=semantic_sha256(read_benchmark,attempt)
            state["state"]="complete"
            state["current"]=relpath(attempt,benchmark_root)
            state["identity"]=identity
            _write_toml(state_path,state)
            if previous !== nothing && previous != state["current"]
                rm(joinpath(benchmark_root,previous);recursive=true)
            end
            push!(outcomes,(id=definition.id,state=:complete,result=value,identity))
        catch error
            state["state"]=error isa InterruptException ? "interrupted" : "failed"
            state["message"]=sprint(showerror,error;context=:limit=>true)
            _write_toml(state_path,state)
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
    resume_campaign(directory; recover_solvers=false)

Continue the saved work order using checksummed declarations. Restore captured
declaration code for constructors without consulting the original source files.
Completed numerical operands keep their original execution provenance; unfinished
work uses the current Julia environment. A legacy queue without saved declarations
must first be supplied to `run_campaign` with `resume=true`.
"""
function resume_campaign(directory::AbstractString;recover_solvers::Bool=false)
    root=abspath(directory)
    validate(Base.write,root)
    manifest=TOML.parsefile(joinpath(root,"campaign.toml"))
    manifest["schema"] == 3 || throw(ArgumentError("historical campaigns remain readable; new execution requires a new staging directory"))
    outcomes=NamedTuple[]
    session=execution_record()
    for id in manifest["benchmarks"]
        state_path=joinpath(root,id,"state.toml")
        isfile(state_path) || throw(ArgumentError(
            "legacy campaign has no saved declaration for $id; supply the original definitions to run_campaign(...; resume=true)"))
        state=TOML.parsefile(state_path)
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
            on_error,resume=true,session,recover_solvers))
    end
    return outcomes
end

"""Report current draft attempts and the identity of the last complete result."""
function campaign_status(directory::AbstractString)
    root=abspath(directory)
    manifest=TOML.parsefile(joinpath(root,"campaign.toml"))
    locked=isfile(joinpath(root,"bundle.toml"))
    return map(manifest["benchmarks"]) do id
        path=joinpath(root,id,"state.toml")
        isfile(path) || return (id,state=:missing,message="No completed attempt",identity=nothing,previous=false)
        record=TOML.parsefile(path)
        state=Symbol(record["state"])
        if state === :running && !locked
            lease=open(joinpath(root,id,"execution.lock"),"a+")
            acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
                ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
            close(lease)
            acquired && (state=:interrupted)
        end
        identity=state === :complete ? semantic_sha256(read_benchmark,joinpath(root,id)) : get(record,"identity",nothing)
        return (id,state,message=get(record,"message",""),identity,
            previous=state !== :complete && haskey(record,"current"))
    end
end

export run_campaign,resume_campaign,campaign_status
