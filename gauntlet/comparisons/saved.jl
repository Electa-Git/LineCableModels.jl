import DataFrames: metadata

"""Read a checksummed completed calculation without running a builder or solver."""
function read_calculation(path::AbstractString; sha256_expected = nothing, evidence::Symbol=:strict)
    evidence in (:strict,:numerical) || throw(ArgumentError("evidence must be :strict or :numerical"))
    path = abspath(path)
    isfile(path) && isfile(path * ".sha256") || throw(ArgumentError(
        "completed calculation and checksum required: $path"))
    digest = bytes2hex(open(sha256, path))
    checksum = split(read(path * ".sha256", String))
    !isempty(checksum) && first(checksum) == digest &&
    (sha256_expected === nothing || sha256_expected == digest) || throw(ArgumentError(
        "calculation checksum mismatch: $path"))
    # Read plain numerical and selection records only. In particular, do not
    # deserialize backend-specific computation_details from older implementations.
    document = jldopen(path, "r") do file
        names = ("schema_version", "kind", "status", "case_id", "backend", "problem",
            "formulation", "calculation", "repository", "active_project", "selection", "frequencies", "basis", "domain", "port_order",
            "Z", "Y", "moments", "scientific_result", "comparison_unsupported", "data_sha256",
            "implementation", "session", "point_sessions", "computation_signature", "elapsed_at_completion_seconds",
            "batch_selection_count", "timing", "propagation", "sampling", "parameter_manifest",
            "applied_variation", "correlation", "computation_details", "retained_files", "points", "axes")
        Dict(name => file[name] for name in names if haskey(file, name))
    end
    document["schema_version"] in (1, 2, 3) && document["status"] === :complete ||
        throw(ArgumentError("unsupported or incomplete calculation: $path"))
    evidence_issues=String[]
    for file in get(document, "retained_files", ())
        (isabspath(file.path) || first(splitpath(normpath(file.path))) == "..") &&
            throw(ArgumentError("retained evidence must be inside its calculation directory"))
        evidence_path=joinpath(dirname(path), file.path)
        isfile(evidence_path) && bytes2hex(open(sha256, evidence_path)) == file.sha256 ||
            push!(evidence_issues,"retained backend evidence is missing or changed: $evidence_path")
    end
    evidence === :strict && !isempty(evidence_issues) && throw(ArgumentError(join(evidence_issues,"\n")))
    kind = document["kind"]
    result = if kind === :gauntlet_calculation
        document["domain"] === :PhaseDomain ||
            throw(ArgumentError("only phase-domain calculations are supported"))
        retained=get(document,"computation_details",(;))
        declaration=get(retained,:formulations,document["formulation"])
        selected=ImportExport.deserialize_value(Val(:formulation),declaration)
        ismissing(selected) || (retained=merge(retained,Engine.completed_formulation(selected,declaration)))
        haskey(retained,:shunt_model) && (retained=merge(retained,
            NamedTuple{(:shunt_model,),Tuple{NamedTuple}}((retained.shunt_model,))))
        LineParameters(LineCableModels.PhaseDomain, document["Z"], document["Y"],
            document["frequencies"]; basis = document["basis"],
            details = Grammar.ComputationDetails(merge(retained,
                (gridpoint=get(retained,:gridpoint,(source_id=digest,problem_index=1,formulation_index=1)),
                    formulations=get(retained,:formulations,document["formulation"]),
                    coordinates=get(retained,:coordinates,String.(document["port_order"])),
                    comparison_unsupported = get(document, "comparison_unsupported", (;)),
                    files = [(path = file.path,
                                 source = joinpath(dirname(path), file.path),
                                 sha256 = file.sha256)
                             for file in get(document, "retained_files", ())]))))
    elseif kind === :gauntlet_result_space
        values=map(enumerate(document["points"])) do (index,point)
            retained=point.computation_details
            declaration=retained.formulations
            selected=ImportExport.deserialize_value(Val(:formulation),declaration)
            ismissing(selected) || (retained=merge(retained,Engine.completed_formulation(selected,declaration)))
            haskey(retained,:shunt_model) && (retained=merge(retained,
                NamedTuple{(:shunt_model,),Tuple{NamedTuple}}((retained.shunt_model,))))
            LineParameters(
                    PhaseDomain, point.Z, point.Y, point.frequencies; basis = point.basis,
                    details = Grammar.ComputationDetails(merge(retained,
                        (files = [(path = file.path,
                                      source = joinpath(dirname(path), file.path), sha256 = file.sha256)
                                  for file in get(document, "retained_files", ())
                                  if file.point == index],))))
        end
        for (point, value) in zip(document["points"], values)
            semantic_sha256(value, (; port_order=document["port_order"])) == point.data_sha256 ||
                throw(ArgumentError("result-space point data changed"))
        end
        ParametricResult(nothing, values, document["axes"], Grammar.ComputationDetails())

    elseif kind === :gauntlet_uncertainty
        ImportExport.deserialize_value(document["scientific_result"])
    else
        throw(ArgumentError("not a saved Gauntlet calculation: $path"))
    end
    ports = String.(document["port_order"])
    allunique(ports) &&
    length(ports) == (kind === :gauntlet_calculation ?
     size(result.Z, 1) :
     kind === :gauntlet_result_space ? size(first(result).Z, 1) :
     size(observe(first(result),Z),1)) ||
        throw(DimensionMismatch("invalid saved terminal order: $path"))
    coordinates=(port_order=ports, axes=get(document, "axes", nothing))
    data_digest = kind === :gauntlet_uncertainty ? semantic_sha256((scientific=document["scientific_result"],port_order=ports)) :
        semantic_sha256(result, coordinates)
    get(document, "data_sha256", data_digest) == data_digest ||
        throw(ArgumentError("calculation data digest mismatch: $path"))
    selection=document["selection"]
    if selection isa AbstractDict
        selection=(; (Symbol(k)=>v for (k,v) in selection)...)
    elseif selection isa Union{AbstractString, Symbol}
        selection=(id=Symbol(selection),)
    end
    selection isa NamedTuple || throw(ArgumentError("saved selection must be a named record"))
    formulation=document["formulation"]
    if formulation isa NamedTuple && all(key -> haskey(formulation,key),(:id,:input_sha256,:formulation,:options))
        formulation=formulation.formulation
    end
    # Executable checkpoints belong to execution reuse, never passive inspection.
    recovery=:portable
    metadata = (path, sha256 = digest, data_sha256 = data_digest,
        loaded_data_sha256=result isa AbstractUncertaintyResult ? data_digest : semantic_sha256(result,coordinates), recovery, evidence_issues,
        case_id = string(document["case_id"]),
        backend = string(document["backend"]), selection = selection,
        axes=coordinates.axes,
        formulation, calculation=get(document,"calculation",nothing),
        repository=get(document,"repository",nothing),active_project=get(document,"active_project",nothing), implementation = get(document, "implementation", nothing),
        session=get(document,"session",nothing), point_sessions=get(document,"point_sessions",()),
        computation_signature = get(document, "computation_signature", nothing),
        input_sha256 = semantic_sha256(document["problem"]), port_order = ports,
        frequencies = document["frequencies"], basis = document["basis"], domain = document["domain"],
        propagation = get(document, "propagation", nothing), sampling = get(document, "sampling", nothing),
        uncertainty = (parameters = get(document, "parameter_manifest", nothing),
            variation = get(document, "applied_variation", nothing), correlation = get(
                document, "correlation", nothing)),
        timing = get(document,"timing",(
            scope = document["schema_version"] == 1 ? :batch_elapsed_at_completion :
                    :legacy_execution_elapsed,
            seconds = get(document, "elapsed_at_completion_seconds", missing),
            batch_selections = get(document, "batch_selection_count", missing))))
    return (; result, metadata)
end

"""
Persist comparisons of saved operands; no numerical execution or CI-reference approval is performed.
"""
function compare_saved(benchmark::BenchmarkDefinition; directory::AbstractString)
    validate(Base.write,directory)
    validate(benchmark)
    publication=report(BenchmarkTableDefinition(;benchmark.comparison_settings...),
        (reference=benchmark.reference.problem, candidate=benchmark.candidate.problem,
            context=(id=benchmark.id,case_id=benchmark.case_id,collection=benchmark.collection));
        requests=benchmark.candidate.problem.result isa AbstractUncertaintyResult ?
            (R,L,G,C,filter(request -> request isa Tuple && first(request)===LineCableModels.statistics,benchmark.comparison_settings.requests)...) : ())
    return record_benchmark(benchmark, publication; directory)
end

"""
    record_benchmark(benchmark, comparisons; directory)

Write completed comparisons with their exact settings and checksummed operands.
No comparison or solver is executed. Each invocation writes a distinct output record.
"""
function record_benchmark(benchmark::BenchmarkDefinition,artifact::ReportArtifact;directory::AbstractString)
    validate(Base.write,directory)
    a,b=benchmark.reference.problem.metadata,benchmark.candidate.problem.metadata
    for operand in (a,b)
        bytes2hex(open(sha256,operand.path))==operand.sha256 ||
            throw(ArgumentError("saved operand changed during comparison: $(operand.path)"))
    end
    record_id=string(Grammar.gridpoint_id().source_id)
    path=joinpath(abspath(directory),string(benchmark.id),record_id,"snapshot.jld2")
    points=(artifact.observed isa ObservedResult ? (artifact.observed,) : artifact.observed)
    comparisons=[(request=ImportExport.serialize_value(row.request,Val(:scientific)),
        quantity=Symbol(Units.symbol(row.quantity)),statistic=row.statistic,
        candidate_id=row.candidate_id,reference_id=row.reference_id,
        absolute=row.absolute,relative=row.relative,details=row.settings) for point in points for row in point.errors]
    mkpath(dirname(path))
    temporary=tempname(dirname(path))
    try
        JLD2.jldsave(temporary;schema_version=3,kind=:gauntlet_benchmark,record_id,
            observed_data=ImportExport.serialize_value(artifact),
            summary=metadata(artifact.tables.maxima,"comparison_records"),
            formulations=[NamedTuple(row) for row in eachrow(artifact.tables.formulations)],
            benchmark_id=string(benchmark.id),case_id=string(benchmark.case_id),
            description=benchmark.model isa LoadedCase ? benchmark.model.definition.description : benchmark.model.description,
            collection=string(benchmark.collection),benchmark_source_sha256=benchmark.source_sha256,
            calculations=(reference=merge(a,(path=relpath(a.path,dirname(path)),id=benchmark.reference.id)),
                candidate=merge(b,(path=relpath(b.path,dirname(path)),id=benchmark.candidate.id))),
            comparison_settings=ImportExport.serialize_value(benchmark.comparison_settings,Val(:scientific)),
            tolerances=ImportExport.serialize_value(benchmark.tolerances,Val(:scientific)),
            port_order=a.port_order,frequencies=a.frequencies,basis=a.basis,domain=a.domain,
            reference_comparison=comparisons,timings=(reference=a.timing,candidate=b.timing),
            recorded_at_utc=string(now(UTC)))
        mv(temporary,path;force=false)
        write(path*".sha256",bytes2hex(open(sha256,path))*"  snapshot.jld2\n")
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

"""
    read_benchmark(path; load_results=false)

Read an explicitly ordered benchmark and verify the snapshot and operand hashes.
A benchmark directory returns `(id, reference, candidate, analyses, measurements)`, with loaded
results and metadata. For a snapshot file, the default returns the recorded
dictionary; `load_results=true` returns the same loaded bundle as a directory,
resolving operand paths relative to the snapshot. Loading performs no solve.

The loaded bundle supports `report(BenchmarkTableDefinition(false), benchmark)`
and `plot(benchmark, ydata; ...)` for REPL tables and matrix-cell overlays.
"""
function read_benchmark(path::AbstractString; load_results::Bool = false, previous::Bool=false,evidence::Symbol=:strict)
    if isdir(path)
        state_path=joinpath(path,"state.toml")
        if isfile(state_path)
            state=TOML.parsefile(state_path)
            if haskey(state,"current") || haskey(state,"attempt")
                (state["state"] == "complete" || previous) || throw(ArgumentError(
                    "latest draft is $(state["state"]); use previous=true to inspect the previous complete result"))
                haskey(state,"current") || throw(ArgumentError("no previous completed result"))
                current=state["current"]
                isabspath(current) || first(splitpath(normpath(current))) == ".." ?
                    throw(ArgumentError("invalid draft path")) : nothing
                return read_benchmark(joinpath(path,current);load_results,evidence)
            end
        end
        reference=read_calculation(joinpath(path, "reference", "calculation.jld2");evidence)
        candidate=read_calculation(joinpath(path, "candidate", "calculation.jld2");evidence)
        analyses=[read_benchmark(joinpath(folder, name))
                  for (folder, _, names) in
                      walkdir(joinpath(path, "analyses"))
                  for name in sort(names) if name == "snapshot.jld2"]
        isempty(analyses) &&
            throw(ArgumentError("benchmark has no retained analysis: $path"))
        return read_benchmark((id = Symbol(first(analyses)["benchmark_id"]), reference, candidate, analyses))
    end
    isfile(path) && isfile(path * ".sha256") ||
        throw(ArgumentError("benchmark record and checksum required: $path"))
    checksum = split(read(path * ".sha256", String))
    !isempty(checksum) && first(checksum) == bytes2hex(open(sha256, path)) ||
        throw(ArgumentError("benchmark checksum mismatch: $path"))
    record = JLD2.load(path)
    get(record,"schema_version",nothing)==3 && get(record,"kind",nothing)===:gauntlet_benchmark &&
        haskey(record,"observed_data") || throw(ArgumentError("a current observed benchmark record is required: $path"))
    record["comparison_settings"]=ImportExport.deserialize_value(record["comparison_settings"])
    record["tolerances"]=ImportExport.deserialize_value(record["tolerances"])
    record["reference_comparison"]=[merge(row,(request=ImportExport.deserialize_value(row.request),))
        for row in record["reference_comparison"]]
    operands = record["calculations"]
    keys(operands) == (:reference, :candidate) ||
        throw(ArgumentError("benchmark requires one explicit reference and candidate"))
    for operand in operands
        source = isabspath(operand.path) ? operand.path :
                 normpath(joinpath(dirname(path), operand.path))
        isfile(source) && bytes2hex(open(sha256, source)) == operand.sha256 ||
            throw(ArgumentError("benchmark operand missing or checksum mismatch: $source"))
    end
    if load_results
        reference, candidate = map(operands) do operand
            source = isabspath(operand.path) ? operand.path :
                     normpath(joinpath(dirname(path), operand.path))
            read_calculation(source; sha256_expected = operand.sha256,evidence)
        end
        return read_benchmark((id = Symbol(record["benchmark_id"]), reference, candidate, analyses = [record]))
    end
    return record
end

"""Bind recorded performance evidence to the already loaded numerical operands."""
function read_benchmark(benchmark::NamedTuple{(:id,:reference,:candidate,:analyses)})
    measurements=(execution=(reference=benchmark.reference.metadata,candidate=benchmark.candidate.metadata),
        performance=nothing,checksum_verified=missing,workload_verified=missing)
    reference_directory=dirname(benchmark.reference.metadata.path)
    candidate_directory=dirname(benchmark.candidate.metadata.path)
    if dirname(reference_directory)==dirname(candidate_directory)
        path=joinpath(dirname(reference_directory),"performance.jld2")
        if isfile(path)
            calculations=(reference=get(benchmark.reference.metadata,:calculation,nothing),
                candidate=get(benchmark.candidate.metadata,:calculation,nothing))
            measurements=merge(measurements,read_benchmark(path,Val(:performance);calculations))
        end
    end
    return merge(benchmark,(;measurements))
end

"""Read optional timing evidence with its checksum, original session and exact workloads."""
function read_benchmark(path::AbstractString,::Val{:performance};calculations)
    bound=isfile(path*".sha256")
    if bound
        checksum=split(read(path*".sha256",String))
        !isempty(checksum) && first(checksum)==bytes2hex(open(sha256,path)) ||
            throw(ArgumentError("performance checksum mismatch"))
    end
    retained=JLD2.load(path)
    get(retained,"schema_version",1)==1 || throw(ArgumentError("unsupported performance record schema"))
    performance=retained["performance"]
    verified=Bool[]
    if performance!==nothing
        for role in (:reference,:candidate)
            calculation=getproperty(calculations,role)
            calculation===nothing && (push!(verified,false);continue)
            expected=haskey(calculation,:id) ? _numerical_record(calculation) : calculation
            getproperty(performance,role).calculation==expected ||
                throw(ArgumentError("performance workload differs from selected $role calculation"))
            push!(verified,true)
        end
    end
    return (performance,checksum_verified=bound ? true : missing,
        workload_verified=length(verified)==2 && all(verified) ? true : missing,
        session=get(retained,"session",nothing))
end

"""
Compare explicit file bindings in a TOML benchmark definition; missing operands are errors, never inferred.
"""
function compare_saved(definition::AbstractString; directory::AbstractString)
    source = abspath(definition)
    plan = TOML.parsefile(source)
    get(plan, "schema_version", nothing) == 1 ||
        throw(ArgumentError("unsupported benchmark definition schema"))
    entries = get(plan, "benchmarks", nothing)
    entries isa Vector && !isempty(entries) ||
        throw(ArgumentError("define at least one explicit benchmark"))
    allunique(entry["id"] for entry in entries) ||
        throw(ArgumentError("duplicate benchmark IDs"))
    # Resolve every declared operand before publishing comparisons.
    benchmarks = map(entries) do entry
        calculations = map((:reference, :candidate)) do role
            haskey(entry, string(role)) ||
                throw(ArgumentError("benchmark $(entry["id"]) needs an explicit $role"))
            binding = entry[string(role)]
            value = read_calculation(normpath(joinpath(dirname(source), binding["path"]));
                sha256_expected = binding["sha256"])
            BenchmarkCalculation(role, value, value.metadata.formulation)
        end
        settings = get(entry, "comparison", get(plan, "comparison", Dict{String, Any}()))
        kind=get(settings, "kind", "line_parameters")
        kind in ("line_parameters", "uq_moments") || throw(ArgumentError("unknown comparison kind"))
        isempty(setdiff(keys(settings), ("kind","quantities","statistics","bands","normalizations",
            "atol","fundamental","harmonics","unsupported","pairing"))) ||
            throw(ArgumentError("unknown comparison settings"))
        # Convert only declared file values. ReportBuilder owns all omitted defaults.
        comparison=(; (Symbol(key)=>
            (key in ("quantities","statistics","normalizations") ? Tuple(Symbol.(value)) :
             key == "bands" ? Tuple(band isa String ? Symbol(band) : Tuple(band) for band in value) :
             key == "pairing" ? Tuple(Tuple(pair) for pair in value) :
             key in ("atol","unsupported") && value isa AbstractDict ? (;(Symbol(k)=>v for (k,v) in value)...) : value)
            for (key,value) in settings if key != "kind")...)
        kind == "uq_moments" && !haskey(comparison,:statistics) &&
            (comparison=merge(comparison,(statistics=(:mean,:std),)))
        model = (id = Symbol(entry["case"]),
            description = get(entry, "description", entry["case"]))
        benchmark_definition(
            Symbol(entry["id"]), model.id, Symbol(get(plan, "collection", "manual")),
            source, model, calculations..., comparison, get(entry, "tolerances", (;)))
    end
    foreach(validate, benchmarks)
    return [compare_saved(benchmark; directory) for benchmark in benchmarks]
end

"""Hash the complete retained operands, comparisons and declaration of one benchmark."""
function semantic_sha256(::typeof(read_benchmark),directory::AbstractString)
    root=abspath(directory)
    read_benchmark(root)
    state_path=joinpath(root,"state.toml")
    if isfile(state_path)
        state=TOML.parsefile(state_path)
        haskey(state,"current") && (root=joinpath(root,state["current"]))
    end
    files=[joinpath(root,role,"calculation.jld2") for role in ("reference","candidate")]
    declaration=joinpath(root,"declarations.jld2")
    isfile(declaration) && push!(files,declaration)
    for (folder,_,names) in walkdir(joinpath(root,"analyses"))
        "snapshot.jld2" in names && push!(files,joinpath(folder,"snapshot.jld2"))
    end
    return semantic_sha256(Dict(relpath(path,root)=>bytes2hex(open(sha256,path)) for path in files))
end
