using DataFrames
using JLD2
import TOML
import Pkg.Artifacts
import LineCableModels as LCM
using .Gauntlet

# Gauntlet chooses retained files; ReportBuilder owns grouping and tabulation.
function gauntlet_results(directory::AbstractString)
    root=abspath(directory)
    isfile(joinpath(root,"bundle.toml")) && return read_campaign(root)
    if isfile(joinpath(root,"release.toml"))
        declaration=TOML.parsefile(joinpath(root,"release.toml"))
        return read_collection(root;collection=Symbol(declaration["collection"]))
    end
    if isfile(joinpath(root,"campaign.toml"))
        return [read_benchmark(joinpath(root,row.id)) for row in campaign_status(root) if row.state === :complete]
    end
    if isfile(root)
        return [read_benchmark(root;load_results=true)]
    end
    paths=sort!([joinpath(folder,"snapshot.jld2") for (folder,_,files) in walkdir(root) if "snapshot.jld2" in files])
    isempty(paths) && throw(ArgumentError("no saved benchmarks; references cannot be inferred from calculations"))
    return [read_benchmark(path;load_results=true) for path in paths]
end

"""Render only saved summaries and explicitly retained illustrations; never calculate or plot."""
function render_gauntlet_report(source)
    io=IOBuffer()
    roots=String[]
    if source === nothing
        document=TOML.parsefile(joinpath(@__DIR__,"gauntlet.toml"))
        selected=document["artifacts"]
        isempty(selected) && return """
        No published benchmark artifacts are selected. Add immutable version bindings
        to `docs/gauntlet.toml` after `lcm gauntlet package`, upload and `lcm gauntlet bind`.
        An explicit `LINECABLEMODELS_GAUNTLET_RESULTS` directory enables a local draft preview.
        """
        for name in selected
            bindings=TOML.parsefile(Gauntlet.ARTIFACTS_TOML)
            haskey(bindings,name) || throw(ArgumentError("published artifact binding is missing: $name"))
            occursin(r"_v[0-9]+_[0-9]+_[0-9]+$",name) || throw(ArgumentError("documentation requires version-specific artifact bindings"))
            Artifacts.ensure_artifact_installed(name,Gauntlet.ARTIFACTS_TOML)
            hash=Artifacts.artifact_hash(name,Gauntlet.ARTIFACTS_TOML)
            push!(roots,Artifacts.artifact_path(hash))
            println(io,"Published artifact `",name,"` · tree `",hash,"`.\n")
            release=TOML.parsefile(joinpath(Artifacts.artifact_path(hash),"release.toml"))
            println(io,"Accepted snapshots: ",join(["`$id`" for id in release["bundles"]],", "),".\n")
        end
    else
        append!(roots,unique(abspath.(split(source,Sys.iswindows() ? ';' : ':'))))
        println(io,"**Local preview of explicitly selected results.**\n")
    end
    records=reduce(vcat,gauntlet_results.(roots);init=Any[])
    statuses=NamedTuple[]
    for root in roots
        isfile(joinpath(root,"campaign.toml")) || continue
        append!(statuses,campaign_status(root))
        if isfile(joinpath(root,"bundle.toml"))
            println(io,"Accepted snapshot `",TOML.parsefile(joinpath(root,"bundle.toml"))["identity"],"`.\n")
        end
    end
    println(io,length(records)," complete benchmarks across ",length(unique(first(record.analyses)["case_id"] for record in records))," cases.\n")
    if any(row -> row.state !== :complete,statuses)
        println(io,"Incomplete drafts: ",join(["$(row.id): $(row.state)" for row in statuses if row.state !== :complete],", "),".\n")
    end
    maxima=NamedTuple[]
    labels=NamedTuple[]
    for record in records,analysis in record.analyses
        # Schema 1 is a finite historical reader; it uses saved RMS arrays only.
        if haskey(analysis,"summary")
            append!(maxima,analysis["summary"])
            append!(labels,[(benchmark=analysis["benchmark_id"],role=row.role,index=row.formulation_index,label=row.label)
                for row in analysis["formulations"]])
        else
            loaded=(id=record.id,reference=record.reference,candidate=record.candidate,analyses=[analysis])
            table=LCM.ReportBuilder.report(LCM.ReportBuilder.BenchmarkTableDefinition(),loaded).table
            append!(maxima,NamedTuple.(eachrow(table.maxima)))
            append!(labels,[(benchmark=analysis["benchmark_id"],role=row.role,index=row.formulation_index,label=row.label)
                for row in eachrow(table.formulations)])
        end
    end
    definition=LCM.ReportBuilder.BenchmarkTableDefinition()
    summary=LCM.ReportBuilder.tabulate(definition,nothing,maxima)
    println(io,"Entries show **maxima of per-term RMS discrepancies**, with the terminal pair and unavailable-term count.\n")
    for band in unique(vcat(collect(definition.settings.bands),[row.band for row in maxima]))
        selected=isempty(summary) ? summary : filter(row -> row.band == band,summary)
        isempty(selected) && continue
        println(io,"### ",LCM.description(definition,band),"\n")
        bounds=unique(selected[:,[:samples,:requested_bounds_Hz,:actual_bounds_Hz]])
        bounds[!,:selection]=collect(1:nrow(bounds))
        selected[!,:sample_selection]=[findfirst(other -> isequal(
            (other.samples,other.requested_bounds_Hz,other.actual_bounds_Hz),
            (row.samples,row.requested_bounds_Hz,row.actual_bounds_Hz)),eachrow(bounds)) for row in eachrow(selected)]
        columns=[:snapshot,:case_id,:benchmark,:problem_index,:formulation_index,:reference_point,
            :statistic,:normalization,:sample_selection]
        append!(columns,[q for q in (:Z,:Y,:R,:L,:G,:C) if q in propertynames(selected)])
        compact=selected[:,columns]
        compact.snapshot=[length(string(id))==64 ? first(string(id),12) : string(id) for id in compact.snapshot]
        println(io,"Stored samples (Hz):\n\n```@raw html")
        show(IOContext(io,:limit=>false),MIME"text/html"(),bounds;summary=false,eltypes=false)
        println(io,"\n```\n\n```@raw html")
        show(IOContext(io,:limit=>false),MIME"text/html"(),compact;summary=false,eltypes=false)
        println(io,"\n```\n")
    end
    if !isempty(labels)
        println(io,"<details><summary>Formulation keys</summary>\n\n```@raw html")
        show(IOContext(io,:limit=>false),MIME"text/html"(),DataFrame(unique(labels));summary=false,eltypes=false)
        println(io,"\n```\n</details>\n")
    end
    println(io,"Reopen an accepted bundle with `read_campaign(path)`, select a benchmark, then use `report(BenchmarkTableDefinition(), benchmark)` or explicitly `plot(benchmark, (Z, Y))`. Detailed plots are generated only by that request.\n")
    # Explicit illustrations are already files; rebuilding the page never renders them.
    for root in roots
        bundles=isfile(joinpath(root,"release.toml")) ?
            [joinpath(root,"bundles",id) for id in TOML.parsefile(joinpath(root,"release.toml"))["bundles"]] : [root]
        for bundle in bundles
            isfile(joinpath(bundle,"bundle.toml")) || continue
            document=TOML.parsefile(joinpath(bundle,"bundle.toml"))
            for figure in get(document,"illustrations",[])
                source_path=joinpath(bundle,figure["path"])
                destination=joinpath(@__DIR__,"src","assets","gauntlet",document["identity"],basename(source_path))
                mkpath(dirname(destination))
                cp(source_path,destination;force=true)
                caption=replace(figure["caption"],"["=>"\\[","]"=>"\\]")
                println(io,"![",caption,"](assets/gauntlet/",document["identity"],"/",basename(source_path),")\n")
            end
        end
    end
    return String(take!(io))
end
