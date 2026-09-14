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

"""Render saved compact report tables only; never calculate, time, plot or copy illustrations."""
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
    for record in records, analysis in record.analyses
        selected=merge(record,(analyses=[analysis],))
        artifact=LCM.ReportBuilder.report(LCM.ReportBuilder.BenchmarkTableDefinition(),selected)
        println(io,"## ",record.id,"\n")
        println(io,"Analysis: `",only(unique(feature.snapshot for feature in artifact.table.features)),"`.\n")
        println(io,"```@raw html")
        show(io,MIME"text/html"(),artifact)
        println(io,"\n```\n")
    end
    println(io,"Historical RMS is rendered unchanged. For current comparisons, explicitly report the saved reference and candidate operands with a new definition. This requires no solver run.\n")
    return String(take!(io))
end
