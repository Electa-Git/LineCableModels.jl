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
    println(io,"Entries are maxima of per-term RMS discrepancies. Relative tables use percent; absolute tables retain physical units.\n")
    for record in records, analysis in record.analyses
        selected=merge(record,(analyses=[analysis],))
        artifact=LCM.ReportBuilder.report(LCM.ReportBuilder.BenchmarkTableDefinition(),selected)
        tables=artifact.table
        println(io,"## ",record.id,"\n")
        for feature in tables.features
            println(io,"### ",feature.quantity," · ",feature.statistic===:std ? "standard deviation" : string(feature.statistic),
                " · point ",feature.problem_index," · ",feature.normalization,"\n")
            println(io,"Analysis: `",feature.snapshot,"`.\n")
            for (label,frame) in (("Relative RMS [%]",feature.relative),
                    ("Absolute RMS ["*feature.absolute_unit*"]",feature.absolute))
                println(io,label,"\n\n```@raw html")
                show(IOContext(io,:limit=>false),MIME"text/html"(),frame;summary=false,eltypes=false)
                println(io,"\n```\n")
            end
        end
        for (label,frame) in (
                "Worst terms and comparison counts"=>tables.maxima,
                "Formulations"=>tables.formulations,
                "Scientific formula descriptions"=>tables.formula_details,
                "Retained UQ statistics"=>tables.statistics,
                "MC sampling precision"=>tables.sampling,
                "MC mean sampling errors by term and frequency"=>tables.mean_sampling_precision,
                "Execution timings"=>tables.execution,
                "Native solver timing scopes"=>tables.source_timings,
                "Controlled performance measurements"=>tables.performance,
                "Reference/candidate timing ratio"=>tables.performance_comparison)
            isempty(frame) && continue
            println(io,"<details><summary>",label," (",nrow(frame)," rows)</summary>\n\n```@raw html")
            show(IOContext(io,:limit=>true,:displaysize=>(30,150)),MIME"text/html"(),frame;eltypes=false)
            println(io,"\n```\n</details>\n")
        end
    end
    println(io,"Historical RMS is rendered unchanged. For current comparisons, explicitly report the saved reference and candidate operands with a new definition. This requires no solver run.\n")
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
