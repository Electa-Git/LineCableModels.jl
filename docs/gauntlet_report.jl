using DataFrames
using JLD2
using SHA
import TOML
import LineCableModels as LCM
include(joinpath(@__DIR__, "..", "test", "gauntlet", "fingerprints.jl"))

# Documentation reads completed defaults, never case builders or solver modules.
function gauntlet_defaults(directory::AbstractString)
    root = abspath(directory)
    plan = TOML.parsefile(joinpath(root, "campaign.toml"))
    get(plan, "schema_version", nothing) == 1 || throw(ArgumentError(
        "unsupported Gauntlet campaign schema: $root"))
    records = NamedTuple[]
    pending = NamedTuple[]
    for job in plan["jobs"]
        get(job, "propagation", "deterministic") == "deterministic" || continue
        get(plan, "dielectric", "default") == "default" || continue
        occursin(r"^[a-z][a-z0-9_]*$", job["id"]) || throw(ArgumentError("invalid campaign job ID"))
        for (index, selection) in enumerate(job["selections"])
            get(selection, "id", "") == "default" || length(selection) > 1 || continue
            all(name == "id" || value == "default" for (name, value) in selection) || continue
            path = joinpath(root, job["id"], lpad(string(index), 4, '0') * ".jld2")
            if !isfile(path) || !isfile(path * ".sha256")
                push!(pending, (case=job["case"], backend=job["backend"]))
                continue
            end
            digest = bytes2hex(open(sha256, path))
            checksum = split(read(path * ".sha256", String))
            !isempty(checksum) && first(checksum) == digest || throw(ArgumentError(
                "Gauntlet calculation checksum mismatch: $path"))
            document = jldopen(path, "r") do file
                Dict(key => file[key] for key in ("kind", "schema_version", "status",
                    "selection", "case_id", "backend", "formulation", "domain",
                    "Z", "Y", "frequencies", "basis", "port_order", "problem"))
            end
            get(document, "kind", nothing) === :gauntlet_calculation &&
                get(document, "schema_version", nothing) == 1 &&
                get(document, "status", nothing) === :complete || throw(ArgumentError(
                    "unsupported or incomplete Gauntlet calculation: $path"))
            document["selection"] == selection &&
                string(document["case_id"]) == job["case"] &&
                string(document["backend"]) == job["backend"] || throw(ArgumentError(
                    "calculation does not match its campaign selection: $path"))
            all(value -> (value isa Symbol ? value : LCM.formula_id(value)) === :default,
                document["formulation"].definitions) || continue
            document["domain"] === :PhaseDomain || throw(ArgumentError(
                "Gauntlet summary requires phase-domain results: $path"))
            parameters = LCM.LineParameters(LCM.PhaseDomain, document["Z"], document["Y"],
                document["frequencies"]; basis=document["basis"])
            length(document["port_order"]) == size(document["Z"], 1) || throw(DimensionMismatch(
                "stored terminal order does not match the result matrices: $path"))
            identity = semantic_sha256(document["problem"])
            push!(records, (; case=job["case"], backend=job["backend"],
                title=get(job, "description", job["case"]), path, digest, identity,
                ports=document["port_order"], parameters))
        end
    end
    return (; records, pending)
end

function gauntlet_comparisons(records)
    rows = NamedTuple[]
    for case in sort!(unique(record.case for record in records))
        selected = sort!(filter(record -> record.case == case, records); by=record -> record.backend)
        analytical = findfirst(record -> record.backend == "coaxial", selected)
        analytical === nothing && continue
        reference = selected[analytical]
        for candidate in selected
            candidate.backend == "coaxial" && continue
            reference.identity == candidate.identity && reference.ports == candidate.ports &&
                LCM.basis(reference.parameters) === LCM.basis(candidate.parameters) &&
                LCM.frequencies(reference.parameters) == LCM.frequencies(candidate.parameters) ||
                throw(ArgumentError("stored default inputs or coordinates differ for $case; " *
                    "select matching campaigns, without interpolation or relabelling"))
            comparison = LCM.Engine.compare(reference.parameters, candidate.parameters)
            push!(rows, (;
                case=reference.title, baseline=reference.backend, backend=candidate.backend,
                max_z_percent=100maximum(LCM.observe(comparison, LCM.Z, LCM.Engine.relative_error)),
                max_y_percent=100maximum(LCM.observe(comparison, LCM.Y, LCM.Engine.relative_error))))
        end
    end
    return rows
end

function render_gauntlet_report(source)
    source === nothing && return """
    !!! note "No recorded results selected"
        Set `LINECABLEMODELS_GAUNTLET_RESULTS` to a campaign directory to show
        completed default-formulation comparisons. No calculations run during
        documentation generation.
    """
    directories = unique(abspath.(split(source, Sys.iswindows() ? ';' : ':')))
    records = NamedTuple[]
    pending = NamedTuple[]
    for directory in directories
        loaded = gauntlet_defaults(directory)
        append!(pending, loaded.pending)
        for record in loaded.records
            previous = findfirst(value -> value.case == record.case && value.backend == record.backend, records)
            if previous !== nothing
                records[previous].digest == record.digest || throw(ArgumentError(
                    "multiple stored defaults for $(record.case)/$(record.backend); select one campaign"))
                continue
            end
            push!(records, record)
        end
    end
    comparisons = gauntlet_comparisons(records)
    io = IOBuffer()
    println(io, length(records), " completed defaults across ",
        length(unique(record.case for record in records)), " cases; ",
        length(pending), " default selections unfinished.\n")
    if isempty(comparisons)
        println(io, "No matching completed backend defaults are available for comparison.")
    else
        frame = DataFrame(:case=>getproperty.(comparisons, :case),
            :baseline=>getproperty.(comparisons, :baseline),
            :backend=>getproperty.(comparisons, :backend),
            Symbol("max εZ")=>getproperty.(comparisons, :max_z_percent),
            Symbol("max εY")=>getproperty.(comparisons, :max_y_percent))
        println(io, "```@raw html")
        show(IOContext(io, :limit=>false), MIME"text/html"(), frame)
        println(io, "\n```")
    end
    return String(take!(io))
end
