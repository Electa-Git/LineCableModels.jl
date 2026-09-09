using DataFrames
using JLD2
using SHA
import LineCableModels as LCM
using .Gauntlet

# Rendering consumes persisted comparisons, not independent calculations. It
# does not select operands, recalculate errors, or load backend implementations.
function gauntlet_results(directory::AbstractString)
    root = abspath(directory)
    isdir(root) || throw(ArgumentError("benchmark directory is missing: $root"))
    paths = sort!([joinpath(folder, "snapshot.jld2")
                   for (folder, _, files) in walkdir(root)
                   if "snapshot.jld2" in files])
    isempty(paths) && throw(ArgumentError("no explicit benchmark records in $root; " *
                        "use lcm gauntlet compare --definition FILE.toml --output DIR; references cannot be inferred from calculations"))
    return read_benchmark.(paths)
end

function gauntlet_comparisons(record)
    rows = NamedTuple[]
    for result in record["reference_comparison"]
        indices = findall(!ismissing, result.relative)
        maximum_index = isempty(indices) ? nothing :
                        indices[argmax(result.relative[indices])]
        entry = maximum_index === nothing ? "—" : string(Tuple(maximum_index))
        bounds = result.details.actual_bounds
        push!(rows,
            (quantity = string(result.quantity), statistic = string(result.statistic),
                reference_index = get(result, :reference_index, 1), candidate_index = get(
                    result, :candidate_index, 1),
                band = string(result.details.band), samples = result.details.sample_count,
                Hz = ismissing(first(bounds)) ? "—" :
                     join(round.(bounds; sigdigits = 5), "–"),
                normalization = result.details.normalization,
                percent = maximum_index === nothing ? missing :
                          100result.relative[maximum_index], entry,
                zeros = count(==(:below_tolerance), result.details.status),
                unavailable = count(ismissing, result.relative),
                reason = result.details.reason === nothing ?
                         join(
                    unique(filter(!isnothing,
                        vec(get(result.details,
                            :normalization_reason, fill(nothing, size(result.relative)))))),
                    "; ") :
                         result.details.reason))
    end
    return rows
end

function gauntlet_table(
        records, band, selections; quantities = nothing, statistic = :value)
    available = [row for record in records for row in gauntlet_comparisons(record)
                 if row.statistic == string(statistic)]
    quantities === nothing && (quantities = unique(Symbol(row.quantity) for row in available))
    normalizations = unique(row.normalization for row in available)
    frame = DataFrame()
    for record in records,
        point in unique((row.reference_index, row.candidate_index)
        for row in gauntlet_comparisons(record))

        operands = record["calculations"]
        row = Dict{Symbol, Any}(
            :reference => "$(operands.reference.backend) [$(selections[operands.reference.selection])]",
            :candidate => "$(operands.candidate.backend) [$(selections[operands.candidate.selection])]",
            :reference_point=>point[1], :candidate_point=>point[2])
        comparisons = gauntlet_comparisons(record)
        for quantity in quantities, normalization in normalizations
            label = normalization === :reference_rms ? "NRMSE" : "pointwise"
            metrics = filter(
                r -> r.quantity == string(quantity) &&
                     r.statistic == string(statistic) && r.band == band &&
                     r.normalization === normalization &&
                     (r.reference_index, r.candidate_index)==point,
                comparisons)
            # Missing metrics stay missing; rendering never manufactures a new comparison.
            cell = if isempty(metrics)
                "missing (not recorded)"
            else
                metric = only(metrics)
                if metric.samples == 0
                    "missing (no samples)"
                elseif ismissing(metric.percent)
                    "missing ($(metric.reason))"
                else
                    value = string(round(metric.percent; sigdigits = 4), " ", metric.entry)
                    metric.unavailable == 0 ? value :
                    "$value; $(metric.unavailable) unavailable"
                end
            end
            row[Symbol("$quantity $label")] = cell
        end
        push!(frame, row; cols = :union)
    end
    columns = [Symbol("$quantity $(normalization === :reference_rms ? "NRMSE" : "pointwise")")
               for quantity in quantities for normalization in normalizations]
    return select!(
        frame, :reference, :candidate, :reference_point, :candidate_point, columns...)
end

function render_gauntlet_report(source)
    source === nothing && return """
    !!! note "No recorded benchmarks selected"
        Set `LINECABLEMODELS_GAUNTLET_RESULTS` to a saved benchmark directory.
        Use `lcm gauntlet compare` to bind completed calculations explicitly.
        No calculations run during documentation generation.
    """
    directories = unique(abspath.(split(source, Sys.iswindows() ? ';' : ':')))
    records = reduce(vcat, gauntlet_results.(directories))
    unique_records = Dict{Tuple, Any}()
    for record in records
        key = (
            record["collection"], record["benchmark_id"], repr(record["comparison_settings"]),
            record["calculations"].reference.sha256, record["calculations"].candidate.sha256)
        haskey(unique_records, key) && !isequal(unique_records[key], record) &&
            throw(ArgumentError("conflicting saved benchmark: $key"))
        unique_records[key] = record
    end
    records = sort!(collect(values(unique_records));
        by = r->(
            r["case_id"], r["benchmark_id"], repr(r["comparison_settings"])))
    io = IOBuffer()
    println(io, length(records), " explicitly configured benchmarks across ",
        length(unique(r["case_id"] for r in records)), " cases. ",
        "Only completed, checksum-verified operands are included; no backend is selected as a reference by this page.\n")
    selections = Dict{Any, Int}()
    captions = String[]
    for record in records, operand in record["calculations"]

        haskey(selections, operand.selection) && continue
        caption = repr(operand.formulation)
        caption in captions || push!(captions, caption)
        selections[operand.selection] = findfirst(==(caption), captions)
    end
    println(io,
        "Bracketed numbers identify calculation declarations; their complete formulations and controls are retained in the artifacts:\n")
    println(io, "```@raw html")
    show(IOContext(io, :limit=>false), MIME"text/html"(),
        DataFrame(:selection=>eachindex(captions), :formulation=>captions);
        summary = false, eltypes = false)
    println(io, "\n```\n")

    # Escape text in the few structural HTML elements. DataFrames owns cell escaping.
    escape = value -> replace(string(value), '&'=>"&amp;", '<'=>"&lt;", '>'=>"&gt;",
        '"'=>"&quot;", '\''=>"&#39;")
    for (kind, statistic, title) in
        ((:line_parameters, :value, "Full-band comparisons"),
        (:uq_moments, :mean, "UQ means"),
        (:uq_moments, :std, "UQ standard deviations (std)"))
        group = filter(r -> statistic in r["comparison_settings"].statistics, records)
        isempty(group) && continue
        visible = [row for record in group
                   for row in gauntlet_comparisons(record)
                   if row.statistic == string(statistic)]
        quantities = unique(Symbol(row.quantity) for row in visible)
        bands = unique(row.band for row in visible)
        for (band_index, band) in enumerate(bands)
            if band == "all"
                println(io, "### ", title, "\n")
            else
                band_index == 2 && println(io, "### Frequency slices — ",
                    string(statistic), "\n")
                println(io, "#### Band `", band, "`\n")
            end
            samples = unique((row.Hz, row.samples)
            for row in visible
            if row.band == band && row.samples > 0)
            if isempty(samples)
                println(io, "No stored samples or comparisons for this band; errors are missing.\n")
                continue
            end
            if length(samples) == 1
                bounds, count = only(samples)
                println(io, "Stored range: **", bounds, " Hz**, **", count, " samples**.\n")
            end
            for case in unique(r["case_id"] for r in group)
                selected = filter(r -> r["case_id"] == case, group)
                println(io, "```@raw html\n<h4>",
                    escape(first(selected)["description"]), "</h4>")
                frame = gauntlet_table(selected, band, selections; quantities, statistic)
                show(IOContext(io, :limit=>false), MIME"text/html"(), frame;
                    summary = false, eltypes = false)
                if band == "all" && statistic !== :std
                    println(io, "\n<details><summary>Benchmark identities and terminal order</summary>")
                    println(io, "<p>Case: <code>", escape(case), "</code>.</p><ol>")
                    for record in selected
                        println(io, "<li><code>", escape(record["collection"]), "/",
                            escape(record["benchmark_id"]), "</code></li>")
                    end
                    println(io, "</ol>")
                    for ports in unique(r["port_order"] for r in selected)
                        indices = findall(r -> r["port_order"] == ports, selected)
                        println(io, "<p>Terminal order (benchmark rows ",
                            join(indices, ", "), "): ",
                            join(
                                ("<code>$(index)=$(escape(port))</code>"
                                for (index, port) in enumerate(ports)),
                                ", "), ".</p>")
                    end
                    println(io, "</details>")
                end
                unavailable = [(index, row.quantity, row.reason)
                               for (index, record) in enumerate(selected)
                               for row in gauntlet_comparisons(record)
                               if row.band == band && Symbol(row.quantity) in quantities &&
                                  row.statistic == string(statistic) && row.unavailable > 0]
                for (index, quantity, reason) in unique(unavailable)
                    println(io, "<p>Benchmark row ", index, ", ", quantity,
                        " unavailable: ", escape(reason), ".</p>")
                end
                if length(samples) > 1
                    println(io, "\n<p>Stored ranges by benchmark row: ",
                        join(
                            (string(index,
                                 ": ",
                                 join(
                                     unique(
                                         "$(row.Hz) Hz ($(row.samples) samples)"
                                     for row in gauntlet_comparisons(record)
                                     if row.band == band &&
                                        Symbol(row.quantity) in quantities &&
                                        row.statistic == string(statistic)),
                                     ", "))
                            for (index, record) in enumerate(selected)),
                            "; "), ".</p>")
                end
                println(io, "\n```\n")
            end
        end
    end
    timings=[(benchmark = record["benchmark_id"], operand = string(role),
                 scope = value.scope, seconds = value.seconds)
             for record in records for (role, value) in pairs(record["timings"])]
    println(io,
        "### Recorded computation timings\n\nEach timing retains its declared scope. Legacy batch timings are not per-formula measurements.\n")
    println(io, "```@raw html")
    show(IOContext(io, :limit=>false), MIME"text/html"(),
        DataFrame(timings); summary = false, eltypes = false)
    println(io, "\n```\n")
    return String(take!(io))
end
