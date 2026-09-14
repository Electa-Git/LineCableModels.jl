_report_name(value) = String(nameof(typeof(value)))

function Base.summary(io::IO, definition::AbstractReportDefinition)
    print(io, _report_name(definition), " report definition")
end
function Base.show(io::IO, definition::AbstractReportDefinition)
    print(io, _report_name(definition), "()")
end
Base.show(io::IO, ::MIME"text/plain", definition::AbstractReportDefinition) =
    show(io, definition)

Base.summary(io::IO, definition::TableReportDefinition) =
    print(io, "Table report with $(length(definition.requests)) observations")
function Base.show(io::IO, definition::TableReportDefinition)
    print(io, "TableReportDefinition(requests=", length(definition.requests), ")")
end
function Base.show(io::IO, ::MIME"text/plain", definition::TableReportDefinition)
    get(io, :compact, false) && return show(io, definition)
    children = Any[
        (label = "requests      $(length(definition.requests))", noun = "fields"),
        (label = "clip          $(definition.clip)", noun = "fields"),
    ]
    definition.illustration === nothing || push!(children,
        (label = "illustration  $(definition.illustration)", noun = "fields"))
    return TextDisplay.tree(io, "Table report definition", Tuple(children))
end

Base.summary(io::IO, ::CableConstantsTableDefinition) =
    print(io, "Cable-constants table definition")
Base.show(io::IO, definition::CableConstantsTableDefinition) =
    print(io, "CableConstantsTableDefinition(clip=", definition.clip, ")")
Base.show(io::IO, ::MIME"text/plain", definition::CableConstantsTableDefinition) =
    show(io, definition)

Base.summary(io::IO, definition::LineParametersTableDefinition) =
    print(io, "Line-parameters table with $(length(definition.requests)) observations")
function Base.show(io::IO, definition::LineParametersTableDefinition)
    print(io, "LineParametersTableDefinition(requests=", length(definition.requests),
        ", length_unit=:", definition.length_unit, ")")
end
Base.show(io::IO, ::MIME"text/plain", definition::LineParametersTableDefinition) =
    show(io, definition)

Base.summary(io::IO, ::BenchmarkTableDefinition) = print(io, "Benchmark table definition")
Base.show(io::IO, definition::BenchmarkTableDefinition) =
    print(io, "BenchmarkTableDefinition(clip=", definition.clip, ")")
Base.show(io::IO, ::MIME"text/plain", definition::BenchmarkTableDefinition) =
    show(io, definition)

Base.summary(io::IO, ::MonteCarloTableDefinition) = print(io, "Monte Carlo table definition")
function Base.show(io::IO, definition::MonteCarloTableDefinition)
    print(io, "MonteCarloTableDefinition(length_unit=:", definition.length_unit,
        ", clip=", definition.clip, ")")
end
Base.show(io::IO, ::MIME"text/plain", definition::MonteCarloTableDefinition) =
    show(io, definition)

Base.summary(io::IO, ::XLSXReportDefinition) = print(io, "XLSX report definition")
function Base.show(io::IO, definition::XLSXReportDefinition)
    destination = definition.file_name === nothing ? "default path" : repr(definition.file_name)
    print(io, "XLSXReportDefinition(", destination, "; clip=", definition.clip, ")")
end
Base.show(io::IO, ::MIME"text/plain", definition::XLSXReportDefinition) =
    show(io, definition)

Base.summary(io::IO, artifact::ReportArtifact) = print(io, "Completed report artifact")
function Base.show(io::IO, artifact::ReportArtifact)
    dimensions = applicable(size, artifact.table) ? join(size(artifact.table), '×') : "one table"
    print(io, "ReportArtifact(table=", dimensions,
        ", illustration=", artifact.illustration === nothing ? "none" : "present",
        ", output=", artifact.output === nothing ? "none" : "present", ")")
end
function Base.show(io::IO, ::MIME"text/plain", artifact::ReportArtifact)
    get(io, :compact, false) && return show(io, artifact)
    dimensions = applicable(size, artifact.table) ? join(size(artifact.table), '×') : "one table"
    return TextDisplay.tree(io, "Report artifact", (
        (label = "table         $dimensions", noun = "fields"),
        (label = "illustration  $(artifact.illustration === nothing ? "none" : "present")", noun = "fields"),
        (label = "output        $(artifact.output === nothing ? "none" : "present")", noun = "fields"),
    ))
end

Base.summary(io::IO, sheet::XLSXSheet) = print(io, "XLSX sheet \"", sheet.name, "\"")
function Base.show(io::IO, sheet::XLSXSheet)
    print(io, "XLSXSheet(\"", sheet.name, "\"; cells=", join(size(sheet.cells), '×'), ")")
end
Base.show(io::IO, ::MIME"text/plain", sheet::XLSXSheet) = show(io, sheet)

Base.summary(io::IO, workbook::XLSXWorkbook) =
    print(io, "XLSX workbook with $(length(workbook.sheets)) sheets")
function Base.show(io::IO, workbook::XLSXWorkbook)
    print(io, "XLSXWorkbook(", repr(workbook.destination), "; sheets=",
        length(workbook.sheets), ")")
end
function Base.show(io::IO, ::MIME"text/plain", workbook::XLSXWorkbook)
    get(io, :compact, false) && return show(io, workbook)
    sheets = Tuple((label = "$(sheet.name) · $(join(size(sheet.cells), '×')) cells",
        noun = "sheets") for sheet in workbook.sheets)
    return TextDisplay.tree(io, "XLSX workbook · $(length(workbook.sheets)) sheets", sheets;
        noun = "sheets")
end

"""
$(TYPEDSIGNATURES)

Display a benchmark's compact tables, without rendering plots or recalculating
comparisons. Each quantity and statistic has its own formula-by-band table.
`metric=:relative` shows RMS percentages; `:absolute` shows native physical units.
`problem` optionally selects a parameter configuration, not a frequency or MC trial.
Recorded workload timings remain whole-call measurements even when a configuration
is selected. `native_timings=true` also displays backend-owned timing scopes.

The same tables are rendered in plain text and HTML. Full term coordinates,
availability reasons, statistics and timing observations remain in `artifact.table`.
"""
function Base.show(io::IO,mime::MIME"text/plain",artifact::ReportArtifact{P};kwargs...) where
        {P<:NamedTuple{(:reference,:candidate,:context,:settings,:comparisons,:measurements)}}
    return _show_benchmark(io,mime,artifact;kwargs...)
end
function Base.show(io::IO,mime::MIME"text/html",artifact::ReportArtifact{P};kwargs...) where
        {P<:NamedTuple{(:reference,:candidate,:context,:settings,:comparisons,:measurements)}}
    return _show_benchmark(io,mime,artifact;kwargs...)
end

# Both MIME methods share table selection and labels. Separate entry methods
# avoid ambiguity with the ordinary ReportArtifact text/plain display.
function _show_benchmark(io::IO,mime,artifact::ReportArtifact;
        metric::Symbol=:relative,problem=nothing,native_timings::Bool=false)
    get(io,:compact,false) && return show(io,artifact)
    metric in (:relative,:absolute) || throw(ArgumentError("metric must be :relative or :absolute"))
    tables=artifact.table
    features=filter(feature -> problem===nothing || feature.problem_index==problem,tables.features)
    isempty(features) && throw(ArgumentError("Selected parameter configuration is not present in this report"))
    multiple_points=length(unique(feature.problem_index for feature in features))>1
    show_reference_point=length(unique(feature.reference_point for feature in features))>1 ||
        any(feature -> feature.reference_point!=feature.problem_index,features)
    multiple_snapshots=length(unique(feature.snapshot for feature in features))>1
    sections=Pair{String,DataFrame}[]
    # Quantity groups stay together, including distinct mean/std or other retained
    # statistics. The original feature matrices and their missing masks are reused.
    for quantity in unique(feature.quantity for feature in features), feature in features
        feature.quantity==quantity || continue
        statistic=feature.statistic===:value ? "" : feature.statistic===:std ?
            " · std (propagated uncertainty)" : " · $(feature.statistic)"
        title=string(quantity,statistic,
            multiple_points || problem!==nothing ? " · configuration $(feature.problem_index)" : "",
            show_reference_point ? " · reference configuration $(feature.reference_point)" : "",
            multiple_snapshots ? " · analysis $(feature.snapshot)" : "",
            metric===:relative ? " — Relative RMS [%] · $(feature.normalization)" :
                " — Absolute RMS [$(feature.absolute_unit)] · $(feature.normalization)")
        push!(sections,title=>getproperty(feature,metric))
    end
    coverage=filter(row -> problem===nothing || row.point==problem,tables.overview.coverage)
    hidden=[:benchmark,:case_id,:collection,:formulation_index,:candidate_point]
    multiple_snapshots || push!(hidden,:snapshot)
    multiple_points || problem!==nothing || push!(hidden,:point)
    show_reference_point || push!(hidden,:reference_point)
    push!(sections,"Frequency coverage"=>DataFrames.select(coverage,Not(hidden)))
    append!(sections,[
        "Recorded calculation times — whole workloads, not controlled repetitions"=>tables.overview.execution,
        "Controlled calculation measurements — whole workloads"=>tables.overview.performance,
        "Recorded reference / candidate time ratio"=>tables.overview.timing_ratio])
    for (title,frame) in (
            "MC sampling workload"=>tables.overview.sampling,
            "Simultaneous MC CDF precision — probability units"=>tables.overview.cdf_precision)
        isempty(frame) && continue
        # MC point coordinates are source-owned, not candidate/formulation indices.
        # Keep all workloads visible rather than filtering by a comparison index.
        visible=length(unique(frame.point))==1 && only(unique(frame.point))==1 ?
            DataFrames.select(frame,Not(:point)) : frame
        push!(sections,title=>visible)
    end
    native_timings && push!(sections,"Backend timing records — scopes as recorded"=>tables.overview.source_timings)
    notes=[join(tables.formulations.label[tables.formulations.role .=== :reference]," / "),
        "RMS is the maximum per-term RMS in each band, not an average across matrix entries. References are comparison methods, not a truth designation."]
    any(>(0),tables.maxima.unavailable) && push!(notes,
        "Relative maxima use eligible terms only. Missing means no meaningful relative comparison or unavailable data; term counts and reasons are in .table.maxima and .table.terms.")
    isempty(tables.overview.performance) ? push!(notes,"No controlled performance measurements were saved.") :
        push!(notes,"Allocated MiB are cumulative Julia allocations per call, not peak memory. The timing ratio compares the recorded whole workloads, not per-formula costs.")
    !isempty(tables.overview.performance) && any(==(1),tables.overview.performance.timed_calls) &&
        push!(notes,"At least one method has only one timed call; its timing variability cannot be assessed.")
    !isempty(tables.overview.timing_ratio) && any(value -> value!==true,tables.overview.timing_ratio.comparable) &&
        push!(notes,"A ratio marked comparable=false or missing is not a validated speedup comparison.")
    any(feature -> feature.statistic===:std,features) && push!(notes,
        "Std compares propagated physical uncertainty. MC standard errors describe sampling noise in the mean, not physical spread; retained values are in .table.mean_sampling_precision.")
    !isempty(tables.overview.sampling) && push!(notes,
        "The CDF bound is in cumulative-probability units, not relative error in mean or std. MC std sampling precision is not retained.")
    html=mime isa MIME"text/html"
    for note in notes
        if html
            println(io,"<p>",replace(note,'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'\"'=>"&quot;",'\''=>"&#39;"),"</p>")
        else
            println(io,note)
        end
    end
    for (title,frame) in sections
        isempty(frame) && continue
        if html
            println(io,"<h3>",replace(title,'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'\"'=>"&quot;",'\''=>"&#39;"),"</h3>")
            show(IOContext(io,:limit=>false),mime,frame;summary=false,eltypes=false)
        else
            println(io,"\n",title)
            show(io,mime,frame;allrows=true,allcols=true,truncate=0)
        end
        println(io)
    end
    return nothing
end
