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
    print(io, "XLSXReportDefinition(", destination, "; clip=", definition.clip, ", overwrite=", definition.overwrite, ")")
end
Base.show(io::IO, ::MIME"text/plain", definition::XLSXReportDefinition) =
    show(io, definition)

Base.summary(io::IO, artifact::ReportArtifact) = print(io, "Completed report artifact")
function Base.show(io::IO, artifact::ReportArtifact)
    dimensions = applicable(size, artifact.tables) ? join(size(artifact.tables), '×') : "one table"
    print(io, "ReportArtifact(table=", dimensions,
        ", illustration=", artifact.illustration === nothing ? "none" : "present",
        ", output=", artifact.output === nothing ? "none" : "present", ")")
end
function Base.show(io::IO, mime::MIME"text/plain", artifact::ReportArtifact;kwargs...)
    get(io, :compact, false) && return show(io, artifact)
    artifact.tables isa NamedTuple && haskey(artifact.tables,:features) &&
        return _show_observed_benchmark(io,mime,artifact;kwargs...)
    dimensions = applicable(size, artifact.tables) ? join(size(artifact.tables), '×') : "one table"
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


function _show_observed_benchmark(io,mime,artifact;metric=:relative,problem=nothing,native_timings=false)
    metric in (:relative,:absolute) || throw(ArgumentError("metric must be :relative or :absolute"))
    features=filter(row -> problem===nothing || row.problem_index==problem,artifact.tables.features)
    isempty(features) && throw(ArgumentError("the selected point has no retained comparison features"))
    html=mime isa MIME"text/html"
    table_options=html ? (;) : (;truncate=0)
    escape(value)=replace(string(value),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;")
    heading(title)=html ? println(io,"<h3>",escape(title),"</h3>") : println(io,title)
    note(value)=html ? println(io,"<p>",escape(value),"</p>") : println(io,value)
    for feature in features
        label=metric===:relative ? "Relative RMS [%]" : "Absolute RMS [$(feature.absolute_unit)]"
        heading("$(feature.quantity) — $label · $(feature.statistic) · point $(feature.problem_index)")
        show(IOContext(io,:limit=>false),mime,getproperty(feature,metric);summary=false,eltypes=false,table_options...)
        println(io)
    end
    note("Maxima use eligible terms only; unavailable coefficients remain missing.")
    overview=artifact.tables.overview
    if isempty(overview.performance)
        note("No controlled performance measurements were recorded.")
    else
        heading("Controlled calculation measurements")
        note("Allocated bytes are cumulative Julia allocations, not peak memory.")
        any(==(1),overview.performance.samples) && note("At least one measurement has only one timed call; variability is unavailable.")
        any(!,overview.timing_ratio.comparable) && note("The recorded ratio is not a validated speedup.")
    end
    for (name,table) in pairs(overview)
        name===:source_timings && !native_timings && continue
        isempty(table) && continue
        heading(replace(string(name),'_'=>' '))
        shown=DataFrames.select(table,Not(intersect(propertynames(table),[:session_id,:checksum_verified,:workload_verified,:candidate_source,:reference_source])))
        show(IOContext(io,:limit=>false),mime,shown;summary=false,eltypes=false,table_options...)
        println(io)
    end
    return nothing
end
function Base.show(io::IO,mime::MIME"text/html",artifact::ReportArtifact;kwargs...)
    artifact.tables isa NamedTuple && haskey(artifact.tables,:features) &&
        return _show_observed_benchmark(io,mime,artifact;kwargs...)
    return show(io,artifact)
end
