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
