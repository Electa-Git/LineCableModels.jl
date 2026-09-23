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

function Base.summary(io::IO, artifact::ReportArtifact)
    groups=_reported_tables(artifact)
    point_count=length(_observed_points(artifact.observed))
    table_count=sum(group -> length(group.tables),groups;init=0)
    print(io,"Report")
    point_count>1 && print(io," · ",point_count," gridpoints")
    print(io," · ",table_count,table_count==1 ? " table" : " tables")
end

function Base.show(io::IO, artifact::ReportArtifact)
    groups=_reported_tables(artifact)
    point_count=length(_observed_points(artifact.observed))
    print(io,"ReportArtifact(")
    point_count>1 && print(io,"gridpoints=",point_count,", ")
    print(io,"tables=",sum(group -> length(group.tables),groups;init=0))
    artifact.illustration===nothing || print(io,", illustration=present")
    artifact.output===nothing || print(io,", output=present")
    print(io,")")
end

_report_escape(value)=replace(string(value),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;",'\''=>"&#39;")

function _report_line(io,::MIME"text/plain",text;level=0,first=false)
    first || print(io,'\n')
    line=replace(string(text),'\n'=>' ','\r'=>' ')
    print(io,get(io,:limit,true) ? TextDisplay._fit(line,max(displaysize(io)[2],0)) : line)
end
function _report_line(io,::MIME"text/html",text;level=0,first=false)
    tag=level==0 ? "p" : "h$(level)"
    print(io,"<",tag,">",_report_escape(text),"</",tag,">")
end

function _reported_heading(table;families=false,limited=false)
    quantity=metadata(table,"quantity",nothing)
    unit=metadata(table,"unit",nothing)
    title=quantity===nothing ? "Table" : unit===nothing ? Units.label(quantity) : Units.label(quantity,unit)
    statistic=metadata(table,"statistic",:value)
    statistic===:value || (title*=" · "*string(statistic))
    coordinates=metadata(table,"coordinates",(;))
    get(coordinates,:kind,nothing)===:diagonal && (title*=" · diagonal")
    family=metadata(table,"family",nothing)
    families && family in (:Z,:Y) && (title=string(family," · ",title))
    columns=metadata(table,"observation_columns",(;))
    for name in metadata(table,"coordinate_columns",())
        descriptor=get(columns,name,nothing)
        text=descriptor===nothing || descriptor.quantity===nothing ? string(name) :
            Units.label(descriptor.quantity,descriptor.unit)
        title*=" · "*text
    end
    if limited
        rows,columns=size(table)
        title*=" · $rows row$(rows==1 ? "" : "s") × $columns column$(columns==1 ? "" : "s")"
    end
    return title
end

# Native DataFrames writers own cell formatting and whole-row/column omission.
# A limited table is buffered only within its allotted display area, to account
# for its actual line use. Unlimited inspection streams directly to the caller.
function _report_table(io,mime::MIME"text/plain",table,height,width,limited)
    if !limited
        print(io,'\n')
        show(IOContext(io,:limit=>false),mime,table;summary=false,eltypes=false,truncate=0)
        return 0
    end
    # Native horizontal cropping can cut a numeric token. Measure a bounded
    # native preview without cell clipping, dropping whole columns until it fits.
    columns=min(DataFrames.ncol(table),max(width÷3,1))
    while columns>0 || DataFrames.ncol(table)==0
        buffer=IOBuffer()
        context=IOContext(IOContext(buffer,io),:limit=>true,:displaysize=>(height,width))
        show(context,mime,table;summary=false,eltypes=false,truncate=0,
            reserved_display_lines=0,maximum_number_of_rows=max(height-3,1),
            maximum_number_of_columns=columns,allrows=true,allcols=true,
            show_omitted_cell_summary=false)
        rendered=String(take!(buffer))
        lines=split(rendered,'\n')
        # Ignore terminal color escapes for width measurement only.
        if all(line -> textwidth(replace(line,r"\e\[[0-9;:]*m"=>""))<=width,lines)
            print(io,'\n',rendered)
            return length(lines)
        end
        columns-=1
        DataFrames.ncol(table)==0 && break
    end
    _report_line(io,mime,"… no complete column fits")
    return 1
end
function _report_table(io,mime::MIME"text/html",table,height,width,limited)
    show(IOContext(io,:limit=>limited),mime,table;summary=false,eltypes=false,
        maximum_number_of_rows=limited ? max(height-3,1) : -1,
        maximum_number_of_columns=limited ? max(width÷12,1) : -1,
        show_omitted_cell_summary=!limited,new_line_at_end=false)
    return limited ? min(height,DataFrames.nrow(table)+2) : 0
end

function _show_quantity_report(io,mime,artifact)
    groups=_reported_tables(artifact)
    total=sum(group -> length(group.tables),groups;init=0)
    limited=get(io,:limit,true)
    rows,width=displaysize(io)
    point_count=length(_observed_points(artifact.observed))
    table_noun=total==1 ? "table" : "tables"
    context=point_count>1 ? "Report · $point_count gridpoints" : "Report"
    _report_line(io,mime,"$context · $total $table_noun";level=2,first=true)
    limited && rows<=1 && return nothing
    # One final line remains available for honest omission counts.
    remaining=limited ? max(rows-2,0) : typemax(Int)
    points=collect(_observed_points(artifact.observed))
    labels=Grammar.observation_labels(artifact.reference===nothing ? points : [points;artifact.reference];fallback="")
    if artifact.reference!==nothing && remaining>0
        _report_line(io,mime,isempty(last(labels)) ? "Reference" : "Reference · $(last(labels))")
        remaining-=1
    end
    shown=0
    shown_points=0
    for group in groups
        index=group.point_index
        tables=group.tables
        isempty(tables) && continue
        context=index===nothing ? join(unique(filter(!isempty,labels[1:point_count])),"; ") : labels[index]
        if index!==nothing && point_count>1 && (isempty(context) || count(==(context),labels[1:point_count])>1)
            context="[$index]"*(isempty(context) ? "" : " · $context")
        end
        context_lines=isempty(context) ? 0 : 1
        # Heading, column labels, separator, at least one data row, and native
        # omission information when needed. No gridpoint consumes the next one's
        # place unless its own selected quantities have received a fair preview.
        costs=[3+min(DataFrames.nrow(table),1)+(DataFrames.nrow(table)>1) for table in tables]
        visible=limited ? count(<=(remaining-context_lines),cumsum(costs)) : length(tables)
        (visible==0 || limited && width<12) && break
        isempty(context) || _report_line(io,mime,context;level=3)
        remaining-=context_lines
        shown_points+=index===nothing ? point_count : 1
        families=length(unique(metadata(table,"family",nothing) for table in tables))>1
        for position in 1:visible
            table=tables[position]
            remaining-=1
            height=limited ? max(costs[position]-1,remaining÷(visible-position+1)-1) : 0
            _report_line(io,mime,_reported_heading(table;families,limited);level=4)
            used=_report_table(io,mime,table,height,width,limited)
            remaining-=used
            shown+=1
        end
        visible<length(tables) && break
    end
    if shown<total
        progress=point_count>1 ? "; $shown_points/$point_count gridpoints shown" : ""
        _report_line(io,mime,"… $(total-shown) tables omitted$progress")
    elseif remaining>0 && artifact.illustration!==nothing
        _report_line(io,mime,"Illustration retained")
        remaining-=1
    end
    if remaining>0 && artifact.output!==nothing
        output=artifact.output isa AbstractString ? "Written output: $(artifact.output)" :
            artifact.output isa AbstractVector ? "$(length(artifact.output)) written destinations" : "Written output retained"
        _report_line(io,mime,output)
    end
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", artifact::ReportArtifact;kwargs...)
    get(io,:compact,false) && return show(io,artifact)
    artifact.tables isa NamedTuple && haskey(artifact.tables,:features) &&
        return _show_observed_benchmark(io,mime,artifact;kwargs...)
    return _show_quantity_report(io,mime,artifact)
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
    heading(title)=html ? println(io,"<h3>",_report_escape(title),"</h3>") : println(io,title)
    note(value)=html ? println(io,"<p>",_report_escape(value),"</p>") : println(io,value)
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
    get(io,:compact,false) && return show(io,artifact)
    artifact.tables isa NamedTuple && haskey(artifact.tables,:features) &&
        return _show_observed_benchmark(io,mime,artifact;kwargs...)
    return _show_quantity_report(io,mime,artifact)
end
