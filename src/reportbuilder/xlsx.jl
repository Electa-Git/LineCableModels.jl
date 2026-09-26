"""
$(TYPEDEF)

Write one workbook per gridpoint and retained quantity. `file_name` supplies a
path prefix. Numeric values and standard deviations occupy separate sheets.

$(TYPEDFIELDS)
"""
struct XLSXReportDefinition <: AbstractReportDefinition
    "Requested output path prefix, or nothing for observed.xlsx."
    file_name::Union{Nothing,String}
    "Optional system name used only when constructing output filenames."
    system_id::Union{Nothing,String}
    "Engineering recentering for raw-input construction."
    clip::Bool
    "Allow replacing existing destination files after complete preflight."
    overwrite::Bool
end
function XLSXReportDefinition(;file_name=nothing,cable_system=nothing,clip::Bool=true,overwrite::Bool=false)
    return XLSXReportDefinition(file_name===nothing ? nothing : String(file_name),
        cable_system===nothing ? nothing : String(cable_system.system_id),clip,overwrite)
end

"""
$(TYPEDEF)

An encoded sheet with numeric cells and textual headings.

$(TYPEDFIELDS)
"""
struct XLSXSheet
    "Worksheet name."
    name::String
    "Numeric cells, textual metadata, or missing cells."
    cells::Matrix{Any}
end
"""
$(TYPEDEF)

A workbook ready for the XLSX writer.

$(TYPEDFIELDS)
"""
struct XLSXWorkbook
    "Absolute output path."
    destination::String
    "Sheets in output order."
    sheets::Vector{XLSXSheet}
end

"""
$(TYPEDSIGNATURES)

Convert a spreadsheet number to Float64. Reject nonfinite values, overflow, and
nonzero underflow. Native persistence retains original precision and uncertainty
dependencies; spreadsheets contain nominal values and standard deviations.
"""
function encode_cell(::XLSXReportDefinition,value::Real)
    number=Float64(value)
    isfinite(number) && (iszero(number) ? iszero(value) : true) ||
        throw(ArgumentError("value cannot be represented as a finite non-underflowing XLSX number"))
    return number
end
encode_cell(::XLSXReportDefinition,value::Complex) = throw(ArgumentError("XLSX quantities require retained real components; select real-valued primary or statistical products"))
encode_cell(::XLSXReportDefinition,::Missing) = missing
encode_cell(::XLSXReportDefinition,value::AbstractString) = String(value)
encode_cell(::XLSXReportDefinition,value) = string(value)

function _numeric_sheet(definition,table,name,transform)
    _xlsx_sheet_size(size(table,1)+1,size(table,2))
    cells=Matrix{Any}(missing,size(table,1)+1,size(table,2))
    cells[1,:]=names(table)
    coordinate_columns=DataFrames.metadata(table,"coordinate_columns",())
    for j in 1:size(table,2),i in 1:size(table,1)
        value=table[i,j]
        cells[i+1,j]=ismissing(value) ? missing : encode_cell(definition,propertynames(table)[j] in coordinate_columns ? Grammar.nominal(value) : transform(value))
    end
    return XLSXSheet(name,cells)
end

select(::XLSXReportDefinition,observed::ObservedResult;reference=nothing) = observed.quantities
function tabulate(::XLSXReportDefinition,observed,selected;reference=nothing)
    observed isa ObservedResult && return _quantity_tables(selected;gridpoint_id=observed.gridpoint.id)
    return map((point,products) -> _quantity_tables(products;gridpoint_id=point.gridpoint.id),
        observed,selected)
end
function report(definition::XLSXReportDefinition,source::Engine.LineParameters;kwargs...)
    return report(definition,ObservedResult(source;clip=definition.clip,kwargs...))
end
function encode(definition::XLSXReportDefinition,observed,tables,illustration;reference=nothing)
    workbooks=XLSXWorkbook[]
    requested=abspath(something(definition.file_name,"observed.xlsx"))
    stem=splitext(basename(requested))[1]
    definition.system_id===nothing || (stem=definition.system_id*"_"*stem)
    for (point_index,point) in enumerate(_observed_points(observed)), product in point.quantities
        point_tables=observed isa ObservedResult ? tables : tables[point_index]
        table=getproperty(getproperty(point_tables,product.family),_quantity_name(product))
        name=replace(string(_quantity_name(product)),r"[^A-Za-z0-9_-]"=>"_")
        id=get(point.gridpoint,:id,nothing)
        index=id===nothing ? string(point_index) : string(id.source_id,"_",id.problem_index,"_",id.formulation_index)
        destination=joinpath(dirname(requested),stem*"_"*index*"_"*name*".xlsx")
        any(book -> book.destination==destination,workbooks) && throw(ArgumentError("duplicate workbook destination"))
        entries=(quantity=string(product.quantity),unit=Units.label(product.unit),basis=string(product.basis),
            coordinates=repr(product.coordinates),gridpoint=repr(point.gridpoint),
            cutoffs=repr(product.thresholds),missing_reason=repr(product.missing_reason))
        cells=Matrix{Any}(undef,length(entries),2)
        for (i,(key,value)) in enumerate(pairs(entries))
            cells[i,1]=string(key); cells[i,2]=value
        end
        sheets=[_numeric_sheet(definition,table,"values",Grammar.nominal),
            _numeric_sheet(definition,table,"std",Grammar.uncertainty),XLSXSheet("metadata",cells)]
        push!(workbooks,XLSXWorkbook(destination,sheets))
    end
    validate(definition,workbooks)
    return workbooks
end

function _xlsx_sheet_size(rows,columns)
    rows<=1_048_576 && columns<=16_384 || throw(DimensionMismatch(
        "XLSX worksheets support at most 1048576 rows and 16384 columns, including headings"))
    return nothing
end

function validate(definition::XLSXReportDefinition,workbooks::AbstractVector{<:XLSXWorkbook})
    destinations=[abspath(book.destination) for book in workbooks]
    allunique(destinations) || throw(ArgumentError("duplicate workbook destination"))
    for book in workbooks
        ispath(book.destination) && (!definition.overwrite || !isfile(book.destination)) &&
            throw(ArgumentError("XLSX destination already exists: $(book.destination); use overwrite=true to replace files"))
        isdir(dirname(book.destination)) || throw(ArgumentError("XLSX destination directory does not exist"))
        for sheet in book.sheets
            _xlsx_sheet_size(size(sheet.cells)...)
        end
    end
    return nothing
end
