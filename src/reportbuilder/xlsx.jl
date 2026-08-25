"""
$(TYPEDEF)

Write the human-facing line-parameter workbook produced by [`report`](@ref).

`file_name` follows the established XLSX export path rules. When
`cable_system` is supplied, its system identifier prefixes the output name.

$(TYPEDFIELDS)
"""
struct XLSXReport{F <: Union{Nothing, String}, C <: Union{Nothing, DataModel.LineCableSystem}} <:
       AbstractReportDefinition
    "Requested workbook path, or `nothing` for the default path."
    file_name::F
    "Optional cable system used to prefix the workbook name."
    cable_system::C
end

function XLSXReport(;
        file_name::Union{Nothing, AbstractString} = nothing,
        cable_system::Union{Nothing, DataModel.LineCableSystem} = nothing
)
    path = file_name === nothing ? nothing : String(file_name)
    return XLSXReport(path, cable_system)
end

struct XLSXSheet
    name::String
    table::DataFrame
end

struct XLSXWorkbook
    sheets::Vector{XLSXSheet}
end

function _xlsx_table_definition()
    return _line_definition(
        (),
        nothing,
        :base,
        :kilo,
        nothing,
        sqrt(eps(Float64))
    )
end

entitle(::XLSXReport, source::Engine.LineParameters) = source
select(::XLSXReport, source::Engine.LineParameters) =
    select(_xlsx_table_definition(), source)
tabulate(::XLSXReport, source::Engine.LineParameters, selected) =
    tabulate(_xlsx_table_definition(), source, selected)
illustrate(::XLSXReport, source, published, table) = nothing

function _is_diagonal(matrix)
    return isapprox(matrix, Diagonal(diag(matrix)); rtol = 1.0e-8, atol = 1.0e-8)
end

function _xlsx_sheets(prefix::String, tables::Matrix{DataFrame}, diagonal::Bool)
    rows, columns = size(tables)
    indices = diagonal ? ((index, index) for index in 1:min(rows, columns)) :
              ((row, column) for row in 1:rows for column in 1:columns)
    return [XLSXSheet("$prefix($row,$column)", tables[row, column])
            for (row, column) in indices]
end

function encode(
        ::XLSXReport,
        source::Engine.LineParameters,
        published,
        table,
        ::Nothing
)
    isempty(frequencies(source)) && throw(ArgumentError(
        "XLSX reports require at least one frequency sample",
    ))
    series_tables, shunt_tables = table
    series_diagonal = _is_diagonal(Engine.Z(source, :, :, 1))
    shunt_diagonal = _is_diagonal(Engine.Y(source, :, :, 1))
    series_diagonal && @warn(
        "Z is diagonal within the selected tolerance. Exporting Z[i,i] and omitting off-diagonal elements."
    )
    shunt_diagonal && @warn(
        "Y is diagonal within the selected tolerance. Exporting Y[i,i] and omitting off-diagonal elements."
    )
    sheets = vcat(
        _xlsx_sheets("Z", series_tables, series_diagonal),
        _xlsx_sheets("Y", shunt_tables, shunt_diagonal)
    )
    return XLSXWorkbook(sheets)
end

_xlsx_string(value) = string(value)
_xlsx_string(::Missing) = ""
_xlsx_string(value::Real) = @sprintf("%.12g", float(value))

function _xlsx_strings(table::DataFrame)
    return DataFrame(
        (name => _xlsx_string.(table[!, name]) for name in names(table))...;
        copycols = false
    )
end

function _xlsx_units(table::DataFrame)
    "units" in metadatakeys(table) || return nothing
    return metadata(table, "units")
end

function _write_xlsx_sheet!(workbook, sheet::XLSXSheet; first_sheet::Bool)
    worksheet = if first_sheet
        existing = try
            workbook["Sheet1"]
        catch
            nothing
        end
        selected = existing === nothing ? XLSX.addsheet!(workbook, sheet.name) : existing
        try
            XLSX.rename!(selected, sheet.name)
        catch
        end
        selected
    else
        XLSX.addsheet!(workbook, sheet.name)
    end

    row = 1
    units = _xlsx_units(sheet.table)
    if units isa AbstractDict
        for name in names(sheet.table)
            worksheet[row, 1] = String(name)
            unit = get(units, name, get(units, Symbol(name), ""))
            worksheet[row, 2] = String(unit)
            row += 1
        end
        row += 1
    end
    XLSX.writetable!(
        worksheet,
        Tables.columntable(_xlsx_strings(sheet.table));
        anchor_cell = XLSX.CellRef(row, 1)
    )
    return nothing
end

function _xlsx_destination(definition::XLSXReport)
    base_directory = joinpath(dirname(@__DIR__), "importexport")
    file_name = if definition.file_name === nothing
        identifier = definition.cable_system === nothing ? "" :
                     "$(definition.cable_system.system_id)_"
        joinpath(base_directory, "$(identifier)ZY_export.xlsx")
    else
        requested = isabspath(definition.file_name) ? definition.file_name :
                    joinpath(base_directory, definition.file_name)
        if definition.cable_system === nothing
            requested
        else
            joinpath(
                dirname(requested),
                "$(definition.cable_system.system_id)_$(basename(requested))"
            )
        end
    end
    return String(file_name)
end

function write(
        definition::XLSXReport,
        source::Engine.LineParameters,
        published,
        table,
        ::Nothing,
        encoded::XLSXWorkbook
)
    destination = _xlsx_destination(definition)
    XLSX.openxlsx(destination, mode = "w") do workbook
        for (index, sheet) in enumerate(encoded.sheets)
            _write_xlsx_sheet!(workbook, sheet; first_sheet = index == 1)
        end
    end
    return destination
end
