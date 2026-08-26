"""
$(TYPEDEF)

Write the human-facing line-parameter workbook produced by [`report`](@ref).

`file_name` follows the established XLSX export path rules. When
`cable_system` is supplied, its system identifier prefixes the output name.

$(TYPEDFIELDS)
"""
struct XLSXReportDefinition{
    F <: Union{Nothing, String}, C <: Union{Nothing, DataModel.LineCableSystem}} <:
       AbstractReportDefinition
    "Requested workbook path, or `nothing` for the default path."
    file_name::F
    "Optional cable system used to prefix the workbook name."
    cable_system::C
end

function XLSXReportDefinition(;
        file_name::Union{Nothing, AbstractString} = nothing,
        cable_system::Union{Nothing, DataModel.LineCableSystem} = nothing
)
    path = file_name === nothing ? nothing : String(file_name)
    return XLSXReportDefinition(path, cable_system)
end

"""
$(TYPEDEF)

Hold the name and fully encoded cell grid for one workbook sheet.

$(TYPEDFIELDS)
"""
struct XLSXSheet
    "Workbook sheet name."
    name::String
    "Encoded cells indexed by worksheet row and column."
    cells::Matrix{String}
end

"""
$(TYPEDEF)

Describe one XLSX workbook without depending on XLSX.jl.

$(TYPEDFIELDS)
"""
struct XLSXWorkbook
    "Absolute caller-owned output path."
    destination::String
    "Sheets in workbook order."
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

entitle(::XLSXReportDefinition, source::Engine.LineParameters) = source
function select(::XLSXReportDefinition, source::Engine.LineParameters)
    select(_xlsx_table_definition(), source)
end
function tabulate(::XLSXReportDefinition, source::Engine.LineParameters, selected)
    tabulate(_xlsx_table_definition(), source, selected)
end
illustrate(::XLSXReportDefinition, source, published, table) = nothing

function _is_diagonal(matrix)
    return isapprox(matrix, Diagonal(diag(matrix)); rtol = 1.0e-8, atol = 1.0e-8)
end

"""
$(TYPEDSIGNATURES)

Encode one value for a cell in an [`XLSXWorkbook`](@ref).

External value owners may add narrow methods for their scalar types.
"""
encode_cell(::XLSXReportDefinition, value) = string(value)
encode_cell(::XLSXReportDefinition, ::Missing) = ""
encode_cell(::XLSXReportDefinition, value::Real) = @sprintf("%.12g", float(value))

function _encoded_sheet(
        definition::XLSXReportDefinition,
        table::DataFrame,
        name::String,
        family::Symbol,
        row::Int,
        column::Int
)
    selected = table[
        (table.family .== family) .& (table.row .== row) .& (table.column .== column),
        :
    ]
    isempty(selected) && throw(ArgumentError(
        "workbook sheet $name has no line-parameter rows",
    ))
    frequency_values = unique(selected.frequency)
    quantity_keys = unique(selected.quantity)
    quantity_rows = [selected[selected.quantity .== quantity, :]
                     for quantity in quantity_keys]
    unit_rows = Pair{String, String}[
        "frequency" => metadata(table, "frequency_unit"),
    ]
    for (quantity, rows) in zip(quantity_keys, quantity_rows)
        rows.frequency == frequency_values || throw(DimensionMismatch(
            "workbook quantity $quantity does not align with the frequency column",
        ))
        push!(unit_rows, String(quantity) => only(unique(rows.unit)))
    end

    heading_row = length(unit_rows) + 2
    cells = fill(
        "",
        heading_row + length(frequency_values),
        max(2, length(quantity_keys) + 1)
    )
    for (unit_row, entry) in enumerate(unit_rows)
        cells[unit_row, 1] = first(entry)
        cells[unit_row, 2] = last(entry)
    end
    cells[heading_row, 1] = "frequency"
    for (quantity_column, quantity) in enumerate(quantity_keys)
        cells[heading_row, quantity_column + 1] = String(quantity)
    end
    for (frequency_index, frequency) in enumerate(frequency_values)
        data_row = heading_row + frequency_index
        cells[data_row, 1] = encode_cell(definition, frequency)
        for (quantity_column, rows) in enumerate(quantity_rows)
            cells[data_row, quantity_column + 1] = encode_cell(
                definition,
                rows.value[frequency_index]
            )
        end
    end
    return XLSXSheet(name, cells)
end

function _encoded_sheets(
        definition::XLSXReportDefinition,
        prefix::String,
        table::DataFrame,
        family::Symbol,
        diagonal::Bool
)
    selected = table[table.family .== family, :]
    rows = maximum(selected.row)
    columns = maximum(selected.column)
    indices = diagonal ? ((index, index) for index in 1:min(rows, columns)) :
              ((row, column) for row in 1:rows for column in 1:columns)
    return [_encoded_sheet(
                definition,
                table,
                "$prefix($row,$column)",
                family,
                row,
                column
            )
            for (row, column) in indices]
end

function encode(
        definition::XLSXReportDefinition,
        source::Engine.LineParameters,
        published,
        table,
        ::Nothing
)
    isempty(frequencies(source)) && throw(ArgumentError(
        "XLSX reports require at least one frequency sample",
    ))
    series_diagonal = _is_diagonal(Z(source, :, :, 1))
    shunt_diagonal = _is_diagonal(Y(source, :, :, 1))
    series_diagonal &&
        @warn("Z is diagonal within the selected tolerance. Exporting Z[i,i] and omitting off-diagonal elements.")
    shunt_diagonal &&
        @warn("Y is diagonal within the selected tolerance. Exporting Y[i,i] and omitting off-diagonal elements.")
    sheets = vcat(
        _encoded_sheets(definition, "Z", table, :series, series_diagonal),
        _encoded_sheets(definition, "Y", table, :shunt, shunt_diagonal)
    )
    requested = abspath(something(definition.file_name, "ZY_export.xlsx"))
    destination = if definition.cable_system === nothing
        requested
    else
        joinpath(
            dirname(requested),
            "$(definition.cable_system.system_id)_$(basename(requested))"
        )
    end
    return XLSXWorkbook(String(destination), sheets)
end
