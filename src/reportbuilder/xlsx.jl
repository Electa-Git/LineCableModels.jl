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
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end

function XLSXReportDefinition(;
        file_name::Union{Nothing, AbstractString} = nothing,
        cable_system::Union{Nothing, DataModel.LineCableSystem} = nothing,
        clip::Bool = true
)
    path = file_name === nothing ? nothing : String(file_name)
    return XLSXReportDefinition(path, cable_system, clip)
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

entitle(::XLSXReportDefinition, source::Engine.LineParameters) = source
function select(definition::XLSXReportDefinition, source::Engine.LineParameters)
    line_definition = _line_definition(
        (
            @observe(R[:, :, :]),
            @observe(X[:, :, :]),
            @observe(G[:, :, :]),
            @observe(B[:, :, :])
        ),
        :base,
        :kilo,
        nothing,
        definition.clip
    )
    return select(line_definition, source)
end
function tabulate(::XLSXReportDefinition, source::Engine.LineParameters, selected)
    return _publication_table(selected)
end
illustrate(::XLSXReportDefinition, source, published, table) = nothing

function _family_columns(table::DataFrame, family::Val)
    contract = observation_columns(table)
    return Tuple(name
    for (name, entry) in pairs(contract)
    if applicable(Units.family, entry.quantity) &&
        Units.family(entry.quantity) === family)
end

function _is_diagonal(table::DataFrame, quantity_columns::Tuple)
    isempty(quantity_columns) && throw(ArgumentError(
        "workbook family has no observed quantity columns",
    ))
    selected = table
    frequency = first(selected.frequency)
    slice = selected[selected.frequency .== frequency, :]
    off_diagonal = slice.row .!= slice.column
    return all(quantity_columns) do quantity
        all(
            value -> isapprox(value, zero(value); rtol = 1.0e-8, atol = 1.0e-8),
            slice[off_diagonal, quantity]
        )
    end
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
        quantity_columns::Tuple,
        row::Int,
        column::Int
)
    selected = table[
        (table.row .== row) .& (table.column .== column),
        :
    ]
    isempty(selected) && throw(ArgumentError(
        "workbook sheet $name has no line-parameter rows",
    ))
    frequency_values = unique(selected.frequency)
    contract = observation_columns(table)
    unit_rows = Pair{String, String}[
        "frequency" => Units.label(contract.frequency.unit),
    ]
    for quantity in quantity_columns
        push!(unit_rows, String(quantity) => Units.label(contract[quantity].unit))
    end

    heading_row = length(unit_rows) + 2
    cells = fill(
        "",
        heading_row + length(frequency_values),
        max(2, length(quantity_columns) + 1)
    )
    for (unit_row, entry) in enumerate(unit_rows)
        cells[unit_row, 1] = first(entry)
        cells[unit_row, 2] = last(entry)
    end
    cells[heading_row, 1] = "frequency"
    for (quantity_column, quantity) in enumerate(quantity_columns)
        cells[heading_row, quantity_column + 1] = String(quantity)
    end
    for (frequency_index, frequency) in enumerate(frequency_values)
        data_row = heading_row + frequency_index
        cells[data_row, 1] = encode_cell(definition, frequency)
        for (quantity_column, quantity) in enumerate(quantity_columns)
            cells[data_row, quantity_column + 1] = encode_cell(
                definition,
                selected[frequency_index, quantity]
            )
        end
    end
    return XLSXSheet(name, cells)
end

function _encoded_sheets(
        definition::XLSXReportDefinition,
        prefix::String,
        table::DataFrame,
        quantity_columns::Tuple,
        diagonal::Bool
)
    rows = maximum(table.row)
    columns = maximum(table.column)
    indices = diagonal ? ((index, index) for index in 1:min(rows, columns)) :
              ((row, column) for row in 1:rows for column in 1:columns)
    return [_encoded_sheet(
                definition,
                table,
                "$prefix($row,$column)",
                quantity_columns,
                row,
                column
            )
            for (row, column) in indices]
end

function encode(
        definition::XLSXReportDefinition,
        ::Engine.LineParameters,
        published,
        table,
        ::Nothing
)
    isempty(table) && throw(ArgumentError(
        "XLSX reports require at least one frequency sample",
    ))
    series_columns = _family_columns(table, Val(:series))
    shunt_columns = _family_columns(table, Val(:shunt))
    series_diagonal = _is_diagonal(table, series_columns)
    shunt_diagonal = _is_diagonal(table, shunt_columns)
    series_diagonal &&
        @warn("Z is diagonal within the selected tolerance. Exporting Z[i,i] and omitting off-diagonal elements.")
    shunt_diagonal &&
        @warn("Y is diagonal within the selected tolerance. Exporting Y[i,i] and omitting off-diagonal elements.")
    sheets = vcat(
        _encoded_sheets(definition, "Z", table, series_columns, series_diagonal),
        _encoded_sheets(definition, "Y", table, shunt_columns, shunt_diagonal)
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
