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

function _xlsx_wide_table(
        table::DataFrame,
        family::Symbol,
        row::Int,
        column::Int
)
    selected = table[
        (table.family .== family) .& (table.row .== row) .& (table.column .== column),
        :
    ]
    frequency_values = unique(selected.frequency)
    quantity_keys = unique(selected.quantity)
    wide = DataFrame(frequency = frequency_values)
    unit_labels = Dict{Symbol, String}(
        :frequency => metadata(table, "frequency_unit"),
    )
    headings = Dict{Symbol, String}(:frequency => "Frequency")
    for key in quantity_keys
        rows = selected[selected.quantity .== key, :]
        wide[!, key] = rows.value
        unit_labels[key] = only(unique(rows.unit))
        headings[key] = key |> String
    end
    metadata!(wide, "units", unit_labels, style = :note)
    metadata!(wide, "headings", headings, style = :note)
    return wide
end

function _xlsx_sheets(
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
    return [XLSXSheet(
                "$prefix($row,$column)",
                _xlsx_wide_table(table, family, row, column)
            )
            for (row, column) in indices]
end

function encode(
        ::XLSXReportDefinition,
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
        _xlsx_sheets("Z", table, :series, series_diagonal),
        _xlsx_sheets("Y", table, :shunt, shunt_diagonal)
    )
    return XLSXWorkbook(sheets)
end

"""
$(TYPEDSIGNATURES)

Encode one value for a cell in an XLSX report definition.
"""
encode_cell(::XLSXReportDefinition, value) = string(value)
encode_cell(::XLSXReportDefinition, ::Missing) = ""
encode_cell(::XLSXReportDefinition, value::Real) = @sprintf("%.12g", float(value))

function _xlsx_strings(definition::XLSXReportDefinition, table::DataFrame)
    return DataFrame(
        (name => encode_cell.(Ref(definition), table[!, name]) for name in names(table))...;
        copycols = false
    )
end

function _xlsx_units(table::DataFrame)
    "units" in metadatakeys(table) || return nothing
    return metadata(table, "units")
end

function _xlsx_destination(definition::XLSXReportDefinition)
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
