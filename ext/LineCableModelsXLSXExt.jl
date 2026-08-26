"""
    LineCableModelsXLSXExt

Write ReportBuilder workbook descriptions with XLSX.jl.
"""
module LineCableModelsXLSXExt

using LineCableModels
using Tables
using XLSX

const ReportBuilder = LineCableModels.ReportBuilder
const Engine = LineCableModels.Engine

function _write_xlsx_sheet!(
        workbook,
        sheet::ReportBuilder.XLSXSheet,
        definition::ReportBuilder.XLSXReportDefinition;
        first_sheet::Bool
)
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
    units = ReportBuilder._xlsx_units(sheet.table)
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
        Tables.columntable(ReportBuilder._xlsx_strings(definition, sheet.table));
        anchor_cell = XLSX.CellRef(row, 1)
    )
    return nothing
end

function ReportBuilder.write(
        definition::ReportBuilder.XLSXReportDefinition,
        source::Engine.LineParameters,
        published,
        table,
        ::Nothing,
        encoded::ReportBuilder.XLSXWorkbook
)
    destination = ReportBuilder._xlsx_destination(definition)
    XLSX.openxlsx(destination, mode = "w") do workbook
        for (index, sheet) in enumerate(encoded.sheets)
            _write_xlsx_sheet!(
                workbook,
                sheet,
                definition;
                first_sheet = index == 1
            )
        end
    end
    return destination
end

end # module LineCableModelsXLSXExt
