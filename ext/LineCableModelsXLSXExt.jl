"""
    LineCableModelsXLSXExt

Write ReportBuilder workbook descriptions with XLSX.jl.
"""
module LineCableModelsXLSXExt

import LineCableModels.ReportBuilder
import LineCableModels.Engine
import XLSX

function ReportBuilder.write(
        ::ReportBuilder.XLSXReportDefinition,
        encoded::ReportBuilder.XLSXWorkbook
)
    XLSX.openxlsx(encoded.destination, mode = "w") do workbook
        for (index, sheet) in enumerate(encoded.sheets)
            worksheet = if index == 1
                existing = workbook["Sheet1"]
                XLSX.renamesheet!(existing, sheet.name)
                existing
            else
                XLSX.addsheet!(workbook, sheet.name)
            end
            for cell in CartesianIndices(sheet.cells)
                value = sheet.cells[cell]
                ismissing(value) || (worksheet[cell[1], cell[2]] = value)
            end
        end
    end
    return encoded.destination
end

function ReportBuilder.write(definition::ReportBuilder.XLSXReportDefinition, books::AbstractVector{<:ReportBuilder.XLSXWorkbook})
    return [ReportBuilder.write(definition,book) for book in books]
end

end # module LineCableModelsXLSXExt
