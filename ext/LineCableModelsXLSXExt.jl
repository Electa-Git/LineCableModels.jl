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
        ::Engine.LineParameters,
        ::Any,
        ::Any,
        ::Nothing,
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
                isempty(value) || (worksheet[cell[1], cell[2]] = value)
            end
        end
    end
    return encoded.destination
end

end # module LineCableModelsXLSXExt
