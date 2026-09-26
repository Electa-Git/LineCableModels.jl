@testsnippet UseBaseParamsSupport begin
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.Materials
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
end

@testsnippet UseDataModelSupport begin
    using DataFrames
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.Materials
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
    using LineCableModels.Earth
    using LineCableModels.Engine
    using LineCableModels.ImportExport
end

@testsnippet UseEngineSupport begin
    using DataFrames
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
    using LineCableModels.Earth
    using LineCableModels.Engine
    using LineCableModels.ParametricBuilder
    using LineCableModels.UQ
    using LineCableModels.ImportExport
end

@testsnippet UseNativePlotSupport begin
    using DataFrames
    using Measurements
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.Earth
    using LineCableModels.Engine
    using LineCableModels.UQ
end

@testsnippet UseImportExportSupport begin
    using DataFrames
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.Earth
    using LineCableModels.Engine
    using LineCableModels.ImportExport
end
