@testmodule BaseParamsTestSupport begin
end

@testsnippet UseBaseParamsSupport begin
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.Materials
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
end

@testmodule DataModelTestSupport begin
end

@testsnippet UseDataModelSupport begin
    using DataFrames
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.Materials
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
    using LineCableModels.EarthProps
    using LineCableModels.Engine
    using LineCableModels.ImportExport
    import LineCableModels.DataModel: Insulator
end

@testmodule EngineTestSupport begin
end

@testsnippet UseEngineSupport begin
    using DataFrames
    using Measurements
    using Measurements: measurement, uncertainty, value
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.DataModel.BaseParams
    using LineCableModels.EarthProps
    using LineCableModels.Engine
    using LineCableModels.ParametricBuilder
    using LineCableModels.UQ
    using LineCableModels.ImportExport
    import LineCableModels.DataModel: Insulator
end

@testmodule PlotBuilderTestSupport begin
end

@testsnippet UsePlotBuilderSupport begin
    using DataFrames
    using Measurements
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.EarthProps
    using LineCableModels.Engine
    using LineCableModels.PlotBuilder
    using LineCableModels.UQ
    import LineCableModels.DataModel: Insulator
end

@testmodule ImportExportTestSupport begin
end

@testsnippet UseImportExportSupport begin
    using DataFrames
    using LineCableModels
    using LineCableModels.DataModel
    using LineCableModels.EarthProps
    using LineCableModels.Engine
    using LineCableModels.ImportExport
end
