@testmodule BaseParamsTestSupport begin
    using Reexport
    @reexport using Measurements
    @reexport using Measurements: measurement, uncertainty, value
    @reexport using LineCableModels.Commons
    @reexport using LineCableModels.Materials
    @reexport using LineCableModels.DataModel
    @reexport using LineCableModels.DataModel.BaseParams
end

@testsnippet UseBaseParamsSupport begin
    using .BaseParamsTestSupport
end

@testmodule DataModelTestSupport begin
    using Reexport
    @reexport using DataFrames
    @reexport using Measurements
    @reexport using Measurements: measurement, uncertainty, value
    @reexport using LineCableModels.Commons
    @reexport using LineCableModels.Utils
    @reexport using LineCableModels.Materials
    @reexport using LineCableModels.DataModel
    @reexport using LineCableModels.DataModel.BaseParams
    @reexport using LineCableModels.EarthProps
    @reexport using LineCableModels.Engine
    @reexport using LineCableModels.ImportExport
end

@testsnippet UseDataModelSupport begin
    using .DataModelTestSupport
    import LineCableModels.DataModel: Insulator
end

@testmodule EngineTestSupport begin
    using Reexport
    @reexport using DataFrames
    @reexport using Measurements
    @reexport using Measurements: measurement, uncertainty, value
    @reexport using LineCableModels.Commons
    @reexport using LineCableModels.DataModel
    @reexport using LineCableModels.DataModel.BaseParams
    @reexport using LineCableModels.EarthProps
    @reexport using LineCableModels.Engine
    @reexport using LineCableModels.ParametricBuilder
    @reexport using LineCableModels.Computation
    @reexport using LineCableModels.ImportExport
end

@testsnippet UseEngineSupport begin
    using .EngineTestSupport
    import LineCableModels.DataModel: Insulator
end

@testmodule PlotBuilderTestSupport begin
    using Reexport
    @reexport using DataFrames
    @reexport using Measurements
    @reexport using LineCableModels.Commons
    @reexport using LineCableModels.DataModel
    @reexport using LineCableModels.EarthProps
    @reexport using LineCableModels.Engine
    @reexport using LineCableModels.PlotBuilder
    @reexport using LineCableModels.Computation
end

@testsnippet UsePlotBuilderSupport begin
    using .PlotBuilderTestSupport
    import LineCableModels.DataModel: Insulator
end

@testmodule ImportExportTestSupport begin
    using Reexport
    @reexport using DataFrames
    @reexport using LineCableModels.Commons
    @reexport using LineCableModels.DataModel
    @reexport using LineCableModels.EarthProps
    @reexport using LineCableModels.Engine
    @reexport using LineCableModels.ImportExport
end

@testsnippet UseImportExportSupport begin
    using .ImportExportTestSupport
end
