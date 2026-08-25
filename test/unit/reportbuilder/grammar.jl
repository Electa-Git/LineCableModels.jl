@testitem "ReportBuilder / grammar / publication and stage order" tags=[:unit] begin
    using DataFrames

    const RB = LineCableModels.ReportBuilder
    const PB = LineCableModels.PlotBuilder
    const U = LineCableModels.Units

    struct ReportProfile
        values::Vector{Float64}
    end
    profile_response(source::ReportProfile) = source.values

    LineCableModels.basis(::ReportProfile) = :total
    LineCableModels.observe(
        source::ReportProfile,
        ::typeof(profile_response),
        indices...
    ) = isempty(indices) ? source.values : getindex(source.values, indices...)
    LineCableModels.observables(::Type{<:ReportProfile}) = (profile_response,)
    U.quantity(::typeof(profile_response)) = U.Quantity{:report_response}()
    U.native_unit(::U.Quantity{:report_response}) = U.units(:base, :ohm)
    U.display_unit(::U.Quantity{:report_response}) = U.units(:milli, :ohm)
    U.label(::U.Quantity{:report_response}) = "Response"
    U.symbol(::U.Quantity{:report_response}) = "u"

    struct ReportProfilePlot <: PB.AbstractPlotDefinition end
    PB.dispatch_on(::Type{ReportProfilePlot}) = ReportProfile
    function PB.fetch(::Type{ReportProfilePlot}, recipe::PB.PlotRecipe)
        published = LineCableModels.observables(
            recipe.object,
            (response = profile_response,)
        )
        return PB.PlotRecipe(
            ReportProfilePlot,
            recipe.object,
            merge(recipe.input, (; published)),
            recipe.renderer
        )
    end
    PB.axis_payload(
        ::Type{ReportProfilePlot},
        ::Val,
        ::Val{:y},
        recipe::PB.PlotRecipe,
        page_key,
        view_key
    ) = recipe.input.published.response
    PB.series_values(
        ::Type{ReportProfilePlot},
        ::Val,
        ::Val{:x},
        recipe::PB.PlotRecipe,
        page_key,
        view_key,
        series_key
    ) = eachindex(recipe.input.published.response.values)
    PB.series_values(
        ::Type{ReportProfilePlot},
        ::Val,
        ::Val{:y},
        recipe::PB.PlotRecipe,
        page_key,
        view_key,
        series_key
    ) = recipe.input.published.response.values
    PB.default_title(
        ::Type{ReportProfilePlot},
        ::Val,
        recipe::PB.PlotRecipe,
        page_key,
        view_key
    ) = "Report profile"

    source = ReportProfile([1.0, 2.0])
    definition = TableReport(
        (response = profile_response,);
        illustration = ReportProfilePlot
    )
    artifact = report(definition, source)

    @test artifact isa ReportArtifact
    @test parentmodule(typeof(artifact)) === RB
    @test artifact.table.response == [1_000.0, 2_000.0]
    @test DataFrames.metadata(artifact.table, "units")[:response] == "mΩ"
    @test DataFrames.metadata(artifact.table, "headings")[:response] ==
          "Response [mΩ]"
    @test artifact.illustration isa PB.PlotRecipe
    @test artifact.output === nothing
    artifact.table.response[1] = 0.0
    @test source.values == [1.0, 2.0]
    @test_throws ArgumentError report(TableReport((value = identity,)), source)
    @test_throws ArgumentError report(TableReport((value = identity,)), :unsupported)

    struct StageReport <: RB.AbstractReportDefinition end
    const report_stage_calls = Symbol[]
    function RB.entitle(::StageReport, source)
        push!(report_stage_calls, :entitle)
        return source
    end
    function RB.select(::StageReport, source)
        push!(report_stage_calls, :select)
        return :published
    end
    function RB.tabulate(::StageReport, source, published)
        push!(report_stage_calls, :tabulate)
        return :table
    end
    function RB.illustrate(::StageReport, source, published, table)
        push!(report_stage_calls, :illustrate)
        return :illustration
    end
    function RB.encode(::StageReport, source, published, table, illustration)
        push!(report_stage_calls, :encode)
        return :encoded
    end
    function RB.write(
            ::StageReport,
            source,
            published,
            table,
            illustration,
            encoded
    )
        push!(report_stage_calls, :write)
        return :written
    end
    function RB.finish(
            ::StageReport,
            source,
            published,
            table,
            illustration,
            encoded,
            written
    )
        push!(report_stage_calls, :finish)
        return (; table, illustration, encoded, written)
    end

    completed = report(StageReport(), :source)
    @test completed == (
        table = :table,
        illustration = :illustration,
        encoded = :encoded,
        written = :written
    )
    @test report_stage_calls == [
        :entitle,
        :select,
        :tabulate,
        :illustrate,
        :encode,
        :write,
        :finish
    ]

    struct MissingEncodeReport <: RB.AbstractReportDefinition end
    RB.entitle(::MissingEncodeReport, source) = source
    RB.select(::MissingEncodeReport, source) = :published
    RB.tabulate(::MissingEncodeReport, source, published) = :table
    RB.illustrate(::MissingEncodeReport, source, published, table) = nothing
    @test_throws MethodError report(MissingEncodeReport(), :source)

    struct MissingWriteReport <: RB.AbstractReportDefinition end
    RB.entitle(::MissingWriteReport, source) = source
    RB.select(::MissingWriteReport, source) = :published
    RB.tabulate(::MissingWriteReport, source, published) = :table
    RB.illustrate(::MissingWriteReport, source, published, table) = nothing
    RB.encode(::MissingWriteReport, source, published, table, illustration) = :encoded
    @test_throws MethodError report(MissingWriteReport(), :source)
end

@testitem "ReportBuilder / XLSX / workbook pipeline and delegation" tags=[:integration] setup=[
    EngineTestSupport,
    UseEngineSupport,
    TestFixtures
] begin
    using LinearAlgebra
    using Measurements
    using XLSX

    const RB=LineCableModels.ReportBuilder
    const IE=LineCableModels.ImportExport
    xlsx_extension=Base.get_extension(LineCableModels, :LineCableModelsXLSXExt)
    @test xlsx_extension !== nothing
    @test any(method -> method.module === xlsx_extension, methods(RB.write))

    frequency=[50.0, 500.0]
    impedance=Array{ComplexF64}(undef, 2, 2, 2)
    admittance=similar(impedance)
    for index in eachindex(frequency)
        impedance[:, :, index]=[1.0+2.0im 0.2+0.3im; 0.2+0.3im 1.5+2.5im]
        admittance[:, :, index]=[3.0+4.0im 0.4+0.5im; 0.4+0.5im 3.5+4.5im] .* 1.0e-6
    end
    parameters=LineParameters(impedance, admittance, frequency)

    mktempdir() do directory
        path=joinpath(directory, "full.xlsx")
        artifact=report(XLSXReport(; file_name = path), parameters)
        @test artifact isa ReportArtifact
        @test fieldnames(typeof(artifact)) == (:table, :illustration, :output)
        @test artifact.illustration === nothing
        @test artifact.output == path
        @test isfile(artifact.output)
        @test artifact.table isa Tuple
        @test length(artifact.table) == 2
        XLSX.openxlsx(path) do workbook
            @test Set(XLSX.sheetnames(workbook)) == Set([
                "Z(1,1)", "Z(1,2)", "Z(2,1)", "Z(2,2)",
                "Y(1,1)", "Y(1,2)", "Y(2,1)", "Y(2,2)"
            ])
            worksheet=workbook["Z(1,1)"]
            @test worksheet["A1"] == "frequency"
            @test worksheet["B1"] == "Hz"
            @test worksheet["A2"] == "real"
            @test worksheet["B2"] == "Ω/km"
            @test worksheet["A5"] == "frequency"
            @test worksheet["B5"] == "real"
            @test worksheet["C5"] == "imag"
            @test worksheet["A6"] == "50"
            @test worksheet["B6"] == "1000"
            @test worksheet["C6"] == "2000"
        end

        delegated=export_data(
            :xlsx,
            parameters;
            file_name = joinpath(directory, "delegated.xlsx")
        )
        @test delegated == joinpath(directory, "delegated.xlsx")
        @test isfile(delegated)
        @test parentmodule(which(
            IE.export_data,
            (Val{:xlsx}, typeof(parameters))
        )) === IE

        cable_system=TestFixtures.three_phase_system()
        prefixed=report(
            XLSXReport(
                file_name = joinpath(directory, "named.xlsx"),
                cable_system = cable_system
            ),
            parameters
        )
        @test basename(prefixed.output) == "$(cable_system.system_id)_named.xlsx"

        diagonal=LineParameters(
            cat(Diagonal([1.0+2.0im, 2.0+3.0im]); dims = 3),
            cat(Diagonal([3.0+4.0im, 4.0+5.0im]); dims = 3),
            [50.0]
        )
        diagonal_artifact=@test_logs (:warn, r"Z is diagonal") (:warn, r"Y is diagonal") report(
            XLSXReport(file_name = joinpath(directory, "diagonal.xlsx")),
            diagonal
        )
        XLSX.openxlsx(diagonal_artifact.output) do workbook
            @test Set(XLSX.sheetnames(workbook)) ==
                  Set(["Z(1,1)", "Z(2,2)", "Y(1,1)", "Y(2,2)"])
        end

        uncertain=LineParameters(
            reshape([complex(measurement(1.0e-4, 1.0e-5), 2.0e-4)], 1, 1, 1),
            reshape([complex(measurement(3.0e-8, 1.0e-9), 4.0e-8)], 1, 1, 1),
            [50.0]
        )
        uncertain_artifact=@test_logs (:warn, r"Z is diagonal") (:warn, r"Y is diagonal") report(
            XLSXReport(file_name = joinpath(directory, "uncertain.xlsx")),
            uncertain
        )
        XLSX.openxlsx(uncertain_artifact.output) do workbook
            @test occursin("±", workbook["Z(1,1)"]["B6"])
        end

        @test_throws Exception report(XLSXReport(file_name = directory), parameters)
    end

    @test !isdefined(IE, :XLSXWorkbook)
    @test !isdefined(IE, :_write_xlsx_sheet!)
    @test !isdefined(IE, :df_to_strings)
    @test !isdefined(RB, :_write_xlsx_sheet!)
end

@testitem "ReportBuilder / adapters / completed results delegate" tags=[:unit] setup=[
    TestFixtures
] begin
    using DataFrames

    const RB=LineCableModels.ReportBuilder
    constants=LineCableModels.CableConstants(1.0, 2.0, 3.0)
    expected=report(RB.CableConstantsTable(), constants).table
    actual=DataFrame(constants)

    @test actual == expected
    @test parentmodule(which(DataFrame, (typeof(constants),))) === RB
    @test actual.parameter == ["R", "L", "C"]
    @test actual.unit == ["Ω/m", "H/m", "F/m"]

    monte_carlo=TestFixtures.cable_monte_carlo_result()
    @test DataFrame(monte_carlo) == report(
        RB.MonteCarloTable(:kilo, nothing),
        monte_carlo
    ).table
    @test parentmodule(which(DataFrame, (typeof(monte_carlo),))) === RB
end
