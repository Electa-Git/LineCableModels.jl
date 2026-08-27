@testitem "ReportBuilder / grammar / publication and stage order" tags=[:unit] begin
    using DataFrames
    using RequiredInterfaces: NotImplementedError

    const RB = LineCableModels.ReportBuilder
    const PB = LineCableModels.PlotBuilder
    const U = LineCableModels.Units

    struct ReportProfile
        values::Vector{Float64}
    end
    profile_response(source::ReportProfile) = source.values
    const report_observation_calls = Ref(0)

    LineCableModels.basis(::ReportProfile) = :total
    LineCableModels.observe(
        source::ReportProfile,
        ::typeof(profile_response),
        indices...
    ) = begin
        report_observation_calls[] += 1
        isempty(indices) ? source.values : getindex(source.values, indices...)
    end
    LineCableModels.observables(::Type{<:ReportProfile}) = (profile_response,)
    U.quantity(::typeof(profile_response)) = U.Quantity{:report_response}()
    U.native_unit(::U.Quantity{:report_response}) = U.units(:base, :ohm)
    U.display_unit(::U.Quantity{:report_response}) = U.units(:milli, :ohm)
    U.label(::U.Quantity{:report_response}) = "Response"
    U.symbol(::U.Quantity{:report_response}) = "u"

    struct ReportProfilePlot <: PB.AbstractPlotDefinition end
    PB.entitle(::Type{ReportProfilePlot}, published::Tuple) = published
    function PB.resolve(
            ::Type{ReportProfilePlot},
            ::Tuple,
            request::NamedTuple
    )
        return request
    end
    function PB.fetch(
            ::Type{ReportProfilePlot},
            published::Tuple,
            request::NamedTuple
    )
        payload = (;
            published,
            legend = PB.LegendDefinition(enabled = false),
            colorbars = (),
            export_definition = PB.ExportDefinition(name = "report_profile")
        )
        return (PB.PlotPage(
            "Report profile",
            (800, 400),
            (; kind = :report_profile),
            payload
        ),)
    end

    source = ReportProfile([1.0, 2.0])
    definition = TableReportDefinition(
        (profile_response,);
        illustration = ReportProfilePlot
    )
    artifact = report(definition, source)

    @test report_observation_calls[] == 1
    @test artifact isa ReportArtifact
    @test parentmodule(typeof(artifact)) === RB
    @test artifact.table.u == [1_000.0, 2_000.0]
    columns=RB.observation_columns(artifact.table)
    @test U.label(columns.u.unit) == "mΩ"
    @test U.label(columns.u.quantity, columns.u.unit) == "Response [mΩ]"
    @test artifact.illustration isa PB.PlotRecipe
    @test only(only(artifact.illustration.pages).payload.published).values ==
          artifact.table.u
    @test artifact.output === nothing
    artifact.table.u[1] = 0.0
    @test source.values == [1.0, 2.0]
    @test_throws ArgumentError report(
        TableReportDefinition((identity,)),
        source
    )
    @test_throws MethodError report(
        TableReportDefinition((identity,)),
        :unsupported
    )
    publication_error=try
        LineCableModels.observables(source, (identity,))
        nothing
    catch error
        error
    end
    entitlement_error=try
        RB.entitle(TableReportDefinition((identity,)), source)
        nothing
    catch error
        error
    end
    @test typeof(entitlement_error) === typeof(publication_error)
    @test sprint(showerror, entitlement_error) == sprint(showerror, publication_error)

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
    @test_throws NotImplementedError report(MissingEncodeReport(), :source)

    struct MissingWriteReport <: RB.AbstractReportDefinition end
    RB.entitle(::MissingWriteReport, source) = source
    RB.select(::MissingWriteReport, source) = :published
    RB.tabulate(::MissingWriteReport, source, published) = :table
    RB.illustrate(::MissingWriteReport, source, published, table) = nothing
    RB.encode(::MissingWriteReport, source, published, table, illustration) = :encoded
    @test_throws NotImplementedError report(MissingWriteReport(), :source)
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
        artifact=report(XLSXReportDefinition(; file_name = path), parameters)
        @test artifact isa ReportArtifact
        @test fieldnames(typeof(artifact)) == (:table, :illustration, :output)
        @test artifact.illustration === nothing
        @test artifact.output == path
        @test isfile(artifact.output)
        @test artifact.table isa DataFrame
        @test names(artifact.table) == [
            "frequency", "row", "column", "R", "X", "G", "B"
        ]
        XLSX.openxlsx(path) do workbook
            @test XLSX.sheetnames(workbook) == [
                "Z(1,1)", "Z(1,2)", "Z(2,1)", "Z(2,2)",
                "Y(1,1)", "Y(1,2)", "Y(2,1)", "Y(2,2)"
            ]
            worksheet=workbook["Z(1,1)"]
            @test worksheet["A1"] == "frequency"
            @test worksheet["B1"] == "Hz"
            @test worksheet["A2"] == "R"
            @test worksheet["B2"] == "Ω/km"
            @test worksheet["A5"] == "frequency"
            @test worksheet["B5"] == "R"
            @test worksheet["C5"] == "X"
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
            XLSXReportDefinition(
                file_name = joinpath(directory, "named.xlsx"),
                cable_system = cable_system
            ),
            parameters
        )
        @test basename(prefixed.output) == "$(cable_system.system_id)_named.xlsx"

        default_artifact=cd(directory) do
            report(XLSXReportDefinition(), parameters)
        end
        @test default_artifact.output == joinpath(directory, "ZY_export.xlsx")
        @test isfile(default_artifact.output)
        @test !startswith(default_artifact.output, dirname(pathof(LineCableModels)))

        relative_artifact=cd(directory) do
            report(
                XLSXReportDefinition(file_name = "relative.xlsx"),
                parameters
            )
        end
        @test relative_artifact.output == joinpath(directory, "relative.xlsx")
        @test isfile(relative_artifact.output)

        diagonal=LineParameters(
            cat(Diagonal([1.0+2.0im, 2.0+3.0im]); dims = 3),
            cat(Diagonal([3.0+4.0im, 4.0+5.0im]); dims = 3),
            [50.0]
        )
        diagonal_artifact=@test_logs (:warn, r"Z is diagonal") (:warn, r"Y is diagonal") report(
            XLSXReportDefinition(file_name = joinpath(directory, "diagonal.xlsx")),
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
            XLSXReportDefinition(file_name = joinpath(directory, "uncertain.xlsx")),
            uncertain
        )
        XLSX.openxlsx(uncertain_artifact.output) do workbook
            @test occursin("±", workbook["Z(1,1)"]["B6"])
        end

        @test_throws Exception report(
            XLSXReportDefinition(file_name = directory),
            parameters
        )
    end

    @test !isdefined(IE, :XLSXWorkbook)
    @test !isdefined(IE, :_write_xlsx_sheet!)
    @test !isdefined(IE, :df_to_strings)
    @test !isdefined(RB, :_write_xlsx_sheet!)
    for private_name in (
        :_xlsx_string,
        :_xlsx_strings,
        :_xlsx_units,
        :_xlsx_destination
    )
        @test !isdefined(RB, private_name)
    end
    @test Base.ispublic(RB, :XLSXSheet)
    @test Base.ispublic(RB, :XLSXWorkbook)
    @test !isdefined(LineCableModels, :XLSXSheet)
    @test !isdefined(LineCableModels, :XLSXWorkbook)
end

@testitem "ReportBuilder / adapters / completed results delegate" tags=[:unit] setup=[
    TestFixtures
] begin
    using DataFrames

    const RB=LineCableModels.ReportBuilder
    constants=LineCableModels.CableConstants(1.0, 2.0, 3.0)
    expected=report(RB.CableConstantsTableDefinition(), constants).table
    actual=DataFrame(constants)

    @test actual == expected
    @test parentmodule(which(DataFrame, (typeof(constants),))) === RB
    @test names(actual) == ["R", "L", "C"]
    @test only(actual.R) == 1.0
    @test only(actual.L) == 2.0
    @test only(actual.C) == 3.0
    @test keys(RB.observation_columns(actual)) == (:R, :L, :C)

    monte_carlo=TestFixtures.cable_monte_carlo_result()
    @test DataFrame(monte_carlo) == report(
        RB.MonteCarloTableDefinition(:kilo, nothing),
        monte_carlo
    ).table
    @test parentmodule(which(DataFrame, (typeof(monte_carlo),))) === RB

    target=only(LineCableModels.Grammar.unit_targets(
        (R,),
        basis(monte_carlo);
        length_prefix = :kilo,
        overrides = :milli
    ))
    selected=RB.select(
        RB.MonteCarloTableDefinition(:kilo, :milli),
        monte_carlo
    )
    @test first(only(selected).observations).unit == target
    @test only(selected).frequency === nothing
    @test first(only(selected).observations).quantity ==
          LineCableModels.Units.quantity(R)

    for name in (
        :TableReportDefinition,
        :CableConstantsTableDefinition,
        :LineParametersTableDefinition,
        :BenchmarkTableDefinition,
        :MonteCarloTableDefinition,
        :XLSXReportDefinition
    )
        @test isdefined(RB, name)
        @test parentmodule(getproperty(RB, name)) === RB
    end
    for retired in (
        :TableReport,
        :CableConstantsTable,
        :LineParametersTable,
        :BenchmarkTable,
        :MonteCarloTable,
        :XLSXReport
    )
        @test !isdefined(RB, retired)
        @test !isdefined(LineCableModels, retired)
    end
end
