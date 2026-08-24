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
    U.quantity(::typeof(profile_response)) = U.QuantityTag{:report_response}()
    U.native_unit(::U.QuantityTag{:report_response}) = U.units(:base, :ohm)
    U.display_unit(::U.QuantityTag{:report_response}) = U.units(:milli, :ohm)
    U.label(::U.QuantityTag{:report_response}) = "Response"
    U.symbol(::U.QuantityTag{:report_response}) = "u"

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
