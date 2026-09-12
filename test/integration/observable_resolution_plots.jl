@testitem "Makie / physical resolution is shared by ordinary and benchmark plots" tags=[:visual] begin
    using CairoMakie
    using LinearAlgebra: diag
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = [1.0, 1e3, 1e7]
    z = fill(1.0 + im, 1, 1, 3)
    y = fill(1e-13 + 1e-18im, 1, 1, 3)
    reference = LineParameters(PhaseDomain, z, y, f; details=(coordinates=["a"],))
    candidate = LineParameters(PhaseDomain, 2z, 2y, f; details=(coordinates=["a"],))
    options = (backend=:cairo, display_plot=false, controls=false,
        length_unit=:base, quantity_units=:base, open_export=false)
    publication = report(BenchmarkTableDefinition(quantities=(G, B, X)),
        (; reference, candidate))
    curves(page) = filter(plot -> plot isa Makie.Lines, only(page.axes).scene.plots)
    ordinates(page) = [last.(curve[1][]) for curve in curves(page)]
    ordinary = LineCableModels.plot(reference; ydata=(G,), options...)
    clean = LineCableModels.plot(publication; ydata=(G,), options...)
    @test all(iszero, only(ordinates(ordinary)))
    @test length(curves(clean)) == 2
    @test all(values -> all(iszero, values), ordinates(clean))
    @test only(clean.axes).subtitle[] == ""
    @test clean.addon_state.resolution.current_comparison
    raw = LineCableModels.plot(publication; ydata=(G,), clip=false, options...)
    @test first(ordinates(raw)) ≈ vec(Float32.(real.(y)))
    @test last(ordinates(raw)) ≈ 2vec(Float32.(real.(y)))
    @test observe(reference, Y) == y
    tight = report(BenchmarkTableDefinition(quantities=(G,), atol=(G=0.0,)),
        (; reference, candidate))
    inherited = LineCableModels.plot(tight; ydata=(G,), options...)
    @test ordinates(inherited) == ordinates(raw)
    @test !inherited.addon_state.resolution.display_override
    overridden = @test_logs (:warn, r"Plot resolution override") LineCableModels.plot(
        tight; ydata=(G,), atol=(G=1e-12,), options...)
    @test overridden.addon_state.resolution.display_override
    @test all(values -> all(iszero, values), ordinates(overridden))
    historical = merge(publication.published, (comparisons=map(publication.published.comparisons) do row
        old_details = Base.structdiff(row.error.details, (; resolution=row.error.details.resolution))
        old = LineCableModels.Engine.RMSError{Float64}(row.error.absolute, row.error.relative;
            details=old_details)
        merge(row, (error=old,))
    end,))
    history = @test_logs (:warn, r"Historical comparison semantics") LineCableModels.plot(
        historical; ydata=(G,), options...)
    @test !history.addon_state.resolution.current_comparison
    phase = LineCableModels.plot(reference; ydata=((Y, angle, 1, 1, :),), options...)
    @test all(isnan, only(ordinates(phase)))
    @test only(phase.axes).subtitle[] == "Undefined phase"
    standalone = LineCableModels.plot(ShuntAdmittance(y), f;
        ydata=((C, diag, :, :),), options...)
    @test all(iszero, only(ordinates(standalone)))
end
