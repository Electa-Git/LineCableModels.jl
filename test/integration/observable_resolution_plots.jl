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
    phase = LineCableModels.plot(reference; ydata=((Y, angle, 1, 1, :),), options...)
    @test all(isnan, only(ordinates(phase)))
    @test only(phase.axes).subtitle[] == "Undefined phase"
    standalone = LineCableModels.plot(ShuntAdmittance(y), f;
        ydata=((C, diag, :, :),), options...)
    @test all(iszero, only(ordinates(standalone)))
end

@testitem "Makie / uncertainty publication removes residue without changing raw curves" tags=[:visual] begin
    using CairoMakie, Measurements, Logging
    f = 10.0 .^ range(-1, 7; length=13)
    omega = reshape(2pi .* f, 1, 1, :)
    options = (backend=:cairo, display_plot=false, controls=true,
        length_unit=:base, quantity_units=:base, open_export=false,
        signed_ylog=true, fig_size=(900, 500))
    # Independent means/spreads reproduce the noisy-zero and finite-baseline
    # failures. Real uncertainty at zero must still produce correctly centred bars.
    for (means, spreads, clean_means, clean_spreads) in (
            (range(-1e-27, 2e-27; length=13), fill(4e-27, 13), zeros(13), zeros(13)),
            (fill(-7e-10, 13) .+ (0:12) .* 1e-25,
                fill(3e-25, 13), nothing, zeros(13)),
            (fill(1e-27, 13), fill(4e-10, 13), zeros(13), fill(4e-10, 13)))
        c = reshape(measurement.(means, spreads), 1, 1, :)
        source = LineParameters(one.(c) .+ im .* one.(c), im .* omega .* c, f)
        before = deepcopy((Z(source), Y(source), frequencies(source)))
        expected = clean_means === nothing ? nominal.(vec(observe(source, C))) : clean_means
        for raw in (false, true)
            logger = Test.TestLogger(min_level=Logging.Warn)
            page = with_logger(logger) do
                page = LineCableModels.plot(source, source; ydata=(C,), clip=!raw, options...)
                Makie.colorbuffer(page.figure)
                for _ in 1:2
                    page.controls[:ylog].active[] = true
                    Makie.colorbuffer(page.figure)
                    page.controls[:ylog].active[] = false
                    page.controls[:reset].clicks[] += 1
                end
                page
            end
            @test isempty(logger.logs)
            axis = only(page.axes)
            lines = filter(p -> p isa Makie.Lines, axis.scene.plots)
            bars = filter(p -> p isa Makie.Errorbars, axis.scene.plots)
            @test length(lines) == 2
            @test all(line -> last.(line[1][]) ≈ (raw ? nominal.(vec(c)) : expected), lines)
            wanted_spread = raw ? vec(uncertainty.(c)) : clean_spreads
            @test isempty(bars) == all(iszero, wanted_spread)
            for bar in bars
                @test [point[2] for point in bar[1][]] ≈ (raw ? nominal.(vec(c)) : expected)
                @test [point[3] for point in bar[1][]] ≈ wanted_spread
                @test [point[4] for point in bar[1][]] ≈ wanted_spread
            end
            if !raw && all(iszero, expected) && all(iszero, clean_spreads)
                @test axis.finallimits[].origin[2] == -1
                @test axis.finallimits[].widths[2] == 2
            elseif clean_means === nothing
                @test axis.finallimits[].widths[2] >= 0.099abs(first(expected))
            end
        end
        @test isequal((Z(source), Y(source), frequencies(source)), before)
    end
end

@testitem "Makie / a new observation owner inherits uncertainty projection" tags=[:visual] begin
    using CairoMakie, Measurements, Logging
    using LineCableModels.Grammar: observation_resolution
    struct ResolutionSource{T}
        samples::T
    end
    LineCableModels.Grammar.basis(::ResolutionSource) = :pul
    LineCableModels.Grammar.observables(::Type{<:ResolutionSource}) = (G,)
    LineCableModels.Grammar.observe(source::ResolutionSource, ::typeof(G)) = source.samples
    function LineCableModels.Grammar.observation_resolution(source::ResolutionSource, request;
            atol=nothing, frequencies=nothing)
        return observation_resolution(observe(source, request), G; atol, frequencies)
    end
    source = ResolutionSource(measurement.([1e-18, 1e-6], [2e-18, 2e-18]))
    original = copy(source.samples)
    options = (backend=:cairo, display_plot=false, controls=false, open_export=false)
    for clip in (true, false)
        publication = observables(source, (G,); clip, length_unit=:base)
        logger = Test.TestLogger(min_level=Logging.Warn)
        page = with_logger(logger) do
            page = LineCableModels.plot(publication; options...)
            Makie.colorbuffer(page.figure)
            page
        end
        @test isempty(logger.logs)
        axis = only(page.axes)
        line = only(filter(p -> p isa Makie.Lines, axis.scene.plots))
        @test last.(line[1][]) ≈ (clip ? [0, 1e-6] : nominal.(original))
        @test isempty(filter(p -> p isa Makie.Errorbars, axis.scene.plots)) == clip
        @test source.samples == original
    end
end
