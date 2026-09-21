@testitem "Makie / physical resolution is shared by ordinary and benchmark plots" tags=[:visual] begin
    using CairoMakie
    using LinearAlgebra: diag
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = [1.0, 1e3, 1e7]
    z = fill(1.0 + im, 1, 1, 3)
    y = fill(1e-13 + 1e-18im, 1, 1, 3)
    reference = LineParameters(PhaseDomain, z, y, f; details=ComputationDetails(;coordinates=["a"],))
    candidate = LineParameters(PhaseDomain, 2z, 2y, f; details=ComputationDetails(;coordinates=["a"],))
    using LineCableModels.Engine: retain_gridpoint
    using LineCableModels.Grammar: gridpoint_id
    reference=retain_gridpoint(reference,gridpoint_id())
    candidate=retain_gridpoint(candidate,gridpoint_id())
    options=(backend=:cairo,display_plot=false,controls=false,open_export=false)
    curves(page)=filter(plot -> plot isa Makie.Lines,only(page.axes).scene.plots)
    ordinates(page)=[last.(curve[1][]) for curve in curves(page)]
    publication=report(BenchmarkTableDefinition(quantities=(G,B,X),clip=true),
        (;reference,candidate);observation_options=(length_unit=:base,quantity_units=:base))
    ordinary=LineCableModels.plot(reference;ydata=(G,),length_unit=:base,options...)
    clean=LineCableModels.plot(publication;ydata=(G,),options...)
    @test all(iszero,only(ordinates(ordinary)))
    @test length(curves(clean))==2
    @test all(values -> all(iszero,values),ordinates(clean))
    @test_throws ArgumentError LineCableModels.plot(publication;ydata=(G,),clip=false,options...)
    raw=report(BenchmarkTableDefinition(quantities=(G,),clip=false),(;reference,candidate);
        observation_options=(length_unit=:base,quantity_units=:base))
    rawpage=LineCableModels.plot(raw;ydata=(G,),options...)
    @test first(ordinates(rawpage))≈2vec(Float32.(real.(y)))
    @test last(ordinates(rawpage))≈vec(Float32.(real.(y)))
    @test observe(reference,Y)==y
    phase = LineCableModels.plot(reference; ydata=((Y, angle, 1, 1, :),), options...)
    @test all(isnan, only(ordinates(phase)))
    @test only(phase.axes).subtitle[] == "Undefined phase"
    axisscale!(phase,:x,log10)
    axisscale!(phase,:y,log10)
    @test only(phase.axes).xscale[]===log10
    @test only(phase.axes).yscale[]===log10
    @test all(>(0),only(phase.axes).targetlimits[].origin)
    resetview!(phase)
    @test all(isfinite,only(phase.axes).targetlimits[].widths)
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
        errorbar_sampling=:all,fig_size=(900,500))
    # Independent means/spreads reproduce the noisy-zero and finite-baseline
    # failures. Real uncertainty at zero must still produce correctly centred bars.
    for (means, spreads, clean_means, clean_spreads) in (
            (range(-1e-27, 2e-27; length=13), fill(4e-27, 13), zeros(13), fill(4e-27,13)),
            (fill(-7e-10, 13) .+ (0:12) .* 1e-25,
                fill(3e-25, 13), nothing, fill(3e-25,13)),
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
