@testitem "Makie addons / narrow log values and native constructors share the shell" tags=[:visual] begin
    using CairoMakie, Measurements, Logging
    f = collect(range(2.0,8.0;length=9))
    r = reshape(measurement.(reverse(f),0.02),1,1,:)
    source = LineParameters(complex.(r,r),ones(ComplexF64,1,1,9),f)
    before = deepcopy((Z(source),Y(source),frequencies(source)))
    logger = Test.TestLogger(min_level=Logging.Warn)
    with_logger(logger) do
        page = LineCableModels.plot(source; ydata=(R,),length_unit=:base,
            backend=:cairo,display_plot=false,open_export=false,
            figure=(size=(800,500),),axis=(limits=((2.0,8.0),(1.9,8.1)),),linewidth=3)
        axis = only(page.axes)
        line = only(filter(p -> p isa Makie.Lines,axis.scene.plots))
        @test line.linewidth[] == 3
        requested = axis.limits[]
        for _ in 1:2
            page.controls[:xlog].active[] = true
            page.controls[:ylog].active[] = true
            Makie.colorbuffer(page.figure)
            @test page.controls[:xlog].active[] && page.controls[:ylog].active[]
            for dim in (:x,:y)
                native = getproperty(axis,Symbol(dim,:axis))
                @test length(native.tickvalues[]) >= 3
                @test all(label -> label isa AbstractString,native.ticklabels[])
                @test parse.(Float64,native.ticklabels[]) ≈ native.tickvalues[]
            end
            # Raw points staying equal cannot detect a failed transform. Compare
            # rendered geometry with the log of these distinguishable values.
            raw = Makie.Point2d.(line[1][])
            pixels = Makie.transform_and_project(line,:data,:pixel,raw)
            for dim in 1:2
                projected = getindex.(pixels,dim)
                values = getindex.(raw,dim)
                @test (projected.-first(projected))./(last(projected)-first(projected)) ≈
                    (log10.(values).-log10(first(values)))./(log10(last(values))-log10(first(values))) atol=2e-6
            end
            for bars in filter(p -> p isa Makie.Errorbars,axis.scene.plots), child in bars.plots
                @test child.transformation.transform_func[] == line.transformation.transform_func[]
                # Check rendered bar endpoints as well as their parent scale;
                # a pre-logged child could otherwise pass a transform-only check.
                endpoints = filter(point -> all(isfinite,point),Makie.Point2d.(child[1][]))
                projected = Makie.transform_and_project(child,:data,:pixel,endpoints)
                limits = axis.finallimits[]
                viewport = axis.scene.viewport[]
                for dimension in 1:2
                    lower = limits.origin[dimension]
                    upper = lower+limits.widths[dimension]
                    fractions = (log10.(getindex.(endpoints,dimension)).-log10(lower))./
                        (log10(upper)-log10(lower))
                    @test getindex.(projected,dimension)./viewport.widths[dimension] ≈ fractions atol=2e-6
                end
            end
            @test axis.limits[] == requested
            page.controls[:xlog].active[] = false
            page.controls[:ylog].active[] = false
        end
        page.controls[:ylog].active[] = true
        for bounds in ((2.71,2.79),(2e-12,8e-12),(2e100,8e100))
            ylims!(axis,bounds...)
            Makie.colorbuffer(page.figure)
            @test length(axis.yaxis.tickvalues[]) >= 2
            @test allunique(axis.yaxis.ticklabels[])
            @test all(label -> label isa AbstractString,axis.yaxis.ticklabels[])
        end
        ylims!(axis,0.1,1e7)
        Makie.colorbuffer(page.figure)
        @test all(value -> isapprox(log10(Float64(value)),round(log10(Float64(value)));atol=1e-7),axis.yaxis.tickvalues[])
        @test !occursin("× 10",repr(axis.ylabel[]))
        for dim in (:x,:y), bounds in ((0.6,3.0),(2.0,80.0),(2.71,2.79),(8.0,12.0))
            getproperty(axis,Symbol(dim,:scale))[] = log10
            (dim === :x ? xlims! : ylims!)(axis,bounds...)
            Makie.colorbuffer(page.figure)
            native = getproperty(axis,Symbol(dim,:axis))
            @test all(label -> label isa AbstractString,native.ticklabels[])
            # Measure actual strings against their native projected positions.
            # Density may change; collisions or identical labels may not.
            dimension = dim === :x ? 1 : 2
            fontsize = getproperty(axis,Symbol(dim,:ticklabelsize))[]
            probe = text!(axis.blockscene,0,0;visible=false,fontsize,
                font=getproperty(axis,Symbol(dim,:ticklabelfont))[],
                rotation=getproperty(axis,Symbol(dim,:ticklabelrotation))[])
            extents = map(native.ticklabels[]) do label
                probe.text[] = label
                Makie.boundingbox(probe,:data).widths[dimension]
            end
            @test allunique(native.ticklabels[])
            @test all(diff(getindex.(native.tickpositions[],dimension)) .>=
                (extents[1:end-1].+extents[2:end])./2 .+fontsize/2 .-1)
            delete!(axis.blockscene,probe)
        end
        ylims!(axis,2e-12,8e-12)
        axis.yticks[] = ([3e-12,7e-12],["low","high"])
        Makie.colorbuffer(page.figure)
        @test axis.yaxis.ticklabels[] == ["low","high"]
        @test !occursin("× 10",repr(axis.ylabel[]))
        mktempdir() do directory
            @test isfile(export_svg(page;path=joinpath(directory,"narrow.svg"),open_file=false))
            @test axis.yaxis.ticklabels[] == ["low","high"]
        end
    end
    @test isempty(logger.logs)
    @test isequal(before,(Z(source),Y(source),frequencies(source)))
end

@testitem "Makie addons / current visible extents control log eligibility" tags=[:visual] begin
    using CairoMakie
    page = LineCableModels.plotwindow(;title="Native",backend=:cairo,
        display_plot=false,open_export=false) do layout
        axis = Axis(layout[1,1])
        lines!(axis,[1.0,2.0],[2.0,2.0])
        lines!(axis,[1.0,2.0],[-1.0,-2.0])
    end
    axis = only(page.axes)
    @test haskey(page.controls,:xlog) && haskey(page.controls,:ylog)
    previous = axis.targetlimits[]
    page.controls[:ylog].active[] = true
    @test !page.controls[:ylog].active[]
    @test axis.yscale[] === identity
    @test axis.targetlimits[] == previous
    @test occursin(repr(axis.title[]),page.addon_state.shell.status[])
    last(axis.scene.plots).visible[] = false
    page.controls[:ylog].active[] = true
    Makie.colorbuffer(page.figure)
    @test axis.yscale[] === log10
    low = axis.finallimits[].origin[2]
    high = low+axis.finallimits[].widths[2]
    @test 1.8 < low < 2 < high < 2.2
end

@testitem "Makie addons / native override precedence and custom transforms" tags=[:visual] begin
    using CairoMakie
    source = LineParameters(reshape(ComplexF64[2,4,8],1,1,:),ones(ComplexF64,1,1,3),[2.0,4.0,8.0])
    page = LineCableModels.plot(source; ydata=R,backend=:cairo,display_plot=false,
        open_export=false,length_unit=:base,xscale=sqrt,
        figure=(size=(500,700),),ylabel="Shared",axis=(ylabel="Explicit",),
        linewidth=3,series_attributes=((linewidth=5,),))
    axis = only(page.axes)
    @test axis.xscale[] === sqrt
    @test !haskey(page.controls,:xlog)
    @test axis.ylabel[] == "Explicit"
    @test Tuple(page.figure.scene.viewport[].widths) == (500,700)
    @test only(filter(p->p isa Makie.Lines,axis.scene.plots)).linewidth[] == 5
    @test_throws ArgumentError LineCableModels.plot(source;ydata=R,backend=:cairo,
        display_plot=false,not_a_native_attribute=true)
    # Reject a second public selection grammar, without banning legitimate
    # local variables called quantities inside the implementation.
    @test_throws ArgumentError LineCableModels.plot(source;backend=:cairo,
        display_plot=false,quantities=(R,))
    @test_throws ArgumentError LineCableModels.plotwindow(_ -> nothing;
        title="No axis",backend=:cairo,display_plot=false,axis=(xlabel="Unused",))

    native = LineCableModels.plotwindow(;title="Native overrides",backend=:cairo,
        display_plot=false,open_export=false,axis=(xscale=log10,ylabel="Explicit"),
        linewidth=4) do grid
        axis = Axis(grid[1,1];xlabel="Native label",ylabel="Native label",limits=(2,8,2,8))
        lines!(axis,[2.0,4.0,8.0],[2.0,4.0,8.0])
    end
    axis = only(native.axes)
    @test axis.xscale[] === log10
    @test axis.xlabel[] == "Native label"
    @test axis.ylabel[] == "Explicit"
    @test axis.limits[] == (2,8,2,8)
    @test first(axis.scene.plots).linewidth[] == 4
    @test native.controls[:xlog].active[]
    @test !native.controls[:ylog].active[]
    @test !isempty(Makie.colorbuffer(native.figure))
end

@testitem "Makie addons / axis multipliers follow the displayed range" tags=[:visual] begin
    using CairoMakie

    frequency = [1.0, 10.0, 100.0]
    impedance = fill(1.0 + 2.0im, 1, 1, 3)
    susceptance = reshape([-1.0, 0.0, 1.0] .* 1e-14, 1, 1, 3)
    parameters = LineParameters(impedance, complex.(zero.(susceptance), susceptance), frequency)
    plot = Makie.plot(parameters, B; backend=:cairo, display_plot=false,
        controls=false, length_unit=:base, quantity_units=:base, clip=false)
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    # Tiny values are not constant merely because their SI magnitude is small.
    @test axis.finallimits[].widths[2] < 1e-12
    @test axis.ytickformat[]([-1e-14, 0.0, 1e-14]) == ["-10", "0", "10"]

    ylims!(axis, -2e-8, 2e-8)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([-1e-8, 0.0, 1e-8]) == ["-10", "0", "10"]
    @test occursin("−9", repr(axis.ylabel[]))
    @test !occursin("−15", repr(axis.ylabel[]))
    axis.ylabel[] = "Caller response [S/m]"
    @test occursin("Caller response", repr(axis.ylabel[]))
    @test occursin("−9", repr(axis.ylabel[]))
    resize!(plot.figure, 1000, 600)
    Makie.colorbuffer(plot.figure)
    label = repr(axis.ylabel[])
    ticks = axis.ytickformat[]([-1e-8, 0.0, 1e-8])
    mktempdir() do directory
        path = export_svg(plot; path=joinpath(directory, "axis.svg"), open_file=false)
        @test isfile(path)
        @test repr(axis.ylabel[]) == label
        @test axis.ytickformat[]([-1e-8, 0.0, 1e-8]) == ticks
    end
    @test Y(parameters) == complex.(zero.(susceptance), susceptance)
    publication = observables(parameters, ((@observe B[1, 1, :]),);
        length_unit=:base, quantity_units=:base, clip=false)
    published_plot = Makie.plot(publication; backend=:cairo, display_plot=false,
        controls=false)
    published_axis = only(published_plot.axes)
    ylims!(published_axis, -2e-9, 2e-9)
    Makie.colorbuffer(published_plot.figure)
    @test published_axis.ytickformat[]([-1e-9, 0.0, 1e-9]) == ["-1", "0", "1"]
    @test occursin("−9", repr(published_axis.ylabel[]))
end

@testitem "Makie addons / near-constant views and stable signed scales share the lifecycle" tags=[:visual] begin
    using CairoMakie, Logging
    options = (backend=:cairo, display_plot=false, controls=true, open_export=false,
        length_unit=:base, quantity_units=:base, clip=false, signed_ylog=true)
    f = [1.0, nextfloat(1.0), nextfloat(nextfloat(1.0))]
    response = [-3.592027795269888e-10, -3.592027795269885e-10, -3.592027795269883e-10]
    source = LineParameters(reshape(response .+ im, 1, 1, :), ones(ComplexF64, 1, 1, 3), f)
    logger = Test.TestLogger(min_level=Logging.Warn)
    page = with_logger(logger) do
        page = LineCableModels.plot(source; ydata=(R,), xscale=:linear, options...)
        Makie.colorbuffer(page.figure)
        page
    end
    @test isempty(logger.logs)
    axis = only(page.axes)
    @test axis.finallimits[].widths[1] ≈ 0.1
    @test axis.finallimits[].widths[2] ≈ 0.1abs(first(response))
    @test axis.limits[] == (nothing, nothing)
    # Native zoom and explicit narrow views remain possible; padding belongs
    # only to automatic fitting, not to tick updates or every targetlimits write.
    manual = Makie.Rect2d(1.0, minimum(response), 1e-12, 1e-20)
    axis.targetlimits[] = manual
    @test axis.targetlimits[] == manual
    ylims!(axis, minimum(response), maximum(response))
    @test axis.targetlimits[].widths[2] < 1e-22
    axis.limits[] = ((0.9, nothing), (nothing, -3e-10))
    page.controls[:reset].clicks[] += 1
    @test axis.targetlimits[].origin[1] == 0.9
    @test axis.targetlimits[].origin[2] + axis.targetlimits[].widths[2] ≈ -3e-10
    autolimits!(axis)
    with_logger(logger) do
        for _ in 1:2
            page.controls[:ylog].active[] = true
            scale = axis.yscale[]
            inverse = Makie.inverse_transform(scale)
            # Numerical behavior, not identity of a dependency's unstable scale,
            # is the contract. It must work both near zero and over large decades.
            for value in (-1e100, -1.0, -1e-18, 0.0, 1e-18, 1.0, 1e100)
                @test inverse(scale(value)) ≈ value rtol=1e-13
                @test sign(scale(value)) == sign(value)
            end
            for bounds in ((-2e-18, 2e-18), (-1000.0, 1000.0))
                ylims!(axis, bounds...)
                Makie.colorbuffer(page.figure)
                @test length(axis.yaxis.tickvalues[]) >= 2
                @test allunique(axis.yaxis.ticklabels[])
                @test all(point -> all(isfinite, point), axis.yaxis.tickpositions[])
            end
            autolimits!(axis)
            page.controls[:ylog].active[] = false
            page.controls[:reset].clicks[] += 1
        end
        mktempdir() do directory
            @test isfile(export_svg(page; path=joinpath(directory, "resolved.svg"), open_file=false))
        end
    end
    @test isempty(logger.logs)
    @test vec(real.(Z(source))) == response
end

@testitem "Makie addons / rendered benchmark axes share scale and limit behavior" tags=[:visual] begin
    using CairoMakie, Statistics, Measurements
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = 10.0 .^ range(-1, 7; length=101)
    omega = reshape(2pi .* f, 1, 1, :)
    r = [(i + 2j) * 1e-3 * (1 + log10(1 + x)) for i in 1:2, j in 1:2, x in f]
    v = (R=r, L=r .* 1e-4, C=r .* 1e-8, G=r .* 1e-5)
    core = LineParameters(v.R .+ im .* omega .* v.L, v.G .+ im .* omega .* v.C, f)
    metadata = (port_order=["a", "b"], formulation=NamedTuple(Formulation()), axes=nothing)
    deterministic = report(BenchmarkTableDefinition((R,); bands=(:all,)),
        (reference=(result=core, metadata=metadata), candidate=(result=core, metadata=metadata)))
    summaries = map(a -> map(x -> SampleSummary([0.9x, 1.1x]), a), v)
    mc = MonteCarloResult(MonteCarlo(Formulation(); trials=2, seed=7), [core],
        [summaries], nothing, nothing, UInt64(7), UInt64[8], [2])
    m = map(a -> measurement.(a, sqrt(2) * 0.1 .* a), v)
    lep = LinearErrorResult(LinearError(Formulation()), [LineParameters(
        m.R .+ im .* omega .* m.L, m.G .+ im .* omega .* m.C, measurement.(f, 0.0))])
    uq = report(BenchmarkTableDefinition(((statistics, R, mean), (statistics, R, std)); bands=(:all,)),
        (reference=(result=mc, metadata=metadata), candidate=(result=lep, metadata=metadata)))
    for (artifact, request, expected) in ((deterministic, (R,), r),
            (uq, ((statistics, R, mean),), r),
            (uq, ((statistics, R, std),), sqrt(2) * 0.1 .* r))
        page = LineCableModels.plot(artifact; ydata=request, backend=:cairo,
            display_plot=false, controls=true, open_export=false, clip=false,
            length_unit=:base, quantity_units=:base, fig_size=(1100, 750))
        Makie.colorbuffer(page.figure)
        axis = last(page.axes)
        @test length(axis.xaxis.tickvalues[]) >= 8
        page.controls[:xlog].active[] = false
        Makie.colorbuffer(page.figure)
        # Test what Makie rendered, not a formatter evaluated on invented ticks.
        @test 5 <= length(axis.xaxis.tickvalues[]) <= 10
        @test occursin("6", repr(axis.xlabel[]))
        @test parse.(Float64, axis.xaxis.ticklabels[]) .* 1e6 ≈ axis.xaxis.tickvalues[]
        for ((i, j), panel) in pairs(page.addon_state.panel_data)
            for line in filter(p -> p isa Makie.Lines, panel.axis.scene.plots)
                @test first.(line[1][]) ≈ f
                @test last.(line[1][]) ≈ expected[i, j, :]
            end
        end
        xlims!(axis, 1e6, 4e6)
        ylims!(axis, 0.004, 0.012)
        axis.xticks[] = [1e6, 3e6]
        requested = axis.limits[]
        for _ in 1:2
            page.controls[:xlog].active[] = true
            page.controls[:xlog].active[] = false
            Makie.colorbuffer(page.figure)
            @test axis.limits[] == requested
            @test axis.xticks[] == [1e6, 3e6]
            @test axis.finallimits[].origin[2] ≈ 0.004
            @test axis.finallimits[].widths[2] ≈ 0.008
        end
        axis.xticks[] = Makie.automatic
        page.controls[:xlog].active[] = true
        xlims!(axis, 2.0, 8.0)
        Makie.colorbuffer(page.figure)
        @test length(axis.xaxis.tickvalues[]) >= 2
        @test all(x -> 2 <= x <= 8, axis.xaxis.tickvalues[])
        @test !occursin("× 10", repr(axis.xlabel[]))
        # A sub-decade LogTicks locator, like PseudologTicks, is specific to its
        # transform. It must not leak into native linear tick conversion.
        page.controls[:xlog].active[] = false
        @test length(axis.xaxis.tickvalues[]) >= 2
    end
    @test R(core) == r
    @test frequencies(core) == f
    @test isequal(observe(lep, statistics, R, std, 1), sqrt(2) * 0.1 .* r)
end

@testitem "Makie addons / live native overrides and requested limits survive" tags=[:visual] begin
    using CairoMakie
    page = LineCableModels.plotwindow(; title="Live overrides", backend=:cairo,
        display_plot=false, controls=true, open_export=false) do layout
        axis = Axis(layout[1, 1]; xlabel="Frequency [Hz]", ylabel="Response",
            limits=((1e6, 4e6), (nothing, 0.8)))
        lines!(axis, [0.1, 1e7], [0.0, 1.0])
    end
    axis = only(page.axes)
    Makie.colorbuffer(page.figure)
    @test axis.limits[] == ((1e6, 4e6), (nothing, 0.8))
    @test axis.finallimits[].origin[1] ≈ 1e6
    @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) ≈ 0.8
    # A labelled tuple must release the managed formatter BEFORE native tick
    # conversion sees the new value; a late observer cannot repair that error.
    labelled = ([1e6, 3e6], ["low", "high"])
    axis.xticks[] = labelled
    Makie.colorbuffer(page.figure)
    @test axis.xaxis.ticklabels[] == labelled[2]
    @test axis.xlabel[] == "Frequency [Hz]"
    axis.xticks[] = (lo, hi) -> labelled
    Makie.colorbuffer(page.figure)
    @test axis.xaxis.ticklabels[] == labelled[2]
    @test axis.xlabel[] == "Frequency [Hz]"
    axis.xticks[] = Makie.automatic
    axis.xtickformat[] = "{:.1f}"
    @test axis.xlabel[] == "Frequency [Hz]"
    axis.xtickformat[] = Makie.automatic
    @test occursin("× 10", repr(axis.xlabel[]))
    page.controls[:reset].clicks[] += 1
    @test axis.limits[] == ((1e6, 4e6), (nothing, 0.8))
end

@testitem "Makie addons / scale changes validate the page and preserve orthogonal views" tags=[:visual] begin
    using CairoMakie
    using LineCableModels.ReportBuilder: BenchmarkTableDefinition
    f = [1.0, 10.0, 100.0]
    z = [complex(i+j+k, k) for i in 1:2, j in 1:2, k in 1:3]
    data = LineParameters(z, z .* 1e-6, f)
    page = LineCableModels.plot(data; ydata=(R,), backend=:cairo, display_plot=false,
        controls=true, open_export=false, length_unit=:base, xscale=:linear)
    first_axis, last_axis = first(page.axes), last(page.axes)
    # Model a native pan/zoom: target limits change, the stored request does not.
    first_axis.targetlimits[] = Makie.Rect2d(2.0, 3.0, 30.0, 2.0)
    page.controls[:xlog].active[] = true
    @test first_axis.limits[] == (nothing, nothing)
    @test first_axis.targetlimits[].origin[2] == 3.0
    @test first_axis.targetlimits[].widths[2] == 2.0
    page.controls[:xlog].active[] = false
    xlims!(last_axis, -1.0, 100.0)
    old_views = [axis.targetlimits[] for axis in page.axes]
    # The invalid request is on the LAST axis: validating only while applying
    # changes would leave the earlier axes switched when the exception is raised.
    page.controls[:xlog].active[] = true
    @test occursin("positive",page.addon_state.shell.status[])
    @test !page.controls[:xlog].active[]
    @test all(axis -> axis.xscale[] === identity, page.axes)
    @test [axis.targetlimits[] for axis in page.axes] == old_views
    autolimits!(last_axis)
    page.controls[:xlog].active[] = true
    page.controls[:xlog].active[] = false
    first_axis.limits[] = ((5.0, nothing), (nothing, 10.0))
    page.controls[:reset].clicks[] += 1
    @test first_axis.targetlimits[].origin[1] == 5.0
    @test sum((first_axis.targetlimits[].origin[2], first_axis.targetlimits[].widths[2])) ≈ 10.0
    @test Z(data) == z

    negative = LineParameters(z, fill(-1e-6 + 1e-6im, 2, 2, 3), f)
    operand = (result=negative,
        metadata=(port_order=["a", "b"], formulation=NamedTuple(Formulation()), axes=nothing))
    comparison = report(BenchmarkTableDefinition((G,)), (reference=operand, candidate=operand))
    signed = LineCableModels.plot(comparison; ydata=(G,), backend=:cairo,
        display_plot=false, controls=true, open_export=false, yscale=:pseudolog10)
    @test signed.controls[:ylog].active[]
    signed.controls[:ylog].active[] = false
    signed.controls[:ylog].active[] = true
    @test all(axis -> axis.yscale[](1e-18) > 0 &&
        axis.yscale[](-1e-18) < 0, signed.axes)
end

@testitem "Makie addons / native axis density responds to size without losing precision" tags=[:visual] begin
    using CairoMakie
    page = LineCableModels.plotwindow(; title="Numeric ranges", backend=:cairo,
        display_plot=false, controls=false, open_export=false, size=(1200, 600)) do layout
        axis = Axis(layout[1, 1]; xlabel="Input", ylabel="Output")
        lines!(axis, [0.1, 1e7], [1e8, 1e8 + 2])
    end
    axis = only(page.axes)
    Makie.colorbuffer(page.figure)
    large_count = length(axis.xaxis.tickvalues[])
    @test large_count >= 6
    resize!(page.figure, 600, 400)
    Makie.colorbuffer(page.figure)
    @test 3 <= length(axis.xaxis.tickvalues[]) < large_count
    # Subnormal mantissa formatting is exercised below, independently of Makie's
    # camera: its reciprocal viewport scale cannot represent a subnormal span.
    for bounds in ((-2e-18, 2e-18), (1e8, 1e8 + 2), (-2e300, 2e300),
            (-2e-300, 2e-300))
        xlims!(axis, bounds...)
        Makie.colorbuffer(page.figure)
        @test length(axis.xaxis.ticklabels[]) >= 2
        @test allunique(axis.xaxis.ticklabels[])
        @test all(label -> !occursin(r"[eE]", label), axis.xaxis.ticklabels[])
        @test all(point -> all(isfinite, point), axis.xaxis.tickpositions[])
    end
    xlims!(axis, 1e8, 1e8 + 2)
    axis.xticklabelsize[] = 20
    Makie.colorbuffer(page.figure)
    @test allunique(axis.xaxis.ticklabels[])
    before = (axis.limits[], axis.targetlimits[], repr(axis.xlabel[]), copy(axis.xaxis.ticklabels[]))
    mktempdir() do directory
        for theme in (:default, :publication)
            export_svg(page; path=joinpath(directory,"$theme.svg"), theme, open_file=false)
            Makie.colorbuffer(page.figure)
            @test (axis.limits[], axis.targetlimits[], repr(axis.xlabel[]), axis.xaxis.ticklabels[]) == before
        end
    end
end

@testitem "Makie addons / plain mantissas and relative constant limits" tags=[:visual] begin
    using CairoMakie
    using Measurements: measurement
    extension = Base.get_extension(LineCableModels, :LineCableModelsMakieExt)
    for exponent in (-300, -14, -3, 0, 8, 300)
        format = extension._addon_linear_tickformat(exponent)
        @test format([-1.0, 0.0, 1.0] .* 10.0^exponent) == ["-1", "0", "1"]
        @test format([-0.0, 0.25, 0.5] .* 10.0^exponent) == ["0", "0.25", "0.5"]
    end
    @test extension._addon_linear_tickformat(-324)([0.0, nextfloat(0.0)]) ==
          ["0", "4.941"]
    format = extension._addon_linear_tickformat(8)
    ticks = [1e8, 1e8 + 1, 1e8 + 2]
    @test length(unique(format(ticks))) == 3
    @test all(!occursin(r"[eE]", label) for label in format(ticks))

    for (value, bounds) in (
        (0.0, (-1.0, 1.0)),
        (1e-14, (0.95e-14, 1.05e-14)),
        (-1e-14, (-1.05e-14, -0.95e-14)),
        (measurement(0.0, 1e-14), (-2e-14, 2e-14)),
        (measurement(1e-14, 2e-15), (6e-15, 1.4e-14)),
    )
        resistance = fill(value, 1, 1, 3)
        parameters = LineParameters(complex.(resistance, one.(resistance)),
            fill(1.0im, 1, 1, 3), [1.0, 10.0, 100.0])
        plot = Makie.plot(parameters, R; backend=:cairo, display_plot=false,
            controls=false, length_unit=:base, quantity_units=:base, clip=false)
        axis = only(plot.axes)
        Makie.colorbuffer(plot.figure)
        @test axis.finallimits[].origin[2] ≈ bounds[1]
        @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) ≈ bounds[2]
    end
end

@testitem "Makie addons / native conversions and visible uncertainty retain ownership" tags=[:visual] begin
    using CairoMakie, Dates, Measurements
    for input in (Makie.Categorical(["first", "second", "third"]),
            [Date(2026, 1, 1), Date(2026, 2, 1), Date(2026, 3, 1)])
        page = LineCableModels.plotwindow(; title="Native conversion", backend=:cairo,
            display_plot=false, controls=true, open_export=false) do layout
            axis = Axis(layout[1, 1]; xlabel="Native coordinate")
            scatter!(axis, input, [1.0, 2.0, 3.0])
        end
        axis = only(page.axes)
        Makie.colorbuffer(page.figure)
        before = copy(axis.xaxis.ticklabels[])
        page.controls[:reset].clicks[] += 1
        Makie.colorbuffer(page.figure)
        @test axis.xlabel[] == "Native coordinate"
        @test !isempty(before)
        @test axis.xaxis.ticklabels[] == before
    end

    f = [1.0, 10.0, 100.0]
    z = fill(measurement(1e-14, 4e-15) + im * measurement(1e-14, 0.0), 1, 1, 3)
    data = LineParameters(z, z, f)
    page = LineCableModels.plot(data; ydata=(R,), backend=:cairo, display_plot=false,
        controls=true, open_export=false, length_unit=:base, clip=false)
    axis = only(page.axes)
    Makie.colorbuffer(page.figure)
    with_spread = axis.finallimits[].widths[2]
    errorbars = filter(plot -> plot isa Makie.Errorbars, axis.scene.plots)
    @test !isempty(errorbars)
    foreach(plot -> plot.visible[] = false, errorbars)
    Makie.colorbuffer(page.figure)
    @test 0 < axis.finallimits[].widths[2] < with_spread
    # A caller's additional native series participates in fitting, unless the
    # caller excludes it through the native per-dimension autolimits attribute.
    added = lines!(axis, f, fill(1e-12, 3); yautolimits=true)
    page.controls[:reset].clicks[] += 1
    Makie.colorbuffer(page.figure)
    @test sum((axis.finallimits[].origin[2], axis.finallimits[].widths[2])) >= 1e-12
    added.yautolimits[] = false
    page.controls[:reset].clicks[] += 1
    Makie.colorbuffer(page.figure)
    @test axis.finallimits[].widths[2] < with_spread
    @test isequal(Z(data), z)
end

@testitem "Makie addons / scale controls and visibility keep multipliers synchronized" tags=[:visual] begin
    using CairoMakie
    frequency = [1.0, 10.0, 100.0]
    small = LineParameters(fill(1e-14 + 1.0im, 1, 1, 3), fill(1.0im, 1, 1, 3), frequency)
    large = LineParameters(fill(1e8 + 1.0im, 1, 1, 3), fill(1.0im, 1, 1, 3), frequency)
    plot = Makie.plot((; small, large), R; backend=:cairo, display_plot=false,
        controls=true, length_unit=:base, quantity_units=:base, clip=false)
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([0.0, 1e8]) == ["0", "100"]
    for item in filter(item -> item isa Makie.Lines, axis.scene.plots)
        last(first(item[1][])) > 1.0 && (item.visible[] = false)
    end
    Makie.colorbuffer(plot.figure)
    @test axis.finallimits[].widths[2] < 1e-12
    @test occursin("−15", repr(axis.ylabel[]))
    for _ in 1:2
        plot.controls[:ylog].active[] = true
        Makie.colorbuffer(plot.figure)
        @test axis.yscale[] === log10
        # A constant's modest log view is value-labelled, with the same single
        # engineering multiplier as other narrow numeric views.
        @test length(findall("× 10",repr(axis.ylabel[]))) == 1
        @test all(label -> label isa AbstractString, axis.yaxis.ticklabels[])
        plot.controls[:ylog].active[] = false
        Makie.colorbuffer(plot.figure)
        @test axis.yscale[] === identity
        @test axis.ytickformat[]([1e-14]) == ["10"]
        @test length(findall("× 10", repr(axis.ylabel[]))) == 1
    end
end

@testitem "Makie addons / scientific formatting respects native customization" tags=[:visual] begin
    using CairoMakie
    custom = values -> fill("custom", length(values))
    plot = LineCableModels.plotwindow(; title="Native customization", backend=:cairo,
        display_plot=false, controls=false) do layout
        for (column, attributes) in enumerate((
            (; ytickformat=custom),
            (; yticks=([0.0, 1e8], ["low", "high"])),
            (;),
        ))
            axis = Axis(layout[1, column]; ylabel="Response", attributes...)
            lines!(axis, [1.0, 10.0], [1e7, 1e8])
        end
    end
    first_axis, labelled_axis, axis = plot.axes
    Makie.colorbuffer(plot.figure)
    @test first_axis.ytickformat[] === custom
    @test first_axis.ylabel[] == "Response"
    @test labelled_axis.yticks[] == ([0.0, 1e8], ["low", "high"])
    @test labelled_axis.ylabel[] == "Response"
    @test occursin("× 10", repr(axis.ylabel[]))
    axis.ytickformat[] = custom
    @test axis.ylabel[] == "Response"
    ylims!(axis, 1e6, 1e9)
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[] === custom
    axis.ytickformat[] = Makie.automatic
    @test axis.ytickformat[]([0.0, 1e9]) == ["0", "1"]
    @test occursin("× 10", repr(axis.ylabel[]))
    for scale in (log10, sqrt)
        axis.yscale[] = scale
        Makie.colorbuffer(plot.figure)
        @test !isempty(axis.yaxis.ticklabels[])
        @test axis.ylabel[] == "Response"
    end
    axis.yscale[] = identity
    Makie.colorbuffer(plot.figure)
    @test axis.ytickformat[]([0.0, 1e9]) == ["0", "1"]
    @test length(findall("× 10", repr(axis.ylabel[]))) == 1
end

@testitem "Makie addons / previews and statistical plots share live multipliers" tags=[:visual] setup=[TestFixtures] begin
    using CairoMakie
    copper = Material(kind=:conductor, rho=1.7241e-8)
    design = @cable "axis-preview" begin
        @terminal :core begin
            core(copper; r=1e-3)
        end
    end
    result = TestFixtures.cable_monte_carlo_result()
    plots = (
        preview(design; backend=:cairo, display_plot=false, controls=false),
        Makie.hist(result, R; backend=:cairo, display_plot=false, controls=false),
        Makie.ecdfplot(result, R; backend=:cairo, display_plot=false, controls=false),
    )
    for plot in plots
        axis = only(plot.axes)
        xlims!(axis, -2e-12, 2e-12)
        ylims!(axis, -2e-12, 2e-12)
        Makie.colorbuffer(plot.figure)
        @test axis.xtickformat[]([-1e-12, 0.0, 1e-12]) == ["-1", "0", "1"]
        @test axis.ytickformat[]([-1e-12, 0.0, 1e-12]) == ["-1", "0", "1"]
        @test occursin("−12", repr(axis.xlabel[]))
        @test occursin("−12", repr(axis.ylabel[]))
    end
end

@testitem "Makie addons / native plot windows share axis formatting" tags=[:visual] begin
    using CairoMakie

    plot = LineCableModels.plotwindow(; title="Axis formatting", backend=:cairo,
        display_plot=false, controls=false) do layout
        axis = Axis(layout[1, 1]; xlabel="Time [s]", ylabel="Response")
        lines!(axis, [0.0, 1e-12, 2e-12], [-2e8, 0.0, 2e8])
    end
    axis = only(plot.axes)
    Makie.colorbuffer(plot.figure)
    @test axis.xtickformat[] isa Function
    @test axis.ytickformat[] isa Function
    if axis.xtickformat[] isa Function && axis.ytickformat[] isa Function
        @test axis.xtickformat[]([0.0, 1e-12, 2e-12]) == ["0", "1", "2"]
        @test axis.ytickformat[]([-2e8, 0.0, 2e8]) == ["-200", "0", "200"]
    end
end
