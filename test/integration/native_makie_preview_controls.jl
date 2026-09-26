@testitem "Makie addons / preview grouping and late horizontal legend retain geometry" tags=[:visual] begin
    using CairoMakie

    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    design = build(CableDesign, "preview-grouping",
        terminal(:core, core(copper; r=1e-3), insulation(dielectric; t=0.2e-3)))
    options = (; backend=:cairo, display_plot=false, controls=false, open_export=false)
    grouped = preview(design; options..., display_colorbars=false,
        legend_group=Dict(:core=>:metal, :insulation=>:dielectric),
        legend_labels=group -> uppercasefirst(String(group)))
    @test Set(values(grouped.addon_state.labels)) == Set(("Metal", "Dielectric"))
    @test isempty(grouped.colorbars)
    @test length(grouped.axes) == 1
    @test_throws ArgumentError preview(design; options..., legend_group=1)
    @test_throws ArgumentError preview(design; options..., legend_labels=1)

    collection = preview([design, design, design]; options..., layout=nothing)
    @test length(collection.axes) == 3
    @test length(collection.colorbars) == 3
    docked = preview(design; options..., display_legend=false, colorbar_position=:top)
    @test docked.legend === nothing
    geometry_plots = [copy(axis.scene.plots) for axis in docked.axes]
    figurelegend!(docked; position=:top, orientation=:horizontal, max_fraction=0.5)
    Makie.colorbuffer(docked.figure)
    @test docked.legend.orientation[] == :horizontal
    bounds = docked.legend.layoutobservables.computedbbox[]
    scale_bounds = [bar.layoutobservables.computedbbox[] for bar in docked.colorbars]
    @test bounds.origin[1] + bounds.widths[1] <=
        minimum(box.origin[1] for box in scale_bounds) + 1
    figurelegend!(docked; position=:right, max_fraction=0.5)
    figurelegend!(docked; position=:top, orientation=:horizontal, max_fraction=0.5)
    Makie.colorbuffer(docked.figure)
    @test [axis.scene.plots for axis in docked.axes] == geometry_plots
    @test length(docked.colorbars) == 3
    @test all(isfinite, docked.legend.layoutobservables.computedbbox[].widths)

    system = build(LineCableSystem, design, (0.0, -0.1); connections=Dict(:core=>1))
    system_plot = preview(system; options..., display_colorbars=false)
    @test isempty(system_plot.colorbars)
    for zoom_factor in ("invalid", -1.0, Inf)
        @test_throws ArgumentError preview(system; options..., zoom_factor)
    end
end

@testitem "Makie addons / constant positive observations retain logarithmic limits" tags=[:visual] begin
    using CairoMakie

    frequency = [10.0, 100.0, 1000.0]
    impedance = fill(ComplexF64(1e-4, 1e-3), 1, 1, 3)
    admittance = fill(ComplexF64(0.0, 1e-8), 1, 1, 3)
    parameters = LineParameters(copy(impedance), copy(admittance), frequency)
    options = (; backend=:cairo, display_plot=false, controls=false,
        length_unit=:base, quantity_units=:base)
    plot = Makie.plot(parameters, (R, 1, 1, :); options..., yscale=:log10)
    Makie.colorbuffer(plot.figure)
    axis = only(plot.axes)
    @test axis.yscale[] === Makie.log10
    limits = axis.finallimits[]
    @test 0 < limits.origin[2] < 1e-4 < limits.origin[2] + limits.widths[2]
    @test all(isfinite, limits.widths)
    line = only(filter(item -> item isa Makie.Lines, axis.scene.plots))
    @test last.(line[1][]) ≈ fill(1e-4, 3)
    @test parameters.Z.values == impedance
    @test parameters.Y.values == admittance
end

@testitem "Makie / callable controls own subscriptions and preserve retained figures" tags=[:visual] begin
    using CairoMakie
    calls=Ref(0)
    constructions=Ref(0)
    controls=Any[]
    function custom(p)
        constructions[]+=1
        push!(controls,addwidget!((p,cell) -> Button(cell;label="Count"),p,:count;
            event=button -> button.clicks,callback=(p,_) -> (calls[]+=1),success="Counted"))
    end
    raw=LineParameters(fill(1+2im,1,1,3),fill(1+2im,1,1,3),[1.,10.,100.])
    options=(backend=:cairo,display_plot=false,open_export=false)
    pages=LineCableModels.plot(raw;ydata=(R,L),widgets=(custom,),options...)
    @test constructions[]==2 && calls[]==0
    p=first(pages)
    @test p.status===p.addon_state.shell.status
    controls[1].clicks[]+=1
    @test calls[]==1 && p.status[]=="Counted"
    for _ in 1:3
        resetview!(p)
        axisscale!(p,:y,:log10)
        axisscale!(p,:y,:linear)
    end
    controls[1].clicks[]+=1
    @test calls[]==2 && constructions[]==2
    before=length(p.figure.content)
    @test_throws ErrorException addwidget!(p,:broken) do p,slot
        Label(slot,"partial")
        error("builder failed")
    end
    @test length(p.figure.content)==before
    @test !haskey(p.controls,:broken)
    @test_throws ArgumentError addwidget!((p,slot) -> Button(slot),p,:reset)
    @test_throws ArgumentError addwidget!((p,slot) -> Button(slot),p,:count)
    @test_throws ArgumentError addwidget!((p,slot) -> Button(slot),p,:unpaired;event=b -> b.clicks)
    @test removewidget!(p,:count)===p
    controls[1].clicks[]+=1
    @test calls[]==2
    @test !haskey(p.controls,:count)
    @test_throws ArgumentError removewidget!(p,:reset)
    @test_throws ArgumentError removewidget!(p,:unknown)
    composite_button=Ref{Any}(nothing)
    group=addwidget!(p,:composite;event=g -> composite_button[].clicks,callback=(p,_) -> (calls[]+=1)) do p,slot
        grid=GridLayout(slot)
        composite_button[]=Button(grid[1,1];label="Composite")
        Label(grid[1,2],"companion")
        grid
    end
    @test group isa GridLayout
    composite_button[].clicks[]+=1
    @test calls[]==3
    removewidget!(p,:composite)
    composite_button[].clicks[]+=1
    @test calls[]==3
    none=LineCableModels.plot(raw;ydata=(R,),controls=false,widgets=(custom,),options...)
    @test constructions[]==2 && isempty(none.controls)
    @test_throws ArgumentError addwidget!((p,slot) -> Button(slot),none,:late)
    @test axisscale!(none,:y,:log10)===none
    @test only(none.axes).yscale[]===log10
    @test resetview!(none;x=false)===none
    @test_throws ArgumentError resetview!(none;panel=(8,8))
    @test_throws ArgumentError axisscale!(none,:z,:linear)
    canvas=LineCableModels.plotwindow(;title="Caller canvas",widgets=(custom,),options...) do grid
        nested=GridLayout(grid[1,1])
        axis=Axis(nested[2,3])
        lines!(axis,[-2.,0.,3.],[-1.,0.,2.])
    end
    @test constructions[]==3
    @test axisscale!(canvas,:x,:log10;panel=1)===canvas
    @test only(canvas.axes).xscale[]!==log10
    @test_throws DomainError axisscale!(canvas,:x,log10;panel=1)
    @test resetview!(canvas;panel=1)===canvas
end

@testitem "Makie / composite widget destruction and SVG share native lifetime" tags=[:visual] begin
    using CairoMakie
    p=LineCableModels.plotwindow(;title="Composite controls",backend=:cairo,display_plot=false) do grid
        axis=Axis(grid[1,1]);lines!(axis,[1.,2.,3.],[2.,4.,8.])
    end
    calls=Ref(0)
    button=Ref{Any}()
    widget=addwidget!(p,:composite;event=grid -> button[].clicks,
        callback=(p,_) -> (calls[]+=1)) do p,slot
        grid=GridLayout(slot)
        button[]=Button(grid[1,1];label="Custom action",height=60)
        Label(grid[1,2],"A companion label")
        grid
    end
    button[].clicks[]+=1
    @test calls[]==1
    before_size=Tuple(p.figure.scene.viewport[].widths)
    before_view=only(p.axes).targetlimits[]
    visible=button[].blockscene.visible[]
    states=Bool[]
    observer=on(p.figure.scene,Makie.events(p.figure).tick) do tick
        tick.state===Makie.OneTimeRenderTick && push!(states,button[].blockscene.visible[])
    end
    mktempdir() do directory
        for theme in (:default,:publication)
            export_svg(p;path=joinpath(directory,"custom-$theme.svg"),theme,open_file=false)
            @test calls[]==1
            @test button[].blockscene.visible[]==visible
            @test Tuple(p.figure.scene.viewport[].widths)==before_size
            @test only(p.axes).targetlimits[]==before_view
        end
    end
    off(observer)
    @test !isempty(states) && all(!,states)
    clicks=button[].clicks
    delete!(button[])
    clicks[]+=1
    @test calls[]==1
    @test removewidget!(p,:composite)===p
    @test !haskey(p.controls,:composite)
end
