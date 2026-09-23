using Test
if "--backends-first" in ARGS
    using GLMakie
    using LineCableModels
else
    using LineCableModels
    using GLMakie
end

const ARTIFACT_DIRECTORY = abspath(get(
    ENV,
    "LINECABLEMODELS_GL_ARTIFACTS",
    joinpath(pwd(), "manual-gl-artifacts")
))
mkpath(ARTIFACT_DIRECTORY)

frequency = [50.0, 100.0, 500.0]
omega = reshape(2π .* frequency, 1, 1, :)
resistance_values = reshape([1.0, 0.2, 0.2, 2.0], 2, 2, 1) .*
                    ones(1, 1, length(frequency)) .* 1.0e-4
inductance_values = fill(2.0e-7, 2, 2, length(frequency))
conductance_values = fill(3.0e-9, 2, 2, length(frequency))
capacitance_values = fill(4.0e-10, 2, 2, length(frequency))
parameters = LineParameters(
    complex.(resistance_values, inductance_values .* omega),
    complex.(conductance_values, capacitance_values .* omega),
    frequency
)

plots = Makie.plot(
    parameters,
    (R, L, G, C);
    backend = :gl,
    display_plot = true,
    open_export = false
)
handle = first(plots)
empty_shell = LineCableModels.plotwindow(; title="Text-only shell", backend=:gl,
        display_plot=false, open_export=false) do canvas
    Label(canvas[1,1], "Native canvas without axes")
end

@testset "manual GL plotting gate / no unloaded SVG action" begin
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing
    @test !haskey(handle.controls, :export_svg)
    @test isempty(empty_shell.controls)
    @test sort!(collect(keys(handle.controls))) == [:reset, :xlog, :ylog]
    mktempdir() do directory
        path = joinpath(directory, "not-created", "unavailable.svg")
        @test_throws r"import CairoMakie" export_svg(handle; path, open_file=false)
        @test isempty(readdir(directory))
    end
end

# Loading a renderer does not retrofit existing toolbars. Explicitly select GL
# again because importing CairoMakie itself activates the upstream backend.
import CairoMakie
GLMakie.activate!()

@testset "manual GL plotting gate / late-loaded renderer preserves the live view" begin
    @test !haskey(handle.controls, :export_svg)
    @test Base.get_extension(
        LineCableModels, :LineCableModelsMakieExt
    ).current_backend_symbol() === :gl
    axis = first(handle.axes)
    area = axis.scene.viewport[]
    mouse = Tuple(area.origin + area.widths / 2)
    Makie.events(axis.scene).mouseposition[] = mouse
    Makie.process_interaction(last(axis.interactions[:scrollzoom]),
        Makie.ScrollEvent(0, 3), axis)
    Makie.process_interaction(last(axis.interactions[:dragpan]),
        Makie.MouseEvent(Makie.MouseEventTypes.rightdrag,
            1.0, Point2d(0), Point2f(mouse .+ (20, 10)),
            0.0, Point2d(0), Point2f(mouse)), axis)
    limits_before = [axis.limits[] for axis in handle.axes]
    targets_before = [axis.targetlimits[] for axis in handle.axes]
    views_before = [axis.finallimits[] for axis in handle.axes]
    background_before = handle.figure.scene.backgroundcolor[]
    rows_before = copy(handle.figure.layout.rowsizes)
    mktempdir() do directory
        cd(directory) do
            path = export_svg(handle; open_file=false)
            @test filesize(path) > 100
            @test occursin("<svg", read(path, String))
            cp(path, joinpath(ARTIFACT_DIRECTORY, "first-gl-export.svg"); force = true)
        end
    end
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) !== nothing
    @test Makie.current_backend() === GLMakie
    @test [axis.limits[] for axis in handle.axes] == limits_before
    @test [axis.targetlimits[] for axis in handle.axes] == targets_before
    @test [axis.finallimits[] for axis in handle.axes] == views_before
    @test handle.figure.scene.backgroundcolor[] == background_before
    @test handle.figure.layout.rowsizes == rows_before
    # An empty toolbar is still interactive chrome. Its status row must not
    # leak into publication output when Cairo was loaded after construction.
    status_label = only(filter(block -> block isa Label && block.text[] == "Ready",
        empty_shell.figure.content))
    rendered_status = Bool[]
    observer = on(empty_shell.figure.scene, Makie.events(empty_shell.figure).tick) do tick
        tick.state === Makie.OneTimeRenderTick && push!(rendered_status, status_label.blockscene.visible[])
    end
    mktempdir() do directory
        export_svg(empty_shell; path=joinpath(directory, "text-only.svg"), open_file=false)
    end
    off(observer)
    @test !isempty(rendered_status) && !any(rendered_status)
    @test status_label.blockscene.visible[]
end

plots = Makie.plot(parameters, (R, L, G, C); backend=:gl,
    display_plot=true, open_export=false)
handle = first(plots)

@testset "manual GL plotting gate" begin
    @test plots isa Vector{UIPlot}
    @test length(plots) == 4
    @test length(handle.axes) == 4
    @test sort!(collect(keys(handle.controls))) ==
          [:export_svg, :reset, :xlog, :ylog]
    @test handle.figure.scene.backgroundcolor[] == Makie.to_color(:grey90)
    @test occursin("\\ue5d5", sprint(show, handle.controls[:reset].label[]))
    @test occursin("\\ue161", sprint(show, handle.controls[:export_svg].label[]))

    handle.controls[:xlog].active[] = true
    handle.controls[:ylog].active[] = true
    @test all(axis -> axis.xscale[] === Makie.log10, handle.axes)
    @test all(axis -> axis.yscale[] === Makie.log10, handle.axes)

    handle.controls[:reset].clicks[] += 1

    GLMakie.save(joinpath(ARTIFACT_DIRECTORY, "gl-ui.png"), handle.figure)
    mktempdir() do directory
        cd(directory) do
            handle.controls[:export_svg].clicks[] += 1
            svg_path = joinpath(directory, only(readdir(directory)))
            @test filesize(svg_path) > 100
            @test occursin("<svg", read(svg_path, String))
            @test handle.addon_state.shell.status[] == "Saved SVG to $svg_path"
            cp(svg_path, joinpath(ARTIFACT_DIRECTORY, "repeated-gl-export.svg"); force = true)
        end
    end
    @test Base.get_extension(
        LineCableModels, :LineCableModelsMakieExt
    ).current_backend_symbol() === :gl

    susceptance = last(Makie.plot(
        parameters;
        backend = :gl,
        display_plot = false
    ))
    susceptance.controls[:ylog].active[] = true
    Makie.colorbuffer(susceptance.figure)
    susceptance_axis = last(susceptance.axes)
    @test susceptance_axis.yscale[] === Makie.log10
    @test occursin("Shunt susceptance [S/km]", string(susceptance_axis.ylabel[]))
    limits = susceptance_axis.finallimits[]
    ymin = limits.origin[2]
    ymax = ymin + limits.widths[2]
    tick_values, tick_labels = Makie.get_ticks(
        susceptance_axis.yticks[],
        susceptance_axis.yscale[],
        susceptance_axis.ytickformat[],
        ymin,
        ymax
    )
    # Short logarithmic spans use readable mantissas, not necessarily decade
    # powers. Check the rendered scale/labels, not a superseded formatter type.
    @test length(tick_values) == length(tick_labels) >= 2
    @test issorted(tick_values) && allunique(tick_values)
    @test all(value -> isfinite(value) && value > 0, tick_values)
    @test allunique(string.(tick_labels))
end

@testset "manual GL retained figure survives screen closure" begin
    clicks=Ref(0)
    control=addwidget!((p,slot) -> Button(slot;label="Count"),handle,:count;
        event=button -> button.clicks,callback=(p,_) -> (clicks[]+=1))
    screen=Makie.getscreen(handle.figure.scene)
    @test screen!==nothing
    close(screen)
    control.clicks[]+=1
    @test clicks[]==1
    figuretitle!(handle,"Retained after closure")
    mktempdir() do directory
        @test isfile(export_svg(handle;path=joinpath(directory,"closed.svg"),open_file=false))
    end
    display(handle.figure)
    control.clicks[]+=1
    @test clicks[]==2
    removewidget!(handle,:count)
end

@testset "manual GL independent guide arrangement and resize" begin
    include(joinpath(@__DIR__, "../../test/support/scenarios.jl"))
    # Previously scales defaulted to a right column. Previews now default to a
    # bottom strip with left-side property labels and 24-pixel plot clearance.
    # Explicit group layout still owns cells independently of bar orientation.
    p=preview(CurrentScenarios.coaxial_design();backend=:gl,display_plot=true,
        size=(1200,900),legend_position=:right,
        legend_attributes=(halign=:left,valign=:top),
        colorbar_attributes=(vertical=false,width=160,height=14),
        colorbar_group_attributes=(layout=(1,3),colgap=16,halign=:center),
        guide_spacing=(rowgap=14,colgap=16))
    guide=p.addon_state.guides[(:colorbars,nothing)]
    bars=copy(p.colorbars)
    axis=only(p.axes)
    limits!(axis,-.024,.032,-.028,.028)
    view=axis.targetlimits[]
    data=copy(axis.scene.plots)
    callbacks=length(guide.subscriptions)
    clicks=Ref(0)
    button=addwidget!((p,slot) -> Button(slot;label="Count"),p,:guide_count;
        event=button -> button.clicks,callback=(p,_) -> (clicks[]+=1))
    GLMakie.save(joinpath(ARTIFACT_DIRECTORY,"guide-bottom.png"),p.figure)
    @test size(guide.layout[])==(1,3)
    # Complete item bounds reserve spacing; bar frames share a baseline/edge
    # independently of the unequal property labels inside those items.
    frames=[bar.layoutobservables.computedbbox[] for bar in bars]
    @test maximum(frame.origin[2] for frame in frames)-minimum(frame.origin[2] for frame in frames)<=1
    @test maximum(frame.origin[2]+frame.widths[2] for frame in frames)-minimum(frame.origin[2]+frame.widths[2] for frame in frames)<=1
    figurecolorbars!(p;position=:right,
        group_attributes=(layout=(3,1),rowgap=12,valign=:bottom),vertical=false)
    figurelegend!(p;position=:right,valign=:top,guide_spacing=14)
    resize!(p.figure,1000,800)
    GLMakie.colorbuffer(p.figure)
    frames=[bar.layoutobservables.computedbbox[] for bar in bars]
    @test maximum(frame.origin[1] for frame in frames)-minimum(frame.origin[1] for frame in frames)<=1
    @test maximum(frame.origin[1]+frame.widths[1] for frame in frames)-minimum(frame.origin[1]+frame.widths[1] for frame in frames)<=1
    legend=p.legend.layoutobservables.computedbbox[]
    scales=guide.layout[].layoutobservables.computedbbox[]
    @test legend.origin[2]-(scales.origin[2]+scales.widths[2])>=13
    @test size(guide.layout[])==(3,1)
    @test all(a===b for (a,b) in zip(bars,p.colorbars))
    @test axis.targetlimits[]==view
    @test axis.scene.plots==data
    @test p.controls[:guide_count]===button
    @test length(guide.subscriptions)==callbacks
    button.clicks[]+=1
    @test clicks[]==1
    GLMakie.save(joinpath(ARTIFACT_DIRECTORY,"guide-right.png"),p.figure)
end

GLMakie.closeall()
