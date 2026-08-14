using Test
using LineCableModels
using GLMakie

const ARTIFACT_DIRECTORY = abspath(get(
    ENV,
    "LINECABLEMODELS_GL_ARTIFACTS",
    joinpath(pwd(), "manual-gl-artifacts")
))
mkpath(ARTIFACT_DIRECTORY)

set_backend!(:gl)

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

plots = Makie.plot(parameters; mode = :RLCG, backend = :gl, display_plot = true)
handle = first(plots)

@testset "manual GL plotting gate without CairoMakie" begin
    @test !LineCableModels.PlotBuilder.BackendHandler.backend_available(:cairo)
    handle.controls[:export_svg].clicks[] += 1
    @test occursin("load CairoMakie", handle.context.status[])
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing
    @test LineCableModels.PlotBuilder.BackendHandler.current_backend_symbol() === :gl
end

using CairoMakie
set_backend!(:gl)

@testset "manual GL plotting gate" begin
    @test plots isa Vector{UIPlot}
    @test length(plots) == 4
    @test all(plot -> plot.context.window !== nothing, plots)
    @test length(unique(objectid(plot.context.window) for plot in plots)) == length(plots)
    @test sort!(collect(keys(handle.controls))) ==
          [:export_svg, :legend, :reset, :xlog, :ylog]
    @test handle.context.backend === :gl
    @test handle.figure.scene.backgroundcolor[] == Makie.to_color(:grey90)
    @test occursin("\\ue5d5", sprint(show, handle.controls[:reset].label[]))
    @test occursin("\\ue161", sprint(show, handle.controls[:export_svg].label[]))

    handle.controls[:xlog].active[] = true
    handle.controls[:ylog].active[] = true
    @test only(handle.panels).axis.xscale[] === Makie.log10
    @test only(handle.panels).axis.yscale[] === Makie.log10

    legend = handle.controls[:legend]
    entry = first(last(first(legend.entrygroups[])))
    Makie.toggle_visibility!(entry)
    @test any(plot_object -> !plot_object.visible[], only(handle.panels).plots)
    Makie.toggle_visibility!(entry)

    handle.controls[:reset].clicks[] += 1
    @test handle.context.status[] == "Axis limits reset"

    GLMakie.save(joinpath(ARTIFACT_DIRECTORY, "gl-ui.png"), handle.figure)
    cd(ARTIFACT_DIRECTORY) do
        before = Set(readdir())
        handle.controls[:export_svg].clicks[] += 1
        after = Set(readdir())
        created = filter(name -> endswith(name, ".svg"), collect(setdiff(after, before)))
        @test length(created) == 1
        svg_path = joinpath(ARTIFACT_DIRECTORY, only(created))
        @test filesize(svg_path) > 100
        @test occursin("<svg", read(svg_path, String))
    end
    @test startswith(handle.context.status[], "Saved SVG to ")
    @test LineCableModels.PlotBuilder.BackendHandler.current_backend_symbol() === :gl

    susceptance = last(Makie.plot(
        parameters;
        mode = :ZY,
        coord = :cart,
        backend = :gl,
        display_plot = false
    ))
    susceptance.controls[:ylog].active[] = true
    susceptance_axis = only(susceptance.panels).axis
    @test susceptance_axis.yscale[] === Makie.log10
    @test susceptance_axis.ytickformat[] === Makie.automatic
    @test susceptance_axis.ylabel[] == "Capacitive susceptance [S/km]"
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
    @test length(tick_values) > 2
    @test length(unique(round.(diff(log10.(tick_values)); digits = 8))) > 1
    @test all(label -> label == "" || label isa Makie.RichText, tick_labels)
end

GLMakie.closeall()
