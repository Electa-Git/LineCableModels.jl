using Test
using LineCableModels
using GLMakie

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

@testset "manual GL plotting gate / first SVG loads its renderer" begin
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing
    @test Base.get_extension(
        LineCableModels, :LineCableModelsMakieExt
    ).current_backend_symbol() === :gl
    limits_before = [axis.limits[] for axis in handle.axes]
    background_before = handle.figure.scene.backgroundcolor[]
    rows_before = copy(handle.figure.layout.rowsizes)
    mktempdir() do directory
        cd(directory) do
            handle.controls[:export_svg].clicks[] += 1
            path = joinpath(directory, only(readdir(directory)))
            @test filesize(path) > 100
            @test occursin("<svg", read(path, String))
            @test handle.addon_state.shell.status[] == "Saved SVG to $path"
            cp(path, joinpath(ARTIFACT_DIRECTORY, "first-gl-export.svg"); force = true)
        end
    end
    @test Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) !== nothing
    @test Makie.current_backend() === GLMakie
    @test [axis.limits[] for axis in handle.axes] == limits_before
    @test handle.figure.scene.backgroundcolor[] == background_before
    @test handle.figure.layout.rowsizes == rows_before
end

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
    susceptance_axis = last(susceptance.axes)
    @test susceptance_axis.yscale[] === Makie.log10
    @test susceptance_axis.ytickformat[] === Makie.automatic
    @test susceptance_axis.ylabel[] == "Shunt susceptance [S/km]"
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
    @test length(tick_values) in 1:4
    @test all(isinteger, log10.(tick_values))
    @test all(isone, round.(diff(log10.(tick_values)); digits = 8))
    @test all(label -> label isa Makie.RichText, tick_labels)
end

GLMakie.closeall()
