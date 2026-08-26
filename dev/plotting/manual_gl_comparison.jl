using Test
using LineCableModels
using GLMakie

set_backend!(:gl)

frequency=10.0 .^ range(0, 6; length = 61)
omega=reshape(2π .* frequency, 1, 1, :)
resistance=[1.0 0.18 0.12; 0.18 1.1 0.16; 0.12 0.16 1.2] .* 1.0e-4
inductance=[2.0 0.45 0.35; 0.45 2.1 0.40; 0.35 0.40 2.2] .* 1.0e-7
conductance=[3.0 0.6 0.4; 0.6 3.2 0.5; 0.4 0.5 3.4] .* 1.0e-9
capacitance=[4.0 0.9 0.7; 0.9 4.2 0.8; 0.7 0.8 4.4] .* 1.0e-10

function comparison_fixture(scale, phase_shift)
    impedance=repeat(scale .* resistance, 1, 1, length(frequency)) .+
              im .* repeat(scale .* inductance, 1, 1, length(frequency)) .* omega
    admittance=repeat(scale .* conductance, 1, 1, length(frequency)) .+
               im .* repeat((scale + phase_shift) .* capacitance, 1, 1,
        length(frequency)) .* omega
    return LineParameters(impedance, admittance, frequency)
end

reference=comparison_fixture(1.0, 0.0)
pollaczek=comparison_fixture(1.015, 0.005)
numerical=comparison_fixture(0.985, -0.004)

labels=(
    "Reference · Wedepohl",
    "LineCableModels · Pollaczek",
    "LineCableModels · direct numerical integration"
)

plots=Makie.plot(
    reference,
    pollaczek,
    numerical;
    legend = labels,
    requests = (R, L, G, C),
    xscale = :log10,
    fig_size = (1400, 900),
    backend = :gl,
    display_plot = true,
    open_export = false
)

@testset "manual GL line-parameter comparison" begin
    @test length(plots) == 4
    @test all(plot -> length(plot.panels) == 9, plots)
    @test all(
        panel -> length(panel.metadata.series) == 3,
        (panel for plot in plots for panel in plot.panels)
    )
    @test all(plot -> plot.context.window !== nothing, plots)
    @test all(plot -> plot.context.status[] == "Ready.", plots)
    @test all(plot -> haskey(plot.controls, :legend), plots)
    @test all(
        plot -> sort!(collect(keys(plot.controls))) ==
                [:export_svg, :legend, :reset, :xlog, :ylog],
        plots
    )
end

println("Opened one 3×3 comparison grid for each of R, L, G, and C.")
println("Resize the windows to inspect the matrix grid and responsive legend.")
println("Toggle a legend entry and confirm that it hides the same result in every panel.")
println("Close all windows to finish, or press Ctrl+C here.")

try
    while any(plot -> isopen(plot.context.window), plots)
        sleep(0.1)
    end
finally
    GLMakie.closeall()
end
