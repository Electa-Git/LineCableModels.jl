@testitem "PlotBuilder / Cairo / responsive outer bounds" tags=[:visual] begin
    if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false") != "true"
        error("set LINECABLEMODELS_TEST_PLOTTING=true to run the visual contract")
    else
        using CairoMakie

        ui_components = Base.get_extension(
            LineCableModels,
            :LineCableModelsMakieExt
        ).UIComponents
        @test ui_components._window_padding((0, 1, 2, 4)) == (3.0, 3.0, 3.0, 4)

        function block_bounds(block)
            box = block.layoutobservables.computedbbox[]
            protrusions = block.layoutobservables.protrusions[]
            return (;
                left = box.origin[1] - protrusions.left,
                right = box.origin[1] + box.widths[1] + protrusions.right,
                bottom = box.origin[2] - protrusions.bottom,
                top = box.origin[2] + box.widths[2] + protrusions.top
            )
        end

        function assert_contained(handle)
            Makie.update_state_before_display!(handle.figure)
            width, height = handle.page.size
            bounds = block_bounds.(filter(
                block -> hasproperty(block, :layoutobservables),
                handle.figure.content
            ))
            @test all(bound -> bound.left >= 2, bounds)
            @test all(bound -> bound.right <= width - 2, bounds)
            @test all(bound -> bound.bottom >= 2, bounds)
            @test all(bound -> bound.top <= height - 2, bounds)
            return nothing
        end

        frequency = 10.0 .^ range(1, 6; length = 21)
        omega = reshape(2π .* frequency, 1, 1, :)
        resistance = [1.0 0.18 0.12; 0.18 1.1 0.16; 0.12 0.16 1.2] .* 1.0e-4
        inductance = [2.0 0.45 0.35; 0.45 2.1 0.4; 0.35 0.4 2.2] .* 1.0e-7

        function fixture(scale)
            impedance = repeat(scale .* resistance, 1, 1, length(frequency)) .+
                        im .* repeat(scale .* inductance, 1, 1,
                length(frequency)) .* omega
            return LineParameters(impedance, impedance .* 1.0e-5, frequency)
        end

        comparison = only(Makie.plot(
            fixture(1.0),
            fixture(1.015),
            fixture(0.985);
            legend = (
                "Reference · Wedepohl",
                "LineCableModels · Pollaczek",
                "LineCableModels · direct numerical integration"
            ),
            requests = (R,),
            xscale = :log10,
            fig_size = (1400, 900),
            backend = :cairo,
            display_plot = false
        ))
        assert_contained(comparison)
        legend = comparison.controls[:legend]
        @test block_bounds(first(comparison.panels).axis).left >= 19
        @test block_bounds(legend).right <= comparison.page.size[1] - 19
        comparison_payload = comparison.page.payload
        @test first(comparison_payload.titles) ==
              "R[1,1] · Series resistance"

        conventional = Makie.plot(
            fixture(1.0),
            (R, L);
            backend = :cairo,
            display_plot = false
        )
        foreach(assert_contained, conventional)
    end
end
