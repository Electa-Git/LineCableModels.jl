module GoldenFixtures

using LineCableModels
using Makie

export custom_layout_plot

function custom_layout_plot(;
        backend = :cairo,
        display_plot::Bool = false,
        controls::Bool = true
)
    return LineCableModels.plotwindow(;
        title = "Nested native layout",
        size = (900, 650),
        backend = backend,
        display_plot = display_plot,
        controls = controls,
        open_export = false,
        export_name = "nested_dashboard"
    ) do canvas
        plots = GridLayout()
        plots.default_rowgap = Fixed(6)
        plots.default_colgap = Fixed(6)
        canvas[1, 1][] = plots

        line_axis = Axis(
            plots[1, 1:2];
            title = "Spanning response",
            xlabel = "x",
            ylabel = "y"
        )
        x = [1.0, 2.0, 3.0]
        response_a = [1.0, 2.0, 1.5]
        response_b = [1.4, 1.6, 2.2]
        line_a = lines!(
            line_axis,
            x,
            response_a;
            label = "response A",
            linewidth = 2,
            color = :steelblue
        )
        line_b = lines!(
            line_axis,
            x,
            response_b;
            label = "response B",
            linewidth = 2,
            color = :darkorange
        )
        axislegend(line_axis; position = :rt)

        scatter_axis = Axis(
            plots[2, 1];
            title = "Samples",
            xlabel = "x",
            ylabel = "y"
        )
        samples = [0.8, 1.7, 2.4]
        scatter_plot = scatter!(
            scatter_axis,
            x,
            samples;
            label = "samples",
            color = :seagreen,
            markersize = 10
        )
        axislegend(scatter_axis; position = :rt)

        heatmap_axis = Axis(
            plots[2, 2];
            title = "Field",
            xlabel = "x",
            ylabel = "y"
        )
        field_x = [1.0, 2.0]
        field_y = [1.0, 2.0]
        field_plot = heatmap!(
            heatmap_axis,
            field_x,
            field_y,
            [0.0 0.5; 0.75 1.0];
            colormap = :viridis
        )
        Colorbar(
            plots[2, 3],
            field_plot;
            label = "field value",
            ticks = ([0.0, 0.5, 1.0], ["0", "0.5", "1"])
        )
        return nothing
    end
end

end
