module GoldenFixtures

using LineCableModels
using Makie

export custom_layout_plot

function custom_layout_plot(;
        backend = :cairo,
        display_plot::Bool = false,
        controls::Bool = true,
        export_mode::Bool = false
)
    PB = LineCableModels.PlotBuilder
    return PB.plotwindow(;
        title = "Nested PlotBuilder layout",
        size = (900, 650),
        backend = backend,
        display_plot = display_plot,
        controls = controls,
        export_mode = export_mode,
        open_export = false,
        export_name = "nested_dashboard",
        colorbars = ((;
            label = "field value",
            colormap = :viridis,
            limits = (0.0, 1.0),
            ticks = ([0.0, 0.5, 1.0], ["0", "0.5", "1"])
        ),)
    ) do ui
        plots = GridLayout()
        plots.default_rowgap = Fixed(6)
        plots.default_colgap = Fixed(6)
        ui.canvas[1, 1][] = plots

        line_axis = PB.axis!(
            ui,
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
        PB.register!(
            ui,
            line_axis;
            xmetadata = (
                label = "x", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            ymetadata = (
                label = "y", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            groups = (
                response_a = (line_a,),
                response_b = (line_b,)
            ),
            labels = (
                response_a = "response A",
                response_b = "response B"
            ),
            data = (
                (; xdata = x, ydata = response_a, group = :response_a),
                (; xdata = x, ydata = response_b, group = :response_b)
            )
        )

        scatter_axis = PB.axis!(
            ui,
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
        PB.register!(
            ui,
            scatter_axis;
            xmetadata = (
                label = "x", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            ymetadata = (
                label = "y", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            groups = (samples = (scatter_plot,),),
            labels = (samples = "samples",),
            data = ((; xdata = x, ydata = samples, group = :samples),)
        )

        heatmap_axis = PB.axis!(
            ui,
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
        PB.register!(
            ui,
            heatmap_axis;
            xmetadata = (
                label = "x", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            ymetadata = (
                label = "y", scale = :linear,
                allowed_scales = (:linear,), exponent = 0
            ),
            groups = (field = (field_plot,),),
            data = ((; xdata = field_x, ydata = field_y, group = :field),)
        )
        return nothing
    end
end

end
