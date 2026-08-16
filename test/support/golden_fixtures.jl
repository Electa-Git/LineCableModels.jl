module GoldenFixtures

using LineCableModels

export custom_layout_render_spec

function custom_layout_render_spec()
    PB = LineCableModels.PlotBuilder
    UH = LineCableModels.UnitHandler

    xaxis = PB.AxisSpec(
        :x,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "x"
    )
    yaxis = PB.AxisSpec(
        :y,
        UH.QuantityTag{:dimensionless}(),
        UH.Units(),
        "y"
    )
    line_view = PB.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        "Spanning response",
        [
            PB.SeriesSpec(
                :line,
                [1.0, 2.0, 3.0],
                [1.0, 2.0, 1.5],
                nothing,
                "response A";
                group = :response_a,
                attributes = (; linewidth = 2, color = :steelblue)
            ),
            PB.SeriesSpec(
                :line,
                [1.0, 2.0, 3.0],
                [1.4, 1.6, 2.2],
                nothing,
                "response B";
                group = :response_b,
                attributes = (; linewidth = 2, color = :darkorange)
            )
        ],
        (; panel = :response);
        placement = PB.PlacementSpec(:canvas, PB.GridArea(1, 1:2))
    )
    scatter_view = PB.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        "Samples",
        [PB.SeriesSpec(
            :scatter,
            [1.0, 2.0, 3.0],
            [0.8, 1.7, 2.4],
            nothing,
            "samples";
            group = :samples,
            attributes = (; color = :seagreen, markersize = 10)
        )],
        (; panel = :samples);
        placement = PB.PlacementSpec(:canvas, PB.GridArea(2, 1))
    )
    heatmap_view = PB.ViewSpec(
        xaxis,
        yaxis,
        nothing,
        "Field",
        [PB.SeriesSpec(
            :heatmap,
            [1.0, 2.0],
            [1.0, 2.0],
            [0.0 0.5; 0.75 1.0],
            nothing;
            attributes = (; colormap = :viridis)
        )],
        (; panel = :field);
        placement = PB.PlacementSpec(:canvas, PB.GridArea(2, 2))
    )

    layout = PB.LayoutSpec(
        :nested_dashboard,
        [
            PB.GridSpec(
                :root;
                rows = PB.AbstractTrackSize[
                    PB.FixedTrack(36), PB.RelativeTrack(), PB.FixedTrack(20)],
                columns = PB.AbstractTrackSize[PB.ContentTrack(), PB.ContentTrack()],
                rowgap = 6,
                columngap = 12,
                padding = (20, 20, 28, 28)
            ),
            PB.GridSpec(
                :plots;
                parent = :root,
                area = PB.GridArea(2, 1),
                rows = PB.AbstractTrackSize[PB.RelativeTrack()],
                columns = PB.AbstractTrackSize[PB.RelativeTrack()]
            ),
            PB.GridSpec(
                :side;
                parent = :root,
                area = PB.GridArea(2, 2),
                rows = PB.AbstractTrackSize[PB.RelativeTrack(), PB.ContentTrack()],
                columns = PB.AbstractTrackSize[PB.ContentTrack()],
                rowgap = 4
            )
        ],
        [
            PB.SlotSpec(:toolbar, :root, PB.GridArea(1, 1:2); halign = :left),
            PB.SlotSpec(:canvas, :plots, PB.GridArea(1, 1)),
            PB.SlotSpec(:status, :root, PB.GridArea(3, 1:2); halign = :left),
            PB.SlotSpec(:legend, :side, PB.GridArea(1, 1); halign = :left, valign = :top),
            PB.SlotSpec(
                :colorbars, :side, PB.GridArea(2, 1); halign = :left, valign = :top)
        ]
    )
    page = PB.PageSpec(
        "Nested PlotBuilder layout",
        (900, 650),
        (; kind = :nested_dashboard),
        layout,
        [line_view, scatter_view, heatmap_view];
        colorbars = [PB.ColorbarSpec(
            "field value",
            :viridis,
            (0.0, 1.0),
            ([0.0, 0.5, 1.0], ["0", "0.5", "1"])
        )],
        export_spec = PB.ExportSpec(name = "nested_dashboard", open_file = false)
    )
    return PB.RenderSpec(LineCableModels.DataModel.MaterialScalePlotSpec, [page])
end

end
