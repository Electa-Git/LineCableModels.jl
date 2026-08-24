function _show_summary(io::IO, name::AbstractString, fields::Pair...)
    print(io, name, "(")
    for (index, (key, value)) in enumerate(fields)
        index == 1 || print(io, ", ")
        print(io, key, "=")
        show(io, value)
    end
    print(io, ")")
end

function _data_shape(data)
    data === nothing ? nothing :
    data isa AbstractArray ? size(data) : nameof(typeof(data))
end

function Base.show(io::IO, value::FixedTrack)
    _show_summary(io, "FixedTrack", :value => value.value)
end
function Base.show(io::IO, value::RelativeTrack)
    _show_summary(io, "RelativeTrack", :weight => value.weight)
end
Base.show(io::IO, ::ContentTrack) = print(io, "ContentTrack()")
function Base.show(io::IO, value::GridArea)
    _show_summary(io, "GridArea", :rows => value.rows, :columns => value.columns)
end
function Base.show(io::IO, value::GridSpec)
    _show_summary(
        io,
        "GridSpec",
        :name => value.name,
        :rows => length(value.rows),
        :columns => length(value.columns)
    )
end
function Base.show(io::IO, value::SlotSpec)
    _show_summary(
        io,
        "SlotSpec",
        :name => value.name,
        :parent => value.parent,
        :area => value.area,
        :alignment => (value.halign, value.valign)
    )
end
function Base.show(io::IO, value::LayoutSpec)
    _show_summary(
        io,
        "LayoutSpec",
        :name => value.name,
        :grids => length(value.grids),
        :slots => length(value.slots)
    )
end
function Base.show(io::IO, value::PlacementSpec)
    _show_summary(io, "PlacementSpec", :slot => value.slot, :area => value.area)
end
function Base.show(io::IO, value::ControlSpec)
    _show_summary(
        io,
        "ControlSpec",
        :reset => value.reset,
        :export_svg => value.export_svg,
        :slot => value.slot
    )
end
function Base.show(io::IO, value::LegendSpec)
    _show_summary(
        io,
        "LegendSpec",
        :enabled => value.enabled,
        :interactive => value.interactive,
        :slot => value.slot,
        :overflow => value.overflow
    )
end
function Base.show(io::IO, value::ColorbarSpec)
    _show_summary(
        io,
        "ColorbarSpec",
        :label => value.label,
        :limits => value.limits,
        :ticks => length(first(value.ticks)),
        :slot => value.slot
    )
end
function Base.show(io::IO, value::StatusSpec)
    _show_summary(io, "StatusSpec", :enabled => value.enabled, :slot => value.slot)
end
function Base.show(io::IO, value::ExportSpec)
    _show_summary(
        io,
        "ExportSpec",
        :theme => value.theme,
        :name => value.name,
        :open_file => value.open_file
    )
end
function Base.show(io::IO, value::AxisSpec)
    _show_summary(
        io,
        "AxisSpec",
        :dimension => value.dim,
        :label => value.label,
        :scale => value.scale,
        :allowed_scales => value.allowed_scales
    )
end
function Base.show(io::IO, value::SeriesSpec)
    _show_summary(
        io,
        "SeriesSpec",
        :kind => value.kind,
        :x => _data_shape(value.xdata),
        :y => _data_shape(value.ydata),
        :z => _data_shape(value.zdata),
        :label => value.label,
        :group => value.group,
        :visible => value.visible
    )
end
function Base.show(io::IO, value::ViewSpec)
    _show_summary(
        io,
        "ViewSpec",
        :title => value.title,
        :series => length(value.series),
        :slot => value.placement.slot
    )
end
function Base.show(io::IO, value::PageSpec)
    _show_summary(
        io,
        "PageSpec",
        :title => value.title,
        :size => value.size,
        :views => length(value.views),
        :layout => value.layout.name
    )
end
function Base.show(io::IO, value::PlotRecipe)
    _show_summary(
        io,
        "PlotRecipe",
        :spec => nameof(value.spec),
        :pages => length(value.figures),
        :object => nameof(typeof(value.object))
    )
end
function Base.show(io::IO, value::UIPlot)
    backend = hasproperty(value.context, :backend) ? getproperty(value.context, :backend) :
              :unknown
    return _show_summary(
        io,
        "UIPlot",
        :title => value.page.title,
        :panels => length(value.panels),
        :backend => backend
    )
end

const _CompactPlotBuilderObject = Union{
    PlotRecipe,
    FixedTrack,
    RelativeTrack,
    ContentTrack,
    GridArea,
    GridSpec,
    SlotSpec,
    LayoutSpec,
    PlacementSpec,
    ControlSpec,
    LegendSpec,
    ColorbarSpec,
    StatusSpec,
    ExportSpec,
    AxisSpec,
    SeriesSpec,
    ViewSpec,
    PageSpec,
    UIPlot
}

Base.show(io::IO, ::MIME"text/plain", value::_CompactPlotBuilderObject) = show(io, value)
