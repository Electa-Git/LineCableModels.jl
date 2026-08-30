function _show_summary(io::IO, name::AbstractString, fields::Pair...)
    print(io, name, "(")
    for (index, (key, value)) in enumerate(fields)
        index == 1 || print(io, ", ")
        print(io, key, "=")
        show(io, value)
    end
    print(io, ")")
end

function Base.show(io::IO, value::LegendDefinition)
    return _show_summary(
        io,
        "LegendDefinition",
        :enabled => value.enabled,
        :interactive => value.interactive,
        :overflow => value.overflow
    )
end
Base.summary(io::IO, ::LegendDefinition) = print(io, "Legend definition")

function Base.show(io::IO, value::ColorbarDefinition)
    return _show_summary(
        io,
        "ColorbarDefinition",
        :label => value.label,
        :limits => value.limits,
        :ticks => length(first(value.ticks))
    )
end
Base.summary(io::IO, ::ColorbarDefinition) = print(io, "Colorbar definition")

function Base.show(io::IO, value::ExportDefinition)
    return _show_summary(
        io,
        "ExportDefinition",
        :theme => value.theme,
        :name => value.name,
        :open_file => value.open_file
    )
end
Base.summary(io::IO, ::ExportDefinition) = print(io, "Export definition")

function Base.show(io::IO, value::PlotPage)
    return _show_summary(
        io,
        "PlotPage",
        :title => value.title,
        :size => value.size,
        :key => value.key,
        :payload => nameof(typeof(value.payload))
    )
end
Base.summary(io::IO, value::PlotPage) = print(io, "Plot page \"", value.title, "\"")

function Base.show(io::IO, value::PlotRecipe)
    return _show_summary(
        io,
        "PlotRecipe",
        :definition => nameof(value.definition),
        :pages => length(value.pages)
    )
end
Base.summary(io::IO, value::PlotRecipe) =
    print(io, "Plot recipe with ", length(value.pages), " pages")

function Base.show(io::IO, value::UIPlot)
    return _show_summary(
        io,
        "UIPlot",
        :title => value.page.title,
        :panels => length(value.context.panels),
        :backend => value.context.backend
    )
end
Base.summary(io::IO, value::UIPlot) = print(io, "Rendered plot \"", value.page.title, "\"")

const _CompactPlotBuilderObject = Union{
    LegendDefinition,
    ColorbarDefinition,
    ExportDefinition,
    PlotPage,
    PlotRecipe,
    UIPlot
}

Base.show(io::IO, ::MIME"text/plain", value::_CompactPlotBuilderObject) = show(io, value)
