function _sanitize_filename(value::AbstractString)
    sanitized = lowercase(strip(value))
    sanitized = replace(sanitized, r"[^0-9a-z]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "linecablemodels_plot" : sanitized
end

function _normalized_path_parts(path::AbstractString)
    parts = collect(splitpath(normpath(realpath(path))))
    return Sys.iswindows() ? lowercase.(parts) : parts
end

function _path_within(path::AbstractString, root::AbstractString)
    path_parts = _normalized_path_parts(path)
    root_parts = _normalized_path_parts(root)
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _export_directory()
    current = abspath(pwd())
    package = abspath(pkgdir(PlotBuilder))
    _path_within(current, package) || return current
    fallback = joinpath(tempdir(), EXPORT_FALLBACK_DIRECTORY)
    mkpath(fallback)
    return fallback
end

function _available_path(definition::ExportDefinition)
    base = _sanitize_filename(definition.name)
    stamp = Dates.format(Dates.now(), EXPORT_TIMESTAMP_FORMAT)
    directory = _export_directory()
    candidate = joinpath(directory, "$(base)_$(stamp).svg")
    index = 2
    while ispath(candidate)
        candidate = joinpath(directory, "$(base)_$(stamp)_$(index).svg")
        index += 1
    end
    return candidate
end

function _open_command(path::AbstractString)
    if Sys.iswindows()
        return Cmd(["cmd", "/c", "start", "", path])
    elseif Sys.isapple()
        executable = Sys.which("open")
        return executable === nothing ? nothing : `$executable $path`
    end
    executable = Sys.which("xdg-open")
    executable !== nothing && return `$executable $path`
    executable = Sys.which("gio")
    return executable === nothing ? nothing : `$executable open $path`
end

function _open_export(path::AbstractString)
    command = _open_command(path)
    command === nothing && return false
    try
        process = run(pipeline(ignorestatus(command); stdout = devnull, stderr = devnull))
        return success(process)
    catch error
        @warn "could not open exported SVG with the system application" path exception = (
            error,
            catch_backtrace()
        )
        return false
    end
end

function _current_scale(scale)
    scale === Makie.log10 && return :log10
    scale === Makie.identity && return :linear
    throw(ArgumentError("SVG export supports linear and log10 axis scales"))
end

function _current_limits(axis)
    limits = axis.finallimits[]
    xlimits = (limits.origin[1], limits.origin[1] + limits.widths[1])
    ylimits = (limits.origin[2], limits.origin[2] + limits.widths[2])
    return xlimits, ylimits
end

_replay_page(plot::UIPlot, ::Type{<:PlotBuilder.AbstractPlotDefinition}) = plot.page

function _current_recipe(plot::UIPlot)
    page = _replay_page(plot, plot.render.definition)
    return PlotRecipe(plot.render.definition, (page,))
end

function _block_vertical_bounds(block)
    layout = block.layoutobservables
    bounding_box = layout.computedbbox[]
    protrusions = layout.protrusions[]
    bottom = bounding_box.origin[2] - protrusions.bottom
    top = bounding_box.origin[2] + bounding_box.widths[2] + protrusions.top
    return Float64(bottom), Float64(top)
end

function _export_dock_growth(figure, page::PlotPage)
    legends = filter(block -> block isa Legend, figure.content)
    isempty(legends) && return 0.0
    legend_bottom = minimum(first(_block_vertical_bounds(legend)) for legend in legends)
    all_colorbars = filter(block -> block isa Colorbar, figure.content)
    length(all_colorbars) == length(page.colorbars) || error(
        "rendered colorbars no longer match the declarative page",
    )
    required_bottom = 0.0
    if !isempty(all_colorbars)
        scale_top = mapreduce(
            block -> last(_block_vertical_bounds(block)),
            max,
            all_colorbars
        )
        required_bottom = scale_top + COLORBAR_ROW_GAP
    end
    return max(0.0, required_bottom - legend_bottom)
end

function _fit_export_content!(figure, page::PlotPage)
    fitted_size = Makie.resize_to_layout!(figure)
    target_size = Tuple(max.(page.size, ceil.(Int, fitted_size)))
    Makie.resize!(figure, target_size...)
    for _ in 1:4
        Makie.update_state_before_display!(figure)
        growth = _export_dock_growth(figure, page)
        growth <= LEGEND_HEIGHT_TOLERANCE && break
        target_size = (target_size[1], target_size[2] + ceil(Int, growth))
        Makie.resize!(figure, target_size...)
    end
    Makie.update_state_before_display!(figure)
    return target_size
end

function PlotBuilder.export_svg(
        plot::UIPlot;
        path::Union{Nothing, AbstractString} = nothing,
        theme::Union{Nothing, Symbol} = nothing,
        open_file::Union{Nothing, Bool} = nothing
)
    PlotBuilder.backend_available(:cairo) || throw(
        ArgumentError(
        "SVG export requires CairoMakie; load CairoMakie first with `using CairoMakie`",
    ),
    )
    definition = plot.context.export_state
    output = path === nothing ? _available_path(definition) : abspath(String(path))
    export_theme = theme === nothing ? definition.theme : theme
    should_open = open_file === nothing ? definition.open_file : open_file
    lowercase(splitext(output)[2]) == ".svg" || throw(
        ArgumentError("SVG export paths must use the .svg extension"),
    )
    ispath(output) && throw(ArgumentError("refusing to overwrite existing file: $output"))
    plot.context.status[] = "Exporting SVG..."
    PlotBuilder.with_backend(:cairo) do
        one_page = _current_recipe(plot)
        exported = build(
            one_page;
            backend = :cairo,
            display = false,
            controls = false,
            export_mode = true,
            export_theme
        )
        exported_plot = only(exported)
        _fit_export_content!(exported_plot.figure, exported_plot.page)
        Makie.save(output, exported_plot.figure)
    end
    opened = should_open && _open_export(output)
    message = if opened
        "Saved SVG to $output and opened it with the system application"
    elseif should_open
        "Saved SVG to $output; automatic opening was unavailable"
    else
        "Saved SVG to $output"
    end
    plot.context.status[] = message
    @info message
    return output
end
