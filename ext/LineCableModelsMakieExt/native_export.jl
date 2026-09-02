function _native_sanitize_filename(value::AbstractString)
    sanitized = lowercase(strip(value))
    sanitized = replace(sanitized, r"[^0-9a-z]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "linecablemodels_plot" : sanitized
end

function _native_path_within(path::AbstractString, root::AbstractString)
    path_parts = collect(splitpath(normpath(realpath(path))))
    root_parts = collect(splitpath(normpath(realpath(root))))
    Sys.iswindows() &&
        (path_parts = lowercase.(path_parts); root_parts = lowercase.(root_parts))
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _native_export_directory()
    current = abspath(pwd())
    package = abspath(pkgdir(LineCableModels))
    _native_path_within(current, package) || return current
    fallback = joinpath(tempdir(), "linecablemodels-exports")
    mkpath(fallback)
    return fallback
end

function _native_available_path(plot::LineCableModels.UIPlot)
    base = _native_sanitize_filename(plot.export_name)
    stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    directory = _native_export_directory()
    candidate = joinpath(directory, "$(base)_$(stamp).svg")
    index = 2
    while ispath(candidate)
        candidate = joinpath(directory, "$(base)_$(stamp)_$(index).svg")
        index += 1
    end
    return candidate
end

function _native_open_command(path::AbstractString)
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

function _native_open_export(path::AbstractString)
    command = _native_open_command(path)
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

function _native_observable_snapshot!(snapshot, object, names)
    for name in names
        hasproperty(object, name) || continue
        value = getproperty(object, name)
        value isa Observable || continue
        push!(snapshot, value => value[])
    end
    return snapshot
end

function _native_set_observable!(snapshot, observable, value)
    push!(snapshot, observable => observable[])
    observable[] = value
    return observable
end

function _native_hide_layout_content!(snapshot, content)
    if content isa GridLayout
        for entry in content.content
            _native_hide_layout_content!(snapshot, entry.content)
        end
    elseif hasproperty(content, :blockscene)
        _native_set_observable!(snapshot, content.blockscene.visible, false)
    end
    return snapshot
end

function _native_hide_interactive_chrome!(snapshot, plot)
    isempty(plot.controls) && return nothing
    root = plot.figure.layout
    row_sizes = copy(root.rowsizes)
    row_gap = root.default_rowgap
    for entry in root.content
        rows = entry.span.rows
        if rows == 1:1 || rows == 3:3
            _native_hide_layout_content!(snapshot, entry.content)
        end
    end
    rowsize!(root, 1, Fixed(0))
    rowsize!(root, 3, Fixed(0))
    root.default_rowgap = Fixed(0)
    return (; root, row_sizes, row_gap)
end

function _native_restore_interactive_chrome!(layout_snapshot)
    layout_snapshot === nothing && return nothing
    for (index, size) in pairs(layout_snapshot.row_sizes)
        rowsize!(layout_snapshot.root, index, size)
    end
    layout_snapshot.root.default_rowgap = layout_snapshot.row_gap
    return nothing
end

function _native_publication_snapshot!(snapshot, plot::LineCableModels.UIPlot, theme::Symbol)
    theme in (:default, :publication) || throw(ArgumentError(
        "theme must be :default or :publication",
    ))
    _native_observable_snapshot!(snapshot, plot.figure.scene, (:backgroundcolor,))
    plot.figure.scene.backgroundcolor[] = Makie.to_color(:white)
    layout_snapshot = _native_hide_interactive_chrome!(snapshot, plot)
    theme === :default && return layout_snapshot

    latex_fonts = Makie.theme_latexfonts().attributes[:fonts][]
    figure_fonts = plot.figure.scene.theme[:fonts]
    for role in (:regular, :italic, :bold)
        _native_set_observable!(snapshot, figure_fonts[role], latex_fonts[role][])
    end
    regular = :regular
    bold = :bold
    for block in plot.figure.content
        if block isa Axis
            _native_observable_snapshot!(snapshot, block,
                (
                    :titlefont,
                    :xlabelfont,
                    :ylabelfont,
                    :xticklabelfont,
                    :yticklabelfont
                ))
            block.titlefont[] = bold
            block.xlabelfont[] = regular
            block.ylabelfont[] = regular
            block.xticklabelfont[] = regular
            block.yticklabelfont[] = regular
        elseif block isa Legend
            _native_observable_snapshot!(snapshot, block, (:labelfont, :titlefont))
            block.labelfont[] = regular
            block.titlefont[] = bold
        elseif block isa Colorbar
            _native_observable_snapshot!(snapshot, block, (:labelfont, :ticklabelfont))
            block.labelfont[] = regular
            block.ticklabelfont[] = regular
        elseif block isa Label
            _native_observable_snapshot!(snapshot, block, (:font,))
            block.font[] = block === plot.title ? bold : regular
        end
    end
    return layout_snapshot
end

function _native_restore_snapshot!(snapshot)
    for (observable, value) in Iterators.reverse(snapshot)
        observable[] = value
    end
    return nothing
end

function _native_restore_backend!(backend)
    backend isa Module || return nothing
    isdefined(backend, :activate!) || return nothing
    Base.invokelatest(getproperty(backend, :activate!))
    return nothing
end

function LineCableModels.export_svg(
        plot::LineCableModels.UIPlot;
        path::Union{Nothing, AbstractString} = nothing,
        theme::Union{Nothing, Symbol} = nothing,
        open_file::Union{Nothing, Bool} = nothing
)
    Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt) === nothing &&
        throw(ArgumentError(
            "SVG export requires CairoMakie; load CairoMakie first with `using CairoMakie`",
        ))
    output = path === nothing ? _native_available_path(plot) : abspath(String(path))
    export_theme = theme === nothing ? plot.export_theme : theme
    should_open = open_file === nothing ? plot.open_export : open_file
    lowercase(splitext(output)[2]) == ".svg" || throw(ArgumentError(
        "SVG export paths must use the .svg extension",
    ))
    ispath(output) && throw(ArgumentError(
        "refusing to overwrite existing file: $output",
    ))
    mkpath(dirname(output))
    previous_backend = Makie.current_backend()
    snapshot = Pair{Any, Any}[]
    layout_snapshot = nothing
    try
        layout_snapshot = _native_publication_snapshot!(
            snapshot,
            plot,
            export_theme
        )
        _addon_activate_backend(:cairo)
        with_theme(_addon_theme(
            export_mode = true,
            export_theme = export_theme
        )) do
            Makie.update_state_before_display!(plot.figure)
            Makie.save(output, plot.figure)
        end
    finally
        _native_restore_interactive_chrome!(layout_snapshot)
        _native_restore_snapshot!(snapshot)
        _native_restore_backend!(previous_backend)
    end
    opened = should_open && _native_open_export(output)
    message = if opened
        "Saved SVG to $output and opened it with the system application"
    elseif should_open
        "Saved SVG to $output; automatic opening was unavailable"
    else
        "Saved SVG to $output"
    end
    @info message
    return output
end
