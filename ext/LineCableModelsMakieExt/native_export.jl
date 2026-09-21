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

function LineCableModels.export_svg(
        plot::LineCableModels.UIPlot;
        path::Union{Nothing, AbstractString} = nothing,
        theme::Union{Nothing, Symbol} = nothing,
        open_file::Union{Nothing, Bool} = nothing
)
    cairo = Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt)
    cairo === nothing && throw(ArgumentError(
        "SVG export requires CairoMakie to be loaded; run `import CairoMakie` first. " *
        "For interactive plots, select backend=:gl after importing both backends.",
    ))
    export_theme = theme === nothing ? plot.export_theme : theme
    export_theme in (:default,:publication) || throw(ArgumentError("theme must be :default or :publication"))
    output = path === nothing ? _native_available_path(plot) : abspath(String(path))
    should_open = open_file === nothing ? plot.open_export : open_file
    lowercase(splitext(output)[2]) == ".svg" || throw(ArgumentError(
        "SVG export paths must use the .svg extension",
    ))
    ispath(output) && throw(ArgumentError(
        "refusing to overwrite existing file: $output",
    ))
    mkpath(dirname(output))
    _addon_export_presentation!(plot,export_theme) do
        # Preserve the current figure and view. Saving must not run the native
        # display preparation that resets automatic axes.
        Makie.save(output,plot.figure;backend=cairo.CairoMakie,update=false)
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
