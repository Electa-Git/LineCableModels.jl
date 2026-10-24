_display_path(path::AbstractString) =
    try
        relpath(abspath(path), pwd())
    catch
        basename(path)
    end

function _json_path(file_name::AbstractString)
    path = isabspath(file_name) ? String(file_name) : abspath(file_name)
    extension = lowercase(splitext(path)[2])
    isempty(extension) && return path * ".json"
    extension == ".json" || throw(ArgumentError(
        "JSON output requires a .json extension; got '$extension'",
    ))
    return path
end
