"""
    identify(config::RemoteConfig)

Read the identity of the PSCAD installation of the station of `config`, without loading a
project. It borrows or starts PSCAD as a run does, and it requires `using PythonCall`.

# Arguments

- `config`: the station and the shared directory.

# Returns

- A `Dict{String, String}`: `schema`, PSCAD's `version`, the caller's `automation_version`
  and `common_version` (the installed `mhi.pscad` and `mhi.common`), and
  `line_constants`, the line-constants program that PSCAD's settings select.

# Notes

Passing this record as the `solver_identity` computation option requires that
installation. A fresh run reads the identity on its own connection before it loads the
project. Completed-run reuse reads the identity of the current server before acceptance.
"""
function identify end

"""
    run_remote_pscad(config, project, output, settings, frequencies; line_name, expected_identity)

Run the line-constants solve of the staged PSCAD `project` on the station of `config`, and
collect its files into the directory `output`. It requires `using PythonCall`.

# Arguments

- `config`: the station and the shared directory.
- `project`: the staged project file, under `config.local_root`.
- `output`: an empty or missing directory for the run's files.
- `settings`: the native settings, `value` and `readback` by field of the `ground`,
  `frequency` and `configuration` rows.
- `frequencies`: the computation frequencies \\[Hz\\].
- `line_name`: the name of the line, and of the detailed files that PSCAD writes into its
  working directory.
- `expected_identity`: a record of [`identify`](@ref) that the server must report before
  the run loads the project.

# Returns

- A named tuple of the compile call's `elapsed_seconds` and `elapsed_scope`, the
  `console_path` of the run's log and the `output_dir`.
"""
function run_remote_pscad end

# The solver identity that a PSCAD automation server reports.
function validate(result::AbstractDict, config::RemoteConfig)
    all(pair -> first(pair) isa AbstractString && last(pair) isa AbstractString, result) ||
        throw(ArgumentError("PSCAD solver identity must contain string fields"))
    get(result, "version", nothing) == PSCAD_VERSION || throw(ArgumentError(
        "PSCAD automation server did not return the PSCAD $PSCAD_VERSION solver identity"))
    get(result, "schema", nothing) == "2" && !isempty(get(result, "line_constants", "")) ||
        throw(ArgumentError("PSCAD automation server returned an incomplete solver identity"))
    return result
end

function _diagnostic_tail(path::AbstractString; count::Integer = 12)
    isfile(path) || return "PSCAD produced no diagnostic log."
    lines=filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return "PSCAD diagnostic log is empty."
    return join(last(lines, min(count, length(lines))), '\n')
end
