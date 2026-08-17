const _HARVEST_ROOT = normpath(joinpath(_GAUNTLET_ROOT, "harvest", "pscad"))
const _PSCAD_LIVE_PROJECT = joinpath(_HARVEST_ROOT, "live")
const _LIVE_COMMANDS = ("inspect", "case", "batch")

function _cli_options(arguments)
    positional = String[]
    options = Dict{String, Any}()
    index = 1
    while index <= length(arguments)
        token = arguments[index]
        if startswith(token, "--")
            key = token[3:end]
            if index == length(arguments) || startswith(arguments[index + 1], "--")
                options[key] = true
            else
                options[key] = arguments[index + 1]
                index += 1
            end
        else
            push!(positional, token)
        end
        index += 1
    end
    return positional, options
end

function _cli_help(io::IO = stdout)
    print(io, """
    linecablebenchmark <command> [arguments]

    Dataset commands:
      pending <dataset.json>
      verify <dataset.json> <amendment-dir>
      ingest <source-dir> <destination> --amendments <dir>
      catalog <donor-dir>
      sources <root> [--verify-only] [--no-extract]

    Live PSCAD commands (activate the PythonCall extension on demand):
      inspect <project.pscx> [--config <benchmark.toml>]
      case <project.pscx> [--output <dir>] [--dry-run]
      batch <donor-dir> [--output <dir>] [--dry-run]
    """)
end

function _activate_live_adapter()
    extension = Base.get_extension(Gauntlet, :GauntletPythonCallExt)
    extension === nothing || return extension
    try
        Base.require(Base.PkgId(
            Base.UUID("6099a3de-0909-46bc-b1f4-468b9a2dfc0d"), "PythonCall"
        ))
    catch error
        error isa InterruptException && rethrow()
        throw(ErrorException(
            "live PSCAD automation could not load PythonCall: $(sprint(showerror, error))",
        ))
    end
    extension = Base.get_extension(Gauntlet, :GauntletPythonCallExt)
    extension === nothing && throw(ErrorException(
        "PythonCall loaded, but the Gauntlet PSCAD extension did not activate",
    ))
    return extension
end

function _delegate_live(arguments)
    script = "using PythonCall; using Gauntlet; exit(Gauntlet.harvest_main(ARGS))"
    command = `$(Base.julia_cmd()) --startup-file=no --project=$(_PSCAD_LIVE_PROJECT) -e $script -- $arguments`
    process = run(addenv(
        ignorestatus(command),
        "JULIA_LOAD_PATH" => "@:@stdlib",
        "LINECABLEBENCHMARK_LIVE" => "1"
    ))
    return process.exitcode
end

function _print_cli(value, options)
    output = get(options, "json", nothing)
    if output === nothing
        JSON3.pretty(stdout, value)
        println()
    else
        _write_pscad_json(String(output), value)
        println(abspath(String(output)))
    end
end

function _report_dict(report::HarvestReport)
    return Dict(
        "expected" => report.expected,
        "verified" => report.verified,
        "valid" => isempty(report.errors),
        "errors" => report.errors,
        "amendments" => [Dict(
             "case_id" => row.case_id,
             "specification_sha256" => row.specification_sha256,
             "manifest" => row.manifest,
             "matrix_dimension" => row.matrix_dimension
         ) for row in report.amendments]
    )
end

function _pending_dict(rows)
    return Dict("count" => length(rows),
        "records" => [Dict(
             "case_id" => row.case_id,
             "specification_sha256" => row.specification_sha256,
             "source_sha256" => row.source_sha256,
             "manifest" => row.manifest,
             "missing" => row.missing
         ) for row in rows])
end

function harvest_main(arguments = ARGS)
    isempty(arguments) && (_cli_help(stderr); return 2)
    first(arguments) in ("-h", "--help", "help") && (_cli_help(); return 0)
    command = first(arguments)
    if command in _LIVE_COMMANDS &&
       !haskey(ENV, "LINECABLEBENCHMARK_LIVE") &&
       Base.get_extension(Gauntlet, :GauntletPythonCallExt) === nothing
        return _delegate_live(arguments)
    end
    positional, options = _cli_options(arguments[2:end])
    result = if command == "pending"
        length(positional) == 1 || throw(ArgumentError("pending requires dataset.json"))
        _pending_dict(pending(PSCAD(), only(positional)))
    elseif command == "verify"
        length(positional) == 2 || throw(ArgumentError(
            "verify requires dataset.json and an amendment directory",
        ))
        _report_dict(verify(PSCAD(), positional[1], positional[2]))
    elseif command == "ingest"
        length(positional) == 2 || throw(ArgumentError(
            "ingest requires a source and destination directory",
        ))
        ingest(PSCAD(), positional[1], positional[2];
            amendments = get(options, "amendments", nothing))
    elseif command == "catalog"
        length(positional) == 1 ||
            throw(ArgumentError("catalog requires a donor directory"))
        catalog(only(positional))
    elseif command == "sources"
        length(positional) == 1 ||
            throw(ArgumentError("sources requires a destination root"))
        harvest_sources(
            only(positional),
            String(get(options, "catalog", joinpath(_HARVEST_ROOT, "sources.toml")));
            download = !haskey(options, "verify-only"),
            extract = !haskey(options, "no-extract")
        )
    elseif command in _LIVE_COMMANDS
        length(positional) == 1 || throw(ArgumentError("$command requires one source path"))
        _activate_live_adapter()
        Base.invokelatest(
            automate,
            PSCAD(),
            Val(Symbol(command)),
            merge(options, Dict("source" => only(positional)))
        )
    else
        throw(ArgumentError("unknown linecablebenchmark command '$command'"))
    end
    _print_cli(result, options)
    return result isa AbstractDict && get(result, "valid", true) === false ? 1 : 0
end
