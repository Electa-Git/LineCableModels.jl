# Compare native values consistently at configuration and result validation.
_same_value(observed::Number, requested::Number) = observed == requested
# Numeric choice fields (Numf, NFP) read back as numbers, named ones as labels.
_same_value(observed::Number, requested::AbstractString) =
    tryparse(Float64, requested) == observed
_same_value(observed, requested::AbstractString) =
    replace(uppercase(string(observed)), r"[^A-Z0-9]" => "") ==
    replace(uppercase(requested), r"[^A-Z0-9]" => "")
_same_value(observed::AbstractDict, requested::AbstractDict) =
    keys(observed) == keys(requested) &&
    all(_same_value(observed[field], value) for (field, value) in requested)
_same_value(observed, requested) = string(observed) == string(requested)

function _output_matches(roots, suffix)
    matches=String[]
    for root in unique(filter(isdir, roots))
        for (directory, _, files) in walkdir(root), file in files

            endswith(lowercase(file), suffix) &&
                push!(matches, realpath(joinpath(directory, file)))
        end
    end
    unique!(matches)
    return matches
end

function _build_messages(messages::AbstractVector, path, case_directory; console::IO = stdout)
    open(path, "w") do io
        for message in messages
            println(io, join((message.status, message.label, message.name,
                replace(message.text, '\n' => ' ', '\r' => ' ')), '\t'))
        end
    end
    output = dirname(path)
    logs = setdiff(_output_matches([case_directory], ".log"),
        _output_matches([output], ".log"))
    for source in logs
        destination = joinpath(output, relpath(source, case_directory))
        try
            mkpath(dirname(destination))
            cp(source, destination; force = true)
        catch log_error
            diagnostic = "Unable to collect $source: $(sprint(showerror, log_error))"
            println(console, diagnostic)
            flush(console)
            open(io -> println(io, "warning\tbuild\t$source\t$diagnostic"), path, "a")
        end
    end
    errors = filter(message -> lowercase(message.status) == "error", messages)
    if !isempty(errors)
        tline = filter(message -> occursin(r"(?i)\b(tline|line.constants)\b",
            join((message.label, message.name, message.text), ' ')), errors)
        if !isempty(tline)
            alerts = [message.text for message in tline]
            references = [lowercase(basename(replace(matched[1], '\\' => '/')))
                for text in alerts for matched in eachmatch(r"['\"]([^'\"]+\.log)['\"]"i, text)]
            matched_logs = filter(source -> lowercase(basename(source)) in references, logs)
            blocks = NamedTuple{(:stage, :errors, :log_path), Tuple{String, Vector{String}, String}}[]
            for source in matched_logs
                stage = ""
                block_stage = ""
                block_lines = String[]
                enclosed = false
                destination = joinpath(output, relpath(source, case_directory))
                local_log = isfile(destination) ?
                    replace(relpath(source, case_directory), '\\' => '/') : ""
                try
                    for raw in eachline(source)
                        text = strip(raw)
                        if !isempty(text) && all(==('*'), text)
                            if enclosed
                                if any(line -> occursin("ERROR", line), block_lines)
                                    push!(blocks, (stage = block_stage,
                                        errors = copy(block_lines), log_path = local_log))
                                end
                                empty!(block_lines)
                            else
                                block_stage = stage
                            end
                            enclosed = !enclosed
                        elseif enclosed
                            push!(block_lines, String(text))
                        end
                        if !isempty(text) && endswith(text, ":") &&
                           !all(character -> character in ('*', '=', '-', '_', ':'), text)
                            stage = String(text)
                        end
                    end
                catch log_error
                    diagnostic = "Unable to read $source: $(sprint(showerror, log_error))"
                    println(console, diagnostic)
                    flush(console)
                    open(io -> println(io, "warning\tbuild\t$source\t$diagnostic"), path, "a")
                end
            end
            native_errors = isempty(blocks) ?
                [alerts; "No Tline ERROR block found in the referenced log."] :
                [line for block in blocks for line in block.errors]
            stage = isempty(blocks) ? "" : first(blocks).stage
            log_path = isempty(blocks) ? "" : first(blocks).log_path
            if isempty(blocks) && !isempty(matched_logs)
                relative = relpath(first(matched_logs), case_directory)
                isfile(joinpath(output, relative)) && (log_path = replace(relative, '\\' => '/'))
            end
            open(joinpath(output, "failure.toml"), "w") do io
                TOML.print(io, Dict("kind" => "line_constants_error", "stage" => stage,
                    "errors" => native_errors, "log_path" => log_path,
                    "blocks" => [Dict("stage" => block.stage, "errors" => block.errors,
                        "log_path" => block.log_path) for block in blocks]); sorted = true)
            end
        end
        error("PSCAD build messages contain error status; see $path\n" *
            join([message.text for message in errors], '\n'))
    end
    return messages
end

function _data_rows(path::AbstractString; fortran::Bool = false)
    return count(eachline(path)) do line
        fields=split(strip(fortran ? replace(line, 'D' => 'E', 'd' => 'e') : line))
        length(fields) >= 2 || return false
        return tryparse(Float64, fields[1]) !== nothing &&
               tryparse(Float64, fields[2]) !== nothing
    end
end

"""
Wait for complete output rows for at most `timeout` \\[s\\], checking every
`poll_interval` \\[s\\]. Require the requested row count before returning a path.
"""
function _wait_output(
        roots,
        suffix::AbstractString,
        expected_rows::Integer;
        timeout::Real = 30,
        poll_interval::Real = 0.1
)
    expected_rows > 0 || throw(ArgumentError("expected_rows must be positive"))
    timeout > 0 || throw(ArgumentError("timeout must be positive"))
    poll_interval > 0 || throw(ArgumentError("poll_interval must be positive"))
    deadline=time() + timeout
    previous=nothing
    observed="no matching file"
    while time() < deadline
        matches=_output_matches(roots, suffix)
        length(matches) <= 1 || throw(ArgumentError(
            "PSCAD emitted $(length(matches)) files ending in $suffix",
        ))
        if length(matches) == 1
            path=only(matches)
            rows=_data_rows(path)
            rows <= expected_rows || throw(ArgumentError(
                "PSCAD output $path contains $rows rows; expected $expected_rows",
            ))
            state=(path, filesize(path), rows)
            observed="$rows of $expected_rows rows in $path"
            state == previous && rows == expected_rows && return path
            previous=state
        end
        sleep(poll_interval)
    end
    throw(ArgumentError(
        "PSCAD output $suffix did not become complete within $timeout seconds; " *
        "last observed $observed",
    ))
end

# INF parsing: PowerImpedance.jl `src/PSCAD/remote/files.jl:23-51`.
# Legacy INF files are Windows-1252.
function _inf_entries(path::AbstractString)
    bytes = read(path)
    special = ('€', '\u0081', '‚', 'ƒ', '„', '…', '†', '‡', 'ˆ', '‰', 'Š', '‹',
        'Œ', '\u008d', 'Ž', '\u008f', '\u0090', '‘', '’', '“', '”', '•', '–', '—',
        '˜', '™', 'š', '›', 'œ', '\u009d', 'ž', 'Ÿ')
    text = join(byte in 0x80:0x9f ? special[Int(byte) - 0x7f] : Char(byte)
    for byte in bytes)
    entries = NamedTuple{(:index, :desc, :group), Tuple{Int, String, String}}[]
    for line in split(text, '\n')
        desc = match(r"\bDesc=\"([^\"]*)\"", line)
        desc === nothing && continue
        isempty(desc[1]) &&
            throw(ArgumentError("PSCAD INF contains an empty channel name: $path"))
        pgb = match(r"PGB\((\d+)\)", line)
        group = match(r"\bGroup=\"([^\"]*)\"", line)
        push!(entries,
            (index = pgb === nothing ? length(entries) + 1 : parse(Int, pgb[1]),
                desc = String(desc[1]), group = group === nothing ? "" : String(group[1])))
    end
    isempty(entries) && throw(ArgumentError("PSCAD INF contains no channels: $path"))
    sort!(entries; by = entry -> entry.index)
    [entry.index for entry in entries] == collect(1:length(entries)) ||
        throw(ArgumentError("PSCAD INF channel indices must be consecutive from 1: $path"))
    return entries
end

# Output waiting: PowerImpedance.jl `src/PSCAD/remote/files.jl:77-137`.
function _wait_outputs(roots, stem::AbstractString, expected_rows::Integer;
        timeout::Real = 60, poll_interval::Real = 0.25)
    expected_rows > 0 || throw(ArgumentError("expected_rows must be positive"))
    isfinite(timeout) && timeout > 0 ||
        throw(ArgumentError("timeout must be positive and finite"))
    isfinite(poll_interval) && poll_interval > 0 ||
        throw(ArgumentError("poll_interval must be positive and finite"))
    deadline = time() + timeout
    previous = nothing
    counted = Set{Tuple}()
    rows = Int[]
    found = String[]
    matches = String[]
    while time() < deadline
        matches = filter(path -> lowercase(basename(path)) == lowercase(stem * ".inf"),
            _output_matches(roots, ".inf"))
        length(matches) <= 1 ||
            throw(ArgumentError("PSCAD emitted multiple $(stem).inf files: $matches"))
        if length(matches) == 1
            inf = only(matches)
            # PSCAD can leave INF metadata incomplete while it writes the file.
            count = try
                length(_inf_entries(inf))
            catch
                0
            end
            outs = [joinpath(dirname(inf), stem * "_" * lpad(string(index), 2, '0') *
                                           ".out")
                    for index in 1:cld(count, 10)]
            found = filter(isfile, outs)
            state = (inf, count, filesize(inf), found, filesize.(found))
            if count > 0 && length(found) == length(outs) && state == previous &&
               state ∉ counted
                # Row counts run once per stable file state to avoid repeated full-file scans.
                push!(counted, state)
                rows = [_data_rows(path; fortran = true) for path in outs]
                any(>(expected_rows), rows) && break
                if state[3] == filesize(inf) && state[5] == filesize.(outs) &&
                   all(==(expected_rows), rows)
                    return (; inf, outs)
                end
                previous = nothing
                sleep(poll_interval)
                continue
            end
            previous = state
        end
        sleep(poll_interval)
    end
    if isempty(counted)
        rows = [_data_rows(path; fortran = true) for path in found]
    end
    throw(ArgumentError("PSCAD outputs did not become complete within $timeout seconds; " *
        "last observed INF files found: $matches; OUT files found: $found; " *
        "rows per file: $rows; expected $expected_rows"))
end
