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

function _data_rows(path::AbstractString)
    return count(eachline(path)) do line
        fields=split(strip(line))
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
