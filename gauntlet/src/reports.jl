function _verdict_name(value::Verdict)
    return lowercase(String(nameof(typeof(value))))
end

_check_name(::MatrixCheck{Q}) where {Q} = "matrix.$Q"
_check_name(::ModalCheck{Q}) where {Q} = "modal.$Q"
_check_name(::FitCheck{Q}) where {Q} = "fit.$Q"
_check_name(::PhysicalCheck{Q}) where {Q} = "physical.$Q"
_check_name(::PerformanceCheck) = "performance"

function Tables.rows(report::Report)
    return (
        (
            case_id = trial.case.id,
            family = lowercase(String(nameof(typeof(trial.case.family)))),
            fidelity = lowercase(String(nameof(typeof(trial.case.fidelity)))),
            check = _check_name(comparison.check),
            verdict = _verdict_name(comparison.verdict),
            max_abs = comparison.metrics === nothing ? missing : comparison.metrics.max_abs,
            max_rel = comparison.metrics === nothing ? missing : comparison.metrics.max_rel,
            worst_frequency = comparison.metrics === nothing ? missing :
                              comparison.metrics.worst_frequency,
            worst_row = comparison.metrics === nothing ? missing :
                        comparison.metrics.worst_entry[1],
            worst_column = comparison.metrics === nothing ? missing :
                           comparison.metrics.worst_entry[2],
            minimum_seconds = trial.performance === nothing ? missing :
                              trial.performance.minimum_seconds,
            median_seconds = trial.performance === nothing ? missing :
                             trial.performance.median_seconds,
            bytes = trial.performance === nothing ? missing : trial.performance.bytes,
            allocations = trial.performance === nothing ? missing :
                          trial.performance.allocations,
            frequency_count = trial.performance === nothing ? missing :
                              trial.performance.frequencies,
            conductor_count = trial.performance === nothing ? missing :
                              trial.performance.conductors,
            frequency_points_per_second = trial.performance === nothing ? missing :
                                          trial.performance.points_per_second,
            julia = trial.performance === nothing ? missing : trial.performance.julia,
            os = trial.performance === nothing ? missing : trial.performance.os,
            arch = trial.performance === nothing ? missing : trial.performance.arch,
            threads = trial.performance === nothing ? missing : trial.performance.threads,
            blas = trial.performance === nothing ? missing : trial.performance.blas,
            bytes_ratio = trial.performance === nothing ? missing :
                          something(trial.performance.bytes_ratio, missing),
            allocations_ratio = trial.performance === nothing ? missing :
                                something(trial.performance.allocations_ratio, missing),
            detail = comparison.detail
        ) for trial in report.trials for comparison in trial.comparisons
    )
end

Tables.istable(::Type{<:Report}) = true
Tables.rowaccess(::Type{<:Report}) = true
Tables.schema(report::Report) = Tables.schema(Tables.rows(report))

function _json_rows(report::Report)
    return collect(Tables.rows(report))
end

"""Write a structured Gauntlet report selected by its filename extension."""
function write_report(path::AbstractString, report::Report)
    extension = lowercase(splitext(path)[2])
    extension == ".json" && return _write_json(path, report)
    extension == ".csv" && return _write_csv(path, report)
    extension in (".md", ".markdown") && return _write_markdown(path, report)
    throw(ArgumentError("report extension must be .json, .csv, or .md"))
end

function _write_json(path, report)
    open(path, "w") do io
        JSON3.pretty(io, _json_rows(report))
        println(io)
    end
    return abspath(path)
end

function _csv_field(value)
    text = ismissing(value) ? "" : string(value)
    return occursin(r"[\",\n]", text) ? "\"$(replace(text, '"' => "\"\""))\"" : text
end

function _write_csv(path, report)
    rows = collect(Tables.rows(report))
    names = isempty(rows) ? Symbol[] : collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(names, ','))
        for row in rows
            println(io, join((_csv_field(getproperty(row, name)) for name in names), ','))
        end
        performance = filter(row -> !ismissing(row.minimum_seconds), rows)
        if !isempty(performance)
            println(io, "\n## Performance\n")
            println(io,
                "| Case | Minimum [ms] | Median [ms] | Bytes | Allocations | Points/s |"
            )
            println(io, "|---|---:|---:|---:|---:|---:|")
            for row in performance
                println(io,
                    @sprintf("| `%s` | %.4f | %.4f | %d | %d | %.3f |",
                        row.case_id,
                        1e3 * row.minimum_seconds,
                        1e3 * row.median_seconds,
                        row.bytes,
                        row.allocations,
                        row.frequency_points_per_second))
            end
        end
    end
    return abspath(path)
end

function _write_markdown(path, report)
    rows = collect(Tables.rows(report))
    open(path, "w") do io
        println(io, "# Gauntlet report: `$(report.suite.name)`\n")
        println(io, "Started: $(report.started_at)  ")
        println(io, "Finished: $(report.finished_at)\n")
        println(io, "| Case | Family | Fidelity | Check | Verdict | Max relative error |")
        println(io, "|---|---|---|---|---|---:|")
        for row in rows
            relative = ismissing(row.max_rel) ? "—" : @sprintf("%.6g", row.max_rel)
            println(io,
                "| `$(row.case_id)` | $(row.family) | $(row.fidelity) | " *
                "$(row.check) | $(row.verdict) | $relative |")
        end
    end
    return abspath(path)
end
