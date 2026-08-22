module CoverageGate

using CoverageTools
using Printf

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const SOURCE_DIRECTORIES = ("src", "ext")
const GAUNTLET_DIRECTORY = joinpath(REPOSITORY_ROOT, "test", "gauntlet")
const GAUNTLET_CASE_DIRECTORY = joinpath(GAUNTLET_DIRECTORY, "cases")
const DEFAULT_MINIMUM = 0.95

function source_files()
    files = String[]
    for relative_directory in SOURCE_DIRECTORIES
        directory = joinpath(REPOSITORY_ROOT, relative_directory)
        for (root, _, names) in walkdir(directory)
            append!(files, joinpath(root, name) for name in names if endswith(name, ".jl"))
        end
    end
    return sort!(normpath.(files))
end

function clean_traces!()
    removed = String[]
    for relative_directory in (SOURCE_DIRECTORIES..., "test")
        directory = joinpath(REPOSITORY_ROOT, relative_directory)
        for (root, _, names) in walkdir(directory)
            for name in names
                if CoverageTools.iscovfile(name)
                    candidate = joinpath(root, name)
                    rm(candidate; force = true)
                    push!(removed, candidate)
                end
            end
        end
    end
    report = joinpath(REPOSITORY_ROOT, "lcov.info")
    if isfile(report)
        rm(report)
        push!(removed, report)
    end
    return removed
end

function collect_coverage()
    coverage = CoverageTools.FileCoverage[]
    for relative_directory in SOURCE_DIRECTORIES
        append!(coverage, CoverageTools.process_folder(
            joinpath(REPOSITORY_ROOT, relative_directory),
        ))
    end
    return coverage
end

function collect_gauntlet_coverage()
    isdir(GAUNTLET_DIRECTORY) || return CoverageTools.FileCoverage[]
    coverage = CoverageTools.process_folder(GAUNTLET_DIRECTORY)
    return filter(coverage) do record
        path = normpath(record.filename)
        !startswith(path, GAUNTLET_CASE_DIRECTORY * Base.Filesystem.path_separator)
    end
end

function assert_complete_inventory(coverage)
    represented = Set(normpath(record.filename) for record in coverage)
    missing = setdiff(Set(source_files()), represented)
    isempty(missing) || error(
        "Coverage report omitted production files:\n" * join(sort!(collect(missing)), "\n"),
    )
    return nothing
end

function check_coverage(; minimum = DEFAULT_MINIMUM, output = "lcov.info")
    coverage = collect_coverage()
    assert_complete_inventory(coverage)
    covered, total = CoverageTools.get_summary(coverage)
    ratio = total == 0 ? 0.0 : covered / total
    @printf("Production line coverage: %d/%d (%.2f%%)\n", covered, total, 100ratio)
    published = vcat(coverage, collect_gauntlet_coverage())
    CoverageTools.LCOV.writefile(joinpath(REPOSITORY_ROOT, output), published)
    ratio >= minimum || error(
        @sprintf("Production line coverage %.2f%% is below the required %.2f%%",
        100ratio,
        100minimum,),
    )
    return ratio
end

function main(arguments)
    command = isempty(arguments) ? "check" : only(arguments)
    if command == "clean"
        removed = clean_traces!()
        println("Removed $(length(removed)) coverage trace files.")
    elseif command == "check"
        minimum = parse(Float64, get(ENV, "LINECABLEMODELS_MIN_COVERAGE", "0.95"))
        check_coverage(; minimum)
    else
        error("Unknown coverage command '$command'; expected 'clean' or 'check'.")
    end
    return nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CoverageGate.main(ARGS)
end
