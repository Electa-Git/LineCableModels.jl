# The displacement of a change, from the git revision `REF` to the working tree. Run
# `julia --startup-file=no test/tools/displacement.jl REF` from the repository root. It
# reads only git and the working tree, and it never loads the package. Its exit status
# reports only a usage error.
#
# For the `.jl` files of `src` and `ext`, it prints the private definitions and private
# types added and removed, and the lines added, removed and net. Then it prints the
# `absence` ceilings of `test/quality/architecture_baseline.toml` per module directory, at
# `REF` and now.
#
# A private definition is a top-level line that defines a function, a constant or a
# struct whose name starts with `_`. A private type is a `struct _X`. A definition that
# moves within a file or to another file counts as removed and added. Untracked files
# count as added.
using TOML

const REPOSITORY = normpath(joinpath(@__DIR__, "..", ".."))
const BASELINE = "test/quality/architecture_baseline.toml"
const DEFINITION = r"^(function |)_[A-Za-z0-9_!]+\(|^const _[A-Za-z0-9_]+ =|^(mutable )?struct _"
const TYPE = r"^(?:mutable )?struct (_[A-Za-z0-9_]+)"

git(arguments...) = Cmd(["git", "-C", REPOSITORY, arguments...])

# The lines of the `.jl` files of `src` and `ext` that the working tree adds and removes
# since `reference`.
function changed_lines(reference)
    added, removed = String[], String[]
    file, hunk = "", false
    for line in eachline(git("diff", "--no-color", "--no-ext-diff", "--unified=0",
            reference, "--", "src", "ext"))
        if startswith(line, "diff --git ")
            file, hunk = line, false
        elseif startswith(line, "@@")
            hunk = true
        elseif hunk && endswith(file, ".jl")
            startswith(line, '+') && push!(added, line[2:end])
            startswith(line, '-') && push!(removed, line[2:end])
        end
    end
    untracked = git("ls-files", "--others", "--exclude-standard", "-z", "--", "src", "ext")
    for name in split(read(untracked, String), '\0'; keepempty = false)
        endswith(name, ".jl") && append!(added, readlines(joinpath(REPOSITORY, name)))
    end
    return (; added, removed)
end

definitions(lines) = count(line -> occursin(DEFINITION, line), lines)
types(lines) = [match(TYPE, line)[1] for line in lines if occursin(TYPE, line)]
listed(names) = isempty(names) ? "" : string(" (", join(names, ", "), ")")
total(table) = sum(values(table); init = 0)

function main(arguments)
    length(arguments) == 1 ||
        error("Usage: julia --startup-file=no test/tools/displacement.jl REF")
    reference = only(arguments)
    success(pipeline(git("rev-parse", "--verify", "--quiet", reference * "^{commit}");
        stdout = devnull)) || error("Unknown git revision: $reference")
    lines = changed_lines(reference)
    println("Displacement from ", reference, " to the working tree, src and ext:")
    println("private definitions: +", definitions(lines.added), " -",
        definitions(lines.removed))
    added, removed = types(lines.added), types(lines.removed)
    println("private types: +", length(added), listed(added), " -", length(removed),
        listed(removed))
    println("lines: +", length(lines.added), " -", length(lines.removed), ", net ",
        length(lines.added) - length(lines.removed))
    now = get(TOML.parsefile(joinpath(REPOSITORY, BASELINE)), "absence", Dict{String, Any}())
    recorded = success(pipeline(git("cat-file", "-e", "$reference:$BASELINE");
        stderr = devnull))
    earlier = recorded ? TOML.parse(read(git("show", "$reference:$BASELINE"), String)) :
        Dict{String, Any}()
    if haskey(earlier, "absence")
        before = earlier["absence"]
        println("absence ceilings at ", reference, " and now: ", total(before), " -> ",
            total(now))
        for key in sort!(collect(union(keys(before), keys(now))))
            println("  ", key, ": ", get(before, key, 0), " -> ", get(now, key, 0))
        end
    else
        println("absence ceilings: table absent at ", reference, ", now ", total(now))
        foreach(key -> println("  ", key, ": ", now[key]), sort!(collect(keys(now))))
    end
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
