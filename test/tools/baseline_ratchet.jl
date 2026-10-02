# Compare the architecture baseline with its version at a git revision. Run
# `julia test/tools/baseline_ratchet.jl REF`. A key added or a count raised since
# `REF` fails. Deleted keys and lowered counts pass. The check passes when `REF`
# has no baseline file, and a table absent at `REF` is not compared. File renames that git detects between `REF` and the working
# tree, and the module renames they imply, are applied to the keys at `REF` first.
using TOML

const BASELINE = "test/quality/architecture_baseline.toml"
const REPOSITORY = normpath(joinpath(@__DIR__, "..", ".."))

git(repository, arguments::AbstractString...) = Cmd(["git", "-C", repository, arguments...])

# Nested tables become keys of the form "table | key | counter".
function entries(document, prefix = "")
    flat = Dict{String, Any}()
    for (key, value) in document
        name = isempty(prefix) ? key : string(prefix, " | ", key)
        value isa AbstractDict ? merge!(flat, entries(value, name)) : (flat[name] = value)
    end
    return flat
end

declares(text, name) = occursin(Regex("(?m)^\\s*(?:bare)?module\\s+" * name * "\\b"), text)

# Renamed files, and the module renames of renamed entry files `<Module>.jl` that
# declare their module on both sides.
function renames(repository, reference)
    fields = split(read(git(repository, "diff", "-M", "--name-status", "-z", reference),
        String), '\0'; keepempty = false)
    files = Dict{String, String}()
    i = 1
    while i <= length(fields)
        if startswith(fields[i], 'R')
            files[fields[i+1]] = fields[i+2]
            i += 3
        else
            i += 2
        end
    end
    modules = Dict{String, String}()
    for (old, new) in files
        (before, extension), after = splitext(basename(old)), first(splitext(basename(new)))
        extension == ".jl" && before != after || continue
        declares(read(git(repository, "show", "$reference:$old"), String), before) &&
            declares(read(joinpath(repository, new), String), after) &&
            (modules[before] = after)
    end
    return (; files, modules)
end

# Path components change only through file renames. Module names change in the
# other components.
function rename_key(key, renamed)
    parts = map(split(key, " | ")) do part
        haskey(renamed.files, part) && return renamed.files[part]
        occursin('/', part) && return part
        for (old, new) in renamed.modules
            part = replace(part, Regex("(?<![\\w!])" * old * "(?![\\w!])") => new)
        end
        return part
    end
    return join(parts, " | ")
end

# The keys added or raised since `reference`, as lines, and the tables new since
# `reference`, or nothing when `reference` has no baseline. A new table belongs
# to a guard introduced after `reference`, and its keys are not compared.
function grown(repository, reference)
    success(pipeline(git(repository, "rev-parse", "--verify", "--quiet",
        reference * "^{commit}"); stdout = devnull)) ||
        error("Unknown git revision: $reference")
    success(pipeline(git(repository, "cat-file", "-e", "$reference:$BASELINE");
        stderr = devnull)) || return nothing
    renamed = renames(repository, reference)
    document = TOML.parse(read(git(repository, "show", "$reference:$BASELINE"), String))
    before = Dict{String, Any}()
    for (key, value) in entries(document)
        mergewith!(+, before, Dict(rename_key(key, renamed) => value))
    end
    current = TOML.parsefile(joinpath(repository, BASELINE))
    tables = sort!([table for table in keys(current) if !haskey(document, table)])
    lines = sort!([haskey(before, key) ?
        string(key, ": ", value, " (", before[key], " at ", reference, ")") :
        string(key, ": ", value, " (absent at ", reference, ")")
        for (key, value) in entries(current)
        if first(split(key, " | ")) ∉ tables && (!haskey(before, key) || value > before[key])])
    return (; lines, tables)
end

function main(arguments)
    length(arguments) == 1 || error("Usage: julia test/tools/baseline_ratchet.jl REF")
    reference = only(arguments)
    result = grown(REPOSITORY, reference)
    if result === nothing
        println("No architecture baseline at $reference. Nothing to compare.")
        return 0
    end
    isempty(result.tables) || println("Tables introduced since ", reference, ": ",
        join(result.tables, ", "), ".")
    if !isempty(result.lines)
        println(stderr, "The architecture baseline may only shrink. Entries added or raised since ",
            reference, ":")
        foreach(line -> println(stderr, "  ", line), result.lines)
        return 1
    end
    count = length(entries(TOML.parsefile(joinpath(REPOSITORY, BASELINE))))
    println("Architecture baseline: ", count, " entries, none added or raised since ",
        reference, ".")
    return 0
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
