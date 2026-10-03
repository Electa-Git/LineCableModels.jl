# Compare the architecture baseline and the preservation locks with their versions at
# a git revision. Run `julia test/tools/baseline_ratchet.jl REF`. Every table of the
# baseline is a ceiling. `preservation.toml` declares the direction of each table in
# `[directions]`. A ceiling key added or raised since `REF` fails, and a floor key
# removed or lowered since `REF` fails. Moves the other way pass. The measured ceilings
# `[jet]` and `[allocations]` may also rise in a change of `[environment] julia`, the
# version they were recorded on. The check passes for
# a file absent at `REF`, and a table absent at `REF` is not compared. File renames that
# git detects between `REF` and the working tree, and the module renames they imply,
# are applied to the keys at `REF` first.
using TOML

const BASELINE = "test/quality/architecture_baseline.toml"
const PRESERVATION = "test/quality/preservation.toml"
# Tables of the preservation file that describe it rather than hold counts.
const METADATA = ("directions", "environment")
# Ceilings measured on one Julia version.
const MEASURED = ("jet", "allocations")
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

# The direction, "floor" or "ceiling", of each table of `document`, a version of `file`.
function directions(file, document)
    file == BASELINE && return Dict(table => "ceiling" for table in keys(document))
    declared = Dict{String, String}(get(document, "directions", Dict{String, Any}()))
    for table in keys(document)
        table in METADATA || haskey(declared, table) ||
            error("$file declares no direction for [$table]")
    end
    all(in(("floor", "ceiling")), values(declared)) ||
        error("Each direction in $file is \"floor\" or \"ceiling\"")
    return declared
end

# The entries of `file` that moved against their table's direction since `reference`,
# as lines, and the tables new since `reference`, or nothing when `reference` has no
# `file`. A new table belongs to a lock introduced after `reference`, and its keys are
# not compared.
function moved(repository, reference, file)
    success(pipeline(git(repository, "rev-parse", "--verify", "--quiet",
        reference * "^{commit}"); stdout = devnull)) ||
        error("Unknown git revision: $reference")
    success(pipeline(git(repository, "cat-file", "-e", "$reference:$file");
        stderr = devnull)) || return nothing
    renamed = renames(repository, reference)
    document = TOML.parse(read(git(repository, "show", "$reference:$file"), String))
    current = TOML.parsefile(joinpath(repository, file))
    direction = merge(directions(file, document), directions(file, current))
    julia(d) = get(get(d, "environment", Dict{String, Any}()), "julia", nothing)
    rerecorded = file == PRESERVATION && julia(document) != julia(current)
    values_of(d) = Dict(key => value for (key, value) in entries(d)
        if first(split(key, " | ")) ∉ METADATA)
    before = Dict{String, Any}()
    for (key, value) in values_of(document)
        mergewith!(+, before, Dict(rename_key(key, renamed) => value))
    end
    now = values_of(current)
    tables = sort!([table for table in keys(current)
        if table ∉ METADATA && !haskey(document, table)])
    lines = String[]
    for key in sort!(collect(union(keys(now), keys(before))))
        table = first(split(key, " | "))
        table in tables && continue
        if direction[table] == "ceiling"
            rerecorded && table in MEASURED && continue
            haskey(now, key) || continue
            haskey(before, key) ||
                (push!(lines, string(key, ": ", now[key], " (absent at ", reference, ")")); continue)
        else
            haskey(before, key) || continue
            haskey(now, key) ||
                (push!(lines, string(key, ": removed (", before[key], " at ", reference, ")")); continue)
        end
        (direction[table] == "ceiling" ? now[key] > before[key] : now[key] < before[key]) &&
            push!(lines, string(key, ": ", now[key], " (", before[key], " at ", reference, ")"))
    end
    return (; lines, tables)
end

# The architecture baseline keys added or raised since `reference`.
grown(repository, reference) = moved(repository, reference, BASELINE)

function main(arguments)
    length(arguments) == 1 || error("Usage: julia test/tools/baseline_ratchet.jl REF")
    reference = only(arguments)
    status = 0
    for file in (BASELINE, PRESERVATION)
        result = moved(REPOSITORY, reference, file)
        if result === nothing
            println("No $file at $reference. Nothing to compare.")
            continue
        end
        isempty(result.tables) || println(file, ": tables introduced since ", reference, ": ",
            join(result.tables, ", "), ".")
        if isempty(result.lines)
            count = length(entries(TOML.parsefile(joinpath(REPOSITORY, file))))
            println(file, ": ", count, " entries, none moved against its direction since ",
                reference, ".")
        else
            println(stderr, file, ": ceilings may only fall and floors only rise. ",
                "Entries moved against their direction since ", reference, ":")
            foreach(line -> println(stderr, "  ", line), result.lines)
            status = 1
        end
    end
    return status
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
