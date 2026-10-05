# Compare the architecture baseline and the preservation locks with their versions at
# a git revision. Run `julia test/tools/baseline_ratchet.jl REF`. Every table of the
# baseline is a ceiling. `preservation.toml` declares the direction of each table in
# `[directions]`. A ceiling key added or raised since `REF` fails, and a floor key
# removed or lowered since `REF` fails. Moves the other way pass. The measured ceilings
# `[jet]` and `[allocations]` may also rise in a change of `[environment]`, the Julia
# version and `Manifest.toml` they were recorded with. The check passes for a file
# absent at `REF`, and a table absent at `REF` is not compared. File renames that
# git detects between `REF` and the working tree, and the module renames they imply,
# are applied to the keys at `REF` first.
#
# An `[inferred]` floor can fall when `@inferred` lines move to other files. Each removed
# line appears again in another file, unchanged apart from indentation, and these lines
# cover the drop.
using TOML

const BASELINE = "test/quality/architecture_baseline.toml"
const PRESERVATION = "test/quality/preservation.toml"
# Tables of the preservation file that describe it rather than hold counts.
const METADATA = ("directions", "environment")
# Ceilings measured in one environment.
const MEASURED = ("jet", "allocations")
# Floors that count, per test file, the lines that contain a marker. These lines can move
# to another file.
const RELOCATABLE = Dict("inferred" => "@inferred")
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

# The lines removed from each file since `reference` and the lines added to each file,
# without indentation. An untracked file counts as added in full.
function changed_lines(repository, reference)
    removed, added = Dict{String, Vector{String}}(), Dict{String, Vector{String}}()
    record!(lines, file, line) =
        push!(get!(Vector{String}, lines, file), String(lstrip(line)))
    path(line, prefix) =
        (name = line[5:end]; name == "/dev/null" ? "" : chopprefix(name, prefix))
    old, new, hunk = "", "", false
    for line in eachline(git(repository, "diff", "-M", "--unified=0", "--no-color",
            "--no-ext-diff", "--src-prefix=a/", "--dst-prefix=b/", reference))
        if startswith(line, "diff --git ")
            hunk = false
        elseif hunk
            startswith(line, '-') && record!(removed, isempty(new) ? old : new, line[2:end])
            startswith(line, '+') && record!(added, new, line[2:end])
        elseif startswith(line, "--- ")
            old = path(line, "a/")
        elseif startswith(line, "+++ ")
            new = path(line, "b/")
        elseif startswith(line, "@@")
            hunk = true
        end
    end
    untracked = git(repository, "ls-files", "--others", "--exclude-standard", "-z")
    for name in split(read(untracked, String), '\0'; keepempty = false)
        file = joinpath(repository, name)
        isfile(file) && foreach(line -> record!(added, String(name), line), eachline(file))
    end
    return (; removed, added)
end

# Whether each line with `marker` removed from `file` appears again as an added line of
# another file. Each added line matches one removed line, and the matched lines contain
# `marker` at least `drop` times. `unmatched` counts the available added lines per file.
function relocated!(unmatched, removed, file, marker, drop)
    found, files = 0, sort!(collect(keys(unmatched)))
    for line in get(removed, file, String[])
        occursin(marker, line) || continue
        index = findfirst(name -> name != file && get(unmatched[name], line, 0) > 0, files)
        index === nothing && return false
        unmatched[files[index]][line] -= 1
        found += count(marker, line)
    end
    return found >= drop
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
    environment(d) = get(d, "environment", Dict{String, Any}())
    rerecorded = file == PRESERVATION && environment(document) != environment(current)
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
    changes = nothing
    unmatched = Dict{String, Dict{String, Int}}()
    function relocated(table, path, drop)
        haskey(RELOCATABLE, table) || return false
        if changes === nothing
            changes = changed_lines(repository, reference)
            for (name, added) in changes.added, line in added
                counts = get!(Dict{String, Int}, unmatched, name)
                counts[line] = get(counts, line, 0) + 1
            end
        end
        return relocated!(unmatched, changes.removed, path, RELOCATABLE[table], drop)
    end
    for key in sort!(collect(union(keys(now), keys(before))))
        table, entry = split(key, " | "; limit = 2)
        table in tables && continue
        if direction[table] == "ceiling"
            rerecorded && table in MEASURED && continue
            haskey(now, key) || continue
            haskey(before, key) ||
                (push!(lines, string(key, ": ", now[key], " (absent at ", reference, ")")); continue)
        else
            haskey(before, key) || continue
            drop = before[key] - get(now, key, 0)
            drop > 0 && relocated(table, entry, drop) && continue
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
