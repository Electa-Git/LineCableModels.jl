# L4. Behaviour of the preservation corpus at a git revision and on the working tree,
# compared bitwise. Run from the repository root:
#
#     julia --project=test test/tools/equivalence.jl REF [ALLOWED]
#
# REF is prepared as for `performance.jl`. Each side records its own revision's corpus
# with its own recording program, `test/tools/fingerprint.jl`, in a fresh process; a
# revision without that program uses the working tree's program and corpus, and the
# report says which copy each side used. The recording holds every node of the inputs
# and the result of each scenario, with its type and the bits of its value, and the
# laziness counts of the parametric scenarios (see `fingerprint.jl`).
#
# `ALLOWED` is an optional file of declarations, one per line (`#` starts a comment):
#
#     rename | Old => New
#     scenario | path prefix | value, type or path
#
# A rename is applied to REF's type names, path segments and Symbol values, matching
# whole identifiers only, before comparing. A difference changes a value under the same
# type, changes a type, or exists on one side only (path); a declaration covers the
# differences of its kind whose path starts with its prefix, and an empty prefix covers the
# whole scenario. Declared differences pass; renames and declarations that match nothing
# are reported; any other difference fails. Not run in CI.
const REPOSITORY = dirname(dirname(@__DIR__))
const RECORDER = joinpath("test", "tools", "fingerprint.jl")
const KINDS = (:value, :type, :path)
const SHOWN = 5

Base.include(@__MODULE__, joinpath(@__DIR__, "worktree.jl"))

# The scenarios in recorded order, and each scenario's (path => (type, value)) nodes in
# recorded order.
function fingerprints(project, recorder, file)
    output(`$(Base.julia_cmd()) --startup-file=no --project=$project $recorder $file`)
    order = String[]
    nodes = Dict{String, Vector{Pair{String, Tuple{String, String}}}}()
    for line in eachline(file)
        scenario, path, type, value = String.(split(line, '\t'; limit = 4))
        haskey(nodes, scenario) || (push!(order, scenario); nodes[scenario] = [])
        push!(nodes[scenario], path => (type, value))
    end
    return (; order, nodes)
end

const IDENTIFIER = r"^[A-Za-z_][\w!]*$"

# The renames and the declared differences of an ALLOWED file.
function declarations(file)
    renames = Pair{String, String}[]
    declared = @NamedTuple{scenario::String, prefix::String, kind::Symbol}[]
    file === nothing && return (; renames, declared)
    for (number, raw) in enumerate(eachline(file))
        text = strip(first(split(raw, '#')))
        isempty(text) && continue
        fields = String.(strip.(split(text, '|')))
        if length(fields) == 2 && fields[1] == "rename" && occursin("=>", fields[2])
            old, new = String.(strip.(split(fields[2], "=>"; limit = 2)))
            occursin(IDENTIFIER, old) && occursin(IDENTIFIER, new) ||
                error("$file:$number: a rename names two identifiers")
            push!(renames, old => new)
        elseif length(fields) == 3 && Symbol(fields[3]) in KINDS
            push!(declared, (; scenario = fields[1], prefix = fields[2], kind = Symbol(fields[3])))
        else
            error("$file:$number: expected `rename | Old => New` or " *
                "`scenario | path prefix | value, type or path`")
        end
    end
    return (; renames, declared)
end

# `text` with every rename applied to whole identifiers; `used[i]` records a match.
function renamed(text, renames, used)
    for (i, (old, new)) in enumerate(renames)
        pattern = Regex("(?<![\\w!])\\Q" * old * "\\E(?![\\w!])")
        occursin(pattern, text) || continue
        used[i] = true
        text = replace(text, pattern => new)
    end
    return text
end

# Nodes whose value is a type name (an `eltype` row or a type met in the walk) or a
# Symbol: renames reach their values too.
const NAMED = ("Type", "DataType", "UnionAll", "Union", "TypeofBottom", "Symbol")

# REF's nodes with the renames applied to paths, type names and Symbol values.
renamed(nodes::Vector, renames, used) = [renamed(path, renames, used) =>
    (renamed(type, renames, used), type in NAMED ? renamed(value, renames, used) : value)
    for (path, (type, value)) in nodes]

# The differing paths of a scenario, with the node on each side: the working tree's
# paths in recorded order, then the paths found only at the revision.
function differences(old, new)
    before, after = Dict(old), Dict(new)
    paths = [first.(new); [path for path in first.(old) if !haskey(after, path)]]
    return [(path, get(before, path, nothing), get(after, path, nothing))
        for path in paths if get(before, path, nothing) != get(after, path, nothing)]
end

# A difference changes a value under the same type, changes the type, or exists on one
# side only.
kind((_, a, b)) = a === nothing || b === nothing ? :path : a[1] == b[1] ? :value : :type

covers(declaration, scenario, change) = declaration.scenario == scenario &&
    declaration.kind === kind(change) && startswith(first(change), declaration.prefix)

shown(x) = x === nothing ? "absent" :
    first(x[2] == "" ? x[1] : string(x[1], " ", x[2]), 150)

laziness(nodes) = (i = findfirst(n -> first(n) == "laziness", nodes);
    i === nothing ? "absent" : last(last(nodes[i])))

function compare(reference, allowed)
    (; renames, declared) = declarations(allowed)
    return at_revision(reference) do project
        recorder = own_copy(dirname(project), RECORDER)
        println(reference, ": ", recorder.own ? "its own recording program and corpus" :
            "the working tree's recording program and corpus ($reference has no $RECORDER)")
        println("working tree: its own recording program and corpus")
        directory = mktempdir()
        old = fingerprints(project, recorder.path, joinpath(directory, "reference.tsv"))
        new = fingerprints(joinpath(REPOSITORY, "test"), joinpath(REPOSITORY, RECORDER),
            joinpath(directory, "working.tsv"))
        used_renames = falses(length(renames))
        used = falses(length(declared))
        unexpected = 0
        for scenario in unique([new.order; old.order])
            before = renamed(get(old.nodes, scenario, []), renames, used_renames)
            nodes = get(new.nodes, scenario, [])
            changes = differences(before, nodes)
            counts = (laziness(before), laziness(nodes))
            lazy = counts == ("absent", "absent") ? "" : first(counts) == last(counts) ?
                "; laziness: " * last(counts) : "; laziness: $(first(counts)) → $(last(counts))"
            if isempty(changes)
                println(rpad(scenario, 16), "identical (", length(nodes), " nodes", lazy, ")")
                continue
            end
            free = filter(changes) do change
                matches = [i for (i, d) in enumerate(declared) if covers(d, scenario, change)]
                used[matches] .= true
                isempty(matches)
            end
            unexpected += length(free)
            tally = [count(c -> kind(c) === k, changes) for k in KINDS]
            println(rpad(scenario, 16), length(changes), " of ", length(nodes), " nodes differ: ",
                tally[1], " values, ", tally[2], " types, ", tally[3], " paths",
                length(free) == length(changes) ? "" : " ($(length(changes) - length(free)) declared)",
                lazy)
            patterns = [replace(first(c), r"\[[^\]]*\]" => "[*]") for c in changes]
            println("    by path: ", join(["$p $(count(==(p), patterns))"
                for p in first(unique(patterns), 12)], ", "),
                length(unique(patterns)) > 12 ? ", …" : "")
            listed = isempty(free) ? changes : free
            for k in KINDS, (path, a, b) in first(filter(c -> kind(c) === k, listed), SHOWN)
                println("    ", k, " ", path, ": ", shown(a), "  →  ", shown(b))
            end
        end
        for (i, (old, new)) in enumerate(renames)
            used_renames[i] || println("Rename that matches nothing: ", old, " => ", new)
        end
        for (i, d) in enumerate(declared)
            used[i] || println("Declared but not found: ", d.scenario, " | ", d.prefix, " | ", d.kind)
        end
        return unexpected
    end
end

function main(arguments)
    1 <= length(arguments) <= 2 ||
        error("Usage: julia --project=test test/tools/equivalence.jl REF [ALLOWED]")
    unexpected = compare(arguments[1], get(arguments, 2, nothing))
    unexpected == 0 && return 0
    println(stderr, unexpected, " undeclared differences.")
    return 1
end

abspath(PROGRAM_FILE) == (@__FILE__) && exit(main(ARGS))
