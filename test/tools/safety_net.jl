# The safety net against an earlier revision. Each revision runs its own tutorials with
# its own code, and this script compares the results. Run
# `julia --project=test test/tools/safety_net.jl BASE [OTHER]`
# from the repository root. OTHER defaults to HEAD.
#
# `worktree.jl` checks out each revision in a temporary git worktree, OTHER included, and
# prepares its test environment. A fresh process of that environment runs the recorder.
# It runs each `examples/tutorial*.jl` of the revision in a fresh module, from the
# revision's `examples` directory.
#
# The recorder keeps every result as numbers or text:
# - the R, L, C and G of cable constants
# - the Z, Y and frequencies of line parameters, and the Z and Y of a modal result
# - the numbers and numeric arrays that a tutorial leaves in its globals
# - the text of every file that a run writes
#
# A failing step leaves its error in the record, and the run goes on. Each step records
# its run time. The outputs go to `.linecablemodels/checkpoint/<revision>/`: `records.jls`,
# the files of the runs, and `log.txt`. The comparison goes to
# `.linecablemodels/checkpoint/comparison.md`. For each key that differs, it gives the
# largest relative difference, the lines that changed, or the side that lacks the key or
# fails.
const REPOSITORY = dirname(dirname(@__DIR__))
const OUTPUT = joinpath(REPOSITORY, ".linecablemodels", "checkpoint")

Base.include(@__MODULE__, joinpath(@__DIR__, "worktree.jl"))
using Printf: @sprintf
using Serialization: deserialize

const RECORDER = raw"""
using LineCableModels, Serialization
root, out = ARGS
records = Dict{String, Any}()
seconds = Dict{String, Float64}()

flat(error) = first(replace(sprint(showerror, error), r"\s+" => " "), 300)
numbers(x::Number) = (kind = :numbers, size = (), values = ComplexF64[nominal(x)])
numbers(x::AbstractArray{<:Number}) =
    (kind = :numbers, size = size(x), values = ComplexF64[nominal(v) for v in vec(x)])
failure(error) = (kind = :failure, message = flat(error))

# The results of one value, by the name of its type, which both revisions share.
function record!(key, value)
    name = String(nameof(typeof(value)))
    try
        if value isa Number || value isa AbstractArray{<:Number}
            records[key] = numbers(value)
        elseif name == "CableConstants"
            for field in (:R, :L, :C, :G)
                records["$key.$field"] = numbers(getproperty(value, field))
            end
        elseif name == "LineParameters" || (occursin("Modal", name) &&
                !occursin("Problem", name) && !occursin("Formulation", name))
            records["$key.Z"] = numbers(observe(value, Z))
            records["$key.Y"] = numbers(observe(value, Y))
            name == "LineParameters" && (records["$key.frequencies"] = numbers(frequencies(value)))
        end
    catch error
        records[key] = failure(error)
    end
    return nothing
end

function record_globals!(prefix, m)
    for name in sort!(names(m; all = true))
        text = String(name)
        (startswith(text, "#") || name in (:eval, :include) || !isdefined(m, name)) && continue
        value = getfield(m, name)
        value isa Union{Module, Function, Type} && continue
        record!("$prefix.$text", value)
    end
end

# Every file under `root`, outside `.git`, with its modification time.
function listing(root)
    found = Dict{String, Float64}()
    for (directory, folders, files) in walkdir(root)
        filter!(!=(".git"), folders)
        for file in files
            path = joinpath(directory, file)
            found[path] = mtime(path)
        end
    end
    return found
end

# The files that a step wrote: new, or changed since `before`. Each is copied to the
# output directory and its text is recorded.
function record_files!(prefix, root, before)
    for (path, time) in listing(root)
        get(before, path, -Inf) == time && continue
        relative = relpath(path, root)
        saved = joinpath(out, "files", prefix, relative)
        mkpath(dirname(saved))
        cp(path, saved; force = true)
        text = read(path, String)
        records["$prefix.file:$relative"] = isvalid(text) ? (kind = :text, text) :
                                            (kind = :binary, bytes = length(text))
    end
end

function timed!(f, key)
    started = time()
    try
        f()
    catch error
        records[key] = failure(error)
    end
    seconds[key] = time() - started
end

examples = joinpath(root, "examples")
for file in sort!(filter(f -> occursin(r"^tutorial\d+\.jl$", f), readdir(examples)))
    key = first(splitext(file))
    m = Module(Symbol(key))
    before = listing(root)
    timed!(key) do
        cd(() -> Base.include(m, joinpath(examples, file)), examples)
    end
    # The tutorial's globals belong to a newer world than this loop.
    Base.invokelatest(record_globals!, key, m)
    record_files!(key, root, before)
    println("recorded ", key, " in ", round(seconds[key]; digits = 1), " s")
end

serialize(joinpath(out, "records.jls"), (; records, seconds))
"""

# Records `reference` into `OUTPUT/<name>` and returns the wall time [s] of the whole run.
function record(reference, name)
    out = joinpath(OUTPUT, name)
    rm(out; force = true, recursive = true)
    mkpath(out)
    started = time()
    at_revision(reference) do project
        open(joinpath(out, "log.txt"), "w") do log
            run(pipeline(ignorestatus(julia(project, RECORDER, dirname(project), out));
                stdout = log, stderr = log))
        end
    end
    return time() - started
end

# The largest difference of `new` from `old`, relative to the largest entry of `old`.
function relative(old, new)
    old.size == new.size || return "sizes $(old.size) → $(new.size)"
    old.values == new.values && return nothing
    scale = maximum(abs, old.values; init = 0.0)
    difference = maximum(abs.(new.values .- old.values); init = 0.0)
    return iszero(scale) ? @sprintf("%.3g absolute", difference) :
           @sprintf("%.3g", difference / scale)
end

function textual(old, new)
    old.text == new.text && return nothing
    a, b = split(old.text, '\n'), split(new.text, '\n')
    changed = [i for i in 1:max(length(a), length(b))
               if get(a, i, nothing) != get(b, i, nothing)]
    shown = join(["line $i: `$(first(get(a, i, ""), 160))` → `$(first(get(b, i, ""), 160))`"
                  for i in first(changed, 3)], "; ")
    return "$(length(changed)) of $(max(length(a), length(b))) lines differ. $shown"
end

describe(entry) = entry === nothing ? "absent" :
    entry.kind === :failure ? "fails: $(entry.message)" : String(entry.kind)

function compare(names, times)
    old, new = (deserialize(joinpath(OUTPUT, name, "records.jls")) for name in names)
    lines = ["# Safety net: $(names[1]) against $(names[2])", "",
        "Run time: $(join(["$(n) $(round(Int, t)) s" for (n, t) in zip(names, times)], ", ")).",
        "", "| Key | Result |", "| --- | --- |"]
    counts = Dict{String, Int}()
    for key in sort!(collect(union(keys(old.records), keys(new.records))))
        a, b = get(old.records, key, nothing), get(new.records, key, nothing)
        result = if a === nothing || b === nothing || a.kind !== b.kind ||
                    a.kind in (:failure, :binary)
            a == b ? nothing : "$(describe(a)) → $(describe(b))"
        else
            a.kind === :numbers ? relative(a, b) : textual(a, b)
        end
        counts[result === nothing ? "identical" : "differs"] =
            get(counts, result === nothing ? "identical" : "differs", 0) + 1
        result === nothing || push!(lines, "| `$key` | $(replace(result, '|' => '/')) |")
    end
    push!(lines, "", "Identical keys: $(get(counts, "identical", 0)). Differing keys: " *
        "$(get(counts, "differs", 0)).", "", "| Step | $(names[1]) [s] | $(names[2]) [s] |",
        "| --- | --- | --- |")
    for key in sort!(collect(union(keys(old.seconds), keys(new.seconds))))
        a, b = (haskey(s, key) ? string(round(s[key]; digits = 1)) : "—"
                for s in (old.seconds, new.seconds))
        push!(lines, "| `$key` | $a | $b |")
    end
    write(joinpath(OUTPUT, "comparison.md"), join(lines, '\n') * '\n')
    println("Identical keys: ", get(counts, "identical", 0), ", differing keys: ",
        get(counts, "differs", 0), ". See ", relpath(joinpath(OUTPUT, "comparison.md"),
            REPOSITORY), ".")
end

function main(arguments)
    1 <= length(arguments) <= 2 ||
        error("Usage: julia --project=test test/tools/safety_net.jl BASE [OTHER]")
    references = (arguments[1], get(arguments, 2, "HEAD"))
    names = Tuple(output(git("rev-parse", "--short=8", r * "^{commit}"))[1:(end - 1)]
                  for r in references)
    times = [record(reference, name) for (reference, name) in zip(references, names)]
    compare(names, times)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS)
