# Compare the architecture baseline with its version at a git revision. Run
# `julia test/tools/baseline_ratchet.jl REF`. A key added or a count raised since
# `REF` fails. Deleted keys and lowered counts pass. The check passes when `REF`
# has no baseline file.
using TOML

const BASELINE = "test/quality/architecture_baseline.toml"
const REPOSITORY = normpath(joinpath(@__DIR__, "..", ".."))

git(arguments::AbstractString...) = Cmd(["git", "-C", REPOSITORY, arguments...])

# Nested tables become keys of the form "table | key | counter".
function entries(document, prefix = "")
    flat = Dict{String, Any}()
    for (key, value) in document
        name = isempty(prefix) ? key : string(prefix, " | ", key)
        value isa AbstractDict ? merge!(flat, entries(value, name)) : (flat[name] = value)
    end
    return flat
end

length(ARGS) == 1 || error("Usage: julia test/tools/baseline_ratchet.jl REF")
reference = only(ARGS)
success(pipeline(git("rev-parse", "--verify", "--quiet", reference * "^{commit}");
    stdout = devnull)) || error("Unknown git revision: $reference")
if !success(pipeline(git("cat-file", "-e", "$reference:$BASELINE"); stderr = devnull))
    println("No architecture baseline at $reference. Nothing to compare.")
    exit(0)
end
before = entries(TOML.parse(read(git("show", "$reference:$BASELINE"), String)))
after = entries(TOML.parsefile(joinpath(REPOSITORY, BASELINE)))
grown = sort!([haskey(before, key) ?
    string(key, ": ", value, " (", before[key], " at ", reference, ")") :
    string(key, ": ", value, " (absent at ", reference, ")")
    for (key, value) in after if !haskey(before, key) || value > before[key]])
if !isempty(grown)
    println(stderr, "The architecture baseline may only shrink. Entries added or raised since ",
        reference, ":")
    foreach(line -> println(stderr, "  ", line), grown)
    exit(1)
end
println("Architecture baseline: ", length(after), " entries, none added or raised since ",
    reference, ".")
