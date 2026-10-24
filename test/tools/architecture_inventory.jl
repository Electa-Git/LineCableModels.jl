# Print the live architecture inventory in the format of
# `test/quality/architecture_baseline.toml`. Run
# `julia --project=test test/tools/architecture_inventory.jl` from the repository root.
# Compare the output with the baseline to find entries to delete or lower.
# A new violation is fixed in the source. It is never added to the baseline.
# No guard ratchets the readings that follow the inventory as comments.
source = normpath(joinpath(@__DIR__, "..", "quality", "architecture.jl"))
document = Meta.parseall(read(source, String); filename = source)
# The quality items and this tool evaluate the same guard module.
definition = only(x for x in document.args
    if Meta.isexpr(x, :macrocall) && x.args[1] === Symbol("@testmodule") &&
       x.args[3] === :ArchitectureGuards)
guards = Module(:ArchitectureGuards)
Core.eval(guards, Expr(:toplevel, definition.args[end].args...))
print(Base.invokelatest(guards.render, Base.invokelatest(guards.live)))

titles = (decisions = "Private functions with two or more methods of one arity",
    named_helpers = "Tiny helpers exempt only because a test file names them",
    owners = "Exported or public Commons names and their owners outside Commons",
    mentions = "Commons names that a Commons unit test names and never calls")
readings = Base.invokelatest(guards.readings, Base.invokelatest(guards.live_tree))
for (reading, entries) in pairs(readings)
    println("\n# ", titles[reading], ": ", length(entries))
    for key in sort!(collect(keys(entries)))
        value = entries[key]
        text = value isa AbstractVector ? (isempty(value) ? "none" : join(value, ", ")) : value
        println("#   ", key, " = ", text)
    end
end
