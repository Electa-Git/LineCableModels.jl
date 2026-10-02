# Print the live architecture inventory in the format of
# `test/quality/architecture_baseline.toml`. Run
# `julia --project=test test/tools/architecture_inventory.jl` from the repository root.
# Compare the output with the baseline to find entries to delete or lower.
# A new violation is fixed in the source. It is never added to the baseline.
source = normpath(joinpath(@__DIR__, "..", "quality", "architecture.jl"))
document = Meta.parseall(read(source, String); filename = source)
# The quality items and this tool evaluate the same guard module.
definition = only(x for x in document.args
    if Meta.isexpr(x, :macrocall) && x.args[1] === Symbol("@testmodule") &&
       x.args[3] === :ArchitectureGuards)
guards = Module(:ArchitectureGuards)
Core.eval(guards, Expr(:toplevel, definition.args[end].args...))
print(Base.invokelatest(guards.render, Base.invokelatest(guards.live)))
