import Pkg

Pkg.Apps.develop(path = @__DIR__)

bindir = joinpath(first(DEPOT_PATH), "bin")
println("Installed the linecablemodels command in $bindir")
println("Add that directory to PATH, then run:")
println("  linecablemodels playground start")
