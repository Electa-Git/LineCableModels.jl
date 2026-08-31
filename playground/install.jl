import Pkg

Pkg.Apps.develop(path = @__DIR__)

bindir = joinpath(first(DEPOT_PATH), "bin")
println("Installed the linecablemodels command in $bindir")
configured_quarto = get(ENV, "QUARTO_PATH", "")
if isempty(configured_quarto) && isnothing(Sys.which("quarto"))
    println("Quarto CLI was not found. Run ./playground/bootstrap.sh before building.")
end
println("Add $bindir to PATH, then run:")
println("  linecablemodels playground start")
