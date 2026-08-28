using LineCableModels
using CairoMakie

CairoMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

state = Observable(design(12.5, 8.0))
plot = cable_figure(state)
output_directory = joinpath(@__DIR__, "build")
mkpath(output_directory)
output_path = joinpath(output_directory, "cable-design.svg")
Makie.save(output_path, plot.figure)
println(abspath(output_path))
