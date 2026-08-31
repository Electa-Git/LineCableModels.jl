using LineCableModels
using CairoMakie

CairoMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))
CairoMakie.activate!()

presentation = find_presentation("live-manual", PRESENTATION_DESCRIPTORS)
cable_deck = find_deck("core-and-insulation", presentation.all_decks)
plot = cable_deck.export_figure()
output_directory = joinpath(@__DIR__, "build")
mkpath(output_directory)
output_path = joinpath(output_directory, "cable-design.svg")
Makie.save(output_path, plot.figure)
println(abspath(output_path))
