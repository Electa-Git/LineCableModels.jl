using Test,LineCableModels,GLMakie,Measurements
include(joinpath(@__DIR__,"../../test/support/scenarios.jl"))
using .CurrentScenarios
smoke=lowercase(get(ENV,"LINECABLEMODELS_GL_GALLERY_SMOKE","false"))=="true"
# Synthetic completed storage exercises presentation without a new sampling run.
result=cable_monte_carlo_result()
views=[verb(result,R;backend=:gl,display_plot=!smoke,open_export=false)
    for verb in (Makie.hist,Makie.stairs,Makie.ecdfplot,Makie.lines,Makie.qqplot)]
@test all(view->view isa UIPlot,views)
@test all(view->length(view.axes)==1,views)
println("Built current histogram, density, empirical/model CDF and QQ views.")
if !smoke
    println("Press Enter to close the figures.")
    readline()
end
GLMakie.closeall()
