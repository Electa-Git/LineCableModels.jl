using Test,LineCableModels,GLMakie
include(joinpath(@__DIR__,"../../test/support/scenarios.jl"))
using .CurrentScenarios
smoke=lowercase(get(ENV,"LINECABLEMODELS_GL_GALLERY_SMOKE","false"))=="true"
design=coaxial_design()
space=Gridspace{CableConstantsProblem}(t->CableConstantsProblem(design;temperature=t),
    (Grid((20.0,60.0),AbsoluteError(1.0)),))
result=compute(ParametricProblem(space),MonteCarlo(CableConstantsFormulation();
    trials=128,seed=2027,distribution=:uniform,return_samples=true,return_histograms=true))
# Real samples are retained once; the native chart calls consume them.
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
