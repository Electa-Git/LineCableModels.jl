using Test,LineCableModels,GLMakie
include(joinpath(@__DIR__,"../../test/support/scenarios.jl"))
using .CurrentScenarios
reference=two_conductor_results()
omega=reshape(2pi.*frequencies(reference),1,1,:)
function comparison(channels)
    r,l,g,c=channels
    LineParameters(r*R(reference).+im.*omega.*l.*L(reference),
        g*G(reference).+im.*omega.*c.*C(reference),frequencies(reference);
        # Before: an untyped details NamedTuple was accepted. Completed results
        # now carry ComputationDetails; the coordinate names keep their meaning.
        details=ComputationDetails((coordinates=["west","east"],)))
end
# Manufactured routing variants carry no author or physical-accuracy label.
plots=Makie.plot(reference,comparison((1.1,1.2,1.3,1.4)),comparison((.9,.8,.7,.6));
    series_labels=("Current routing input","Larger channels","Smaller channels"),
    # Before: result plots accepted size. Now fig_size describes the nominal layout;
    # each quantity still returns its own UIPlot (preview/plotwindow retain size).
    ydata=(R,L,G,C),xscale=:log10,fig_size=(1100,750),backend=:gl,display_plot=true,open_export=false)
@test length(plots)==4
@test all(page->length(page.axes)==4,plots)
@test all(page->haskey(page.controls,:legend),plots)
println("Inspect the four quantity pages, resize them and toggle legend entries.")
println("Press Enter to close the figures.")
readline()
GLMakie.closeall()
