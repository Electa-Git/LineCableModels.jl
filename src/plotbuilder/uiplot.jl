"""
    UIPlot

Hold a completed renderer-independent plot recipe together with one built figure,
its panels, controls, and backend context. Line-parameter plotting returns a
`Vector{UIPlot}`. Previews and statistical plots return one `UIPlot`.
"""
struct UIPlot{S <: AbstractPlotDefinition, F, P, W, C}
    "Completed renderer-independent plot recipe."
    render::PlotRecipe{S}
    "Page represented by this handle."
    page::PageSpec
    "Backend-built figure."
    figure::F
    "Built axes or panels."
    panels::P
    "Interactive control objects keyed by purpose."
    controls::W
    "Active backend and status context."
    context::C
end
