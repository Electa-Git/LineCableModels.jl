"""
$(TYPEDEF)

Hold one built page together with its completed detached recipe, runtime
panels, controls, and backend context.

$(TYPEDFIELDS)
"""
struct UIPlot{D <: AbstractPlotDefinition, G, F, P, W, C}
    "Completed detached plot recipe."
    render::PlotRecipe{D}
    "Page represented by this handle."
    page::G
    "Backend-built figure."
    figure::F
    "Registered runtime axes."
    panels::P
    "Interactive control objects keyed by purpose."
    controls::W
    "Active backend and status context."
    context::C
end
