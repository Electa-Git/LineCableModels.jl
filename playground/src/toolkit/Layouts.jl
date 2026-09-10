"""
    WorkspacePage(title, content; eyebrow="", description=nothing, tools=nothing, fill=false)

Compose the shared engineering-page heading and content region. `fill=true`
reserves the remaining host height for a persistent viewport or `SplitPane`;
otherwise content follows its natural height. The host supplies available size,
not a second heading or padding. Construction owns no application state.
"""
struct WorkspacePage{C,D,T}
    "Visible page title."
    title::String
    "Renderable page content."
    content::C
    "Optional compact heading label."
    eyebrow::String
    "Optional description below the title."
    description::D
    "Optional heading-side content."
    tools::T
    "Whether content fills the available host height."
    fill::Bool
end

WorkspacePage(title::AbstractString, content; eyebrow::AbstractString="",
    description=nothing, tools=nothing, fill::Bool=false) =
    WorkspacePage(string(title), content, string(eyebrow), description, tools, fill)

function Bonito.jsrender(session::Session, page::WorkspacePage)
    node = DOM.section(
        DOM.header(DOM.div(
            isempty(page.eyebrow) ? nothing : DOM.span(page.eyebrow; class="lc-workspace-eyebrow"),
            DOM.h1(page.title; class="lc-workspace-title"),
            page.description === nothing ? nothing : DOM.p(page.description)),
            page.tools; class="lc-workspace-header"),
        DOM.div(page.content; class="lc-workspace-body");
        class="lc-workspace-page" * (page.fill ? " is-fill" : ""))
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, page))
end

function ComponentXRay.inspection(page::WorkspacePage)
    return ComponentXRay.ComponentInspection(page; name="WorkspacePage",
        source=toolkit_source(@__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:title, page.title),
            ComponentXRay.PropertyInspection(:eyebrow, page.eyebrow),
            ComponentXRay.PropertyInspection(:fill, page.fill)],
        css_scopes=[".lc-workspace-page", ".lc-workspace-header", ".lc-workspace-body",
            ".lc-workspace-eyebrow", ".lc-workspace-title"])
end
