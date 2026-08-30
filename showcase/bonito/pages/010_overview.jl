module OverviewPage

using Bonito
using Markdown
using ..PageAuthoring

function page_body(session::Session)
    introduction = article(
        md"""
        # Live technical manual

        This showcase behaves like a live documentation page, with navigation and fullscreen controls owned by the application shell.

        > **Interactive example**
        >
        > The cable on the next page is rebuilt from ordinary LineCableModels constructors whenever either native Bonito slider changes.

        """
    )
    boundary = article(
        md"""

        ## Showcase boundary

        - The application shell mirrors the current Documenter dark theme.
        - The WGLMakie figure is a small showcase-local renderer.
        - No PlotBuilder integration or package API adapter is introduced.

        ```text
        CableDesign
        └── Stack
            ├── Group(:core, Conductor.Solid(...))
            └── Insulator.Shell(...)
        ```
        """;
        hint = "Use the sidebar, top-bar controls, arrow keys, or swipe to move between pages."
    )
    return (introduction, preparation_status(session), boundary)
end

function build(session::Session)
    return (; body = page_body(session))
end

const PAGE = page_descriptor(
    id = "overview",
    group = "Showcase",
    title = "Live technical manual",
    order = 10,
    render = true,
    class = "lc-overview-slide",
    build = build
)

end

OverviewPage.PAGE
