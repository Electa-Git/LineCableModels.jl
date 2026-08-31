const PLAYGROUND_ROOT = normpath(joinpath(@__DIR__, ".."))
const REPOSITORY_ROOT = normpath(joinpath(PLAYGROUND_ROOT, ".."))
const THEME_PATH = joinpath(PLAYGROUND_ROOT, "assets", "theme.css")
const LOGO_PATH = joinpath(REPOSITORY_ROOT, "docs", "src", "assets", "logo.svg")

function navigation(logo)
    return DOM.aside(
        DOM.div(
            DOM.img(; src = logo, alt = "", class = "lc-brand-logo"),
            DOM.div(
                DOM.strong("LineCableModels.jl"),
                DOM.span("Playground");
                class = "lc-brand-copy"
            );
            class = "lc-sidebar-header"
        ),
        DOM.nav(
            DOM.div(
                DOM.p("Playground"; class = "lc-nav-group-title"),
                DOM.ul(
                    DOM.li(
                        DOM.a(
                            "Home";
                            href = "./",
                            class = "lc-nav-entry is-active",
                            ariaCurrent = "page"
                        )
                    );
                    class = "lc-nav-list"
                );
                class = "lc-nav-group"
            );
            class = "lc-sidebar-nav",
            ariaLabel = "Playground pages"
        ),
        DOM.div(
            DOM.span("LineCableModels playground"),
            DOM.span("Bonito · NATS · Quarto");
            class = "lc-sidebar-footer"
        );
        class = "lc-sidebar",
        ariaLabel = "Page navigation"
    )
end

function navbar()
    return DOM.header(
        DOM.nav(
            DOM.span("Playground"; class = "lc-breadcrumb-parent"),
            DOM.span("›"; class = "lc-breadcrumb-separator", ariaHidden = "true"),
            DOM.span("Home"; class = "lc-current-title");
            class = "lc-breadcrumb",
            ariaLabel = "Current page"
        ),
        DOM.span("Static publisher shell"; class = "lc-navbar-status");
        class = "lc-navbar"
    )
end

function foundation_card(title, role, copy; accent)
    return DOM.article(
        DOM.div(
            DOM.span(role; class = "lc-card-role"),
            DOM.span("declared"; class = "lc-card-state");
            class = "lc-card-heading"
        ),
        DOM.h2(title),
        DOM.p(copy);
        class = "lc-foundation-card",
        style = "--lc-card-accent: $accent"
    )
end

function home_page()
    return DOM.section(
        DOM.header(
            DOM.h1("LineCableModels playground"),
            DOM.p(
                "A clean publishing surface for live engineering material, with authoring, presentation and numerical execution kept deliberately separate."
            );
            class = "lc-page-heading"
        ),
        DOM.div(
            DOM.section(
                DOM.p("NEW FOUNDATION"; class = "lc-eyebrow"),
                DOM.h2("The web application publishes. The engine computes."),
                DOM.p(
                    "This home is intentionally static. It establishes the visual shell without loading presentation content, starting numerical workers or hiding execution inside the browser-facing process."
                ),
                DOM.p(
                    "Future collections will be authored independently and dispatched through an explicit broker boundary."
                );
                class = "lc-introduction-panel"
            ),
            DOM.section(
                foundation_card(
                    "Quarto",
                    "AUTHORING",
                    "Owns structured source material and reproducible document generation.";
                    accent = "#9558b2"
                ),
                foundation_card(
                    "Bonito",
                    "PUBLISHER",
                    "Owns the Julia-native web shell, sessions and presentation-facing interaction.";
                    accent = "#4eb5de"
                ),
                foundation_card(
                    "NATS",
                    "BROKER",
                    "Will route typed requests, progress and results to independent engine workers.";
                    accent = "#1abc9c"
                );
                class = "lc-foundation-grid"
            );
            class = "lc-home-canvas"
        );
        class = "lc-page"
    )
end

function playground_home(logo)
    return DOM.div(
        navigation(logo),
        DOM.main(
            navbar(),
            DOM.div(home_page(); class = "lc-page-frame");
            class = "lc-main"
        );
        class = "lc-documenter"
    )
end

function playground_app()
    stylesheet = read(THEME_PATH, String)
    logo = Asset(LOGO_PATH)
    return App(
        _session -> DOM.div(DOM.style(stylesheet), playground_home(logo));
        title = "LineCableModels playground",
        indicator = ConnectionIndicator(; size = 8, top = "20px", right = "14px")
    )
end
