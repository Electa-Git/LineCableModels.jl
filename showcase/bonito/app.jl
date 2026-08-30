using Bonito
using Markdown

include(joinpath(@__DIR__, "PowerImpedanceDiagramExt.jl"))
include(joinpath(@__DIR__, "PageAuthoring.jl"))
using .PageAuthoring

const LATO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/lato-font/3.0.0/css/lato-font.min.css"
const JULIA_MONO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/juliamono/0.050/juliamono.min.css"
const DOCUMENTATION_LOGO = joinpath(
    @__DIR__, "..", "..", "docs", "src", "assets", "logo.svg")
const PAGE_DIRECTORY = joinpath(@__DIR__, "pages")
const REQUIRED_PAGE_FIELDS = (:id, :group, :title, :order, :render, :class, :build)

function validate_page_descriptor(page, source::AbstractString)
    page isa NamedTuple || throw(ArgumentError(
        "page file $(basename(source)) must return a named tuple"
    ))
    missing = filter(field -> !hasproperty(page, field), REQUIRED_PAGE_FIELDS)
    isempty(missing) || throw(ArgumentError(
        "page file $(basename(source)) is missing fields: $(join(missing, ", "))"
    ))
    page.id isa AbstractString && !isempty(page.id) || throw(ArgumentError(
        "page id in $(basename(source)) must be a non-empty string"
    ))
    occursin(r"^[a-z0-9]+(?:-[a-z0-9]+)*$", page.id) || throw(ArgumentError(
        "page id $(repr(page.id)) must contain lowercase words separated by hyphens"
    ))
    page.group isa AbstractString && !isempty(page.group) || throw(ArgumentError(
        "page group in $(basename(source)) must be a non-empty string"
    ))
    page.title isa AbstractString && !isempty(page.title) || throw(ArgumentError(
        "page title in $(basename(source)) must be a non-empty string"
    ))
    page.order isa Real && isfinite(page.order) || throw(ArgumentError(
        "page order in $(basename(source)) must be finite"
    ))
    page.render isa Bool || throw(ArgumentError(
        "page render flag in $(basename(source)) must be true or false"
    ))
    page.class isa AbstractString || throw(ArgumentError(
        "page class in $(basename(source)) must be a string"
    ))
    page.build isa Function || throw(ArgumentError(
        "page build entry in $(basename(source)) must be a function"
    ))
    return merge(page, (; source = abspath(source)))
end

function load_page_descriptors(directory::AbstractString = PAGE_DIRECTORY)
    isdir(directory) || throw(ArgumentError("page directory does not exist: $directory"))
    sources = sort(filter(path -> endswith(path, ".jl"), readdir(directory; join = true)))
    isempty(sources) && throw(ArgumentError("no Julia page files found in $directory"))
    pages = map(sources) do source
        validate_page_descriptor(Base.include(@__MODULE__, source), source)
    end
    ids = getproperty.(pages, :id)
    duplicates = unique(filter(id -> count(==(id), ids) > 1, ids))
    isempty(duplicates) || throw(ArgumentError(
        "duplicate page ids: $(join(duplicates, ", "))"
    ))
    return sort(pages; by = page -> (page.order, basename(page.source)))
end

function select_rendered_pages(pages)
    rendered = filter(page -> page.render, pages)
    isempty(rendered) && throw(ArgumentError(
        "no showcase pages are enabled; set `render = true` in at least one pages/*.jl file"
    ))
    return rendered
end

function find_page(id::AbstractString, pages)
    return only(filter(page -> page.id == id, pages))
end

const PAGE_DESCRIPTORS = load_page_descriptors()
const RENDERED_PAGES = select_rendered_pages(PAGE_DESCRIPTORS)

function app_assets()
    return (
        lato_css = Asset(LATO_STYLESHEET),
        julia_mono_css = Asset(JULIA_MONO_STYLESHEET),
        app_css = read(joinpath(@__DIR__, "assets", "theme.css"), String),
        logo = Asset(DOCUMENTATION_LOGO)
    )
end

page_route(page, index::Integer) = index == 1 ? "/" : "/$(page.id)"
page_href(page, index::Integer) = index == 1 ? "./" : "./$(page.id)"

function page_index(page, pages)
    index = findfirst(candidate -> candidate.id == page.id, pages)
    isnothing(index) && throw(ArgumentError("page $(repr(page.id)) is not rendered"))
    return index
end

function navigation(assets, pages, current_index::Integer)
    groups = unique(page.group for page in pages)
    group_nodes = map(groups) do group
        entries = map(
            filter(pair -> pair[2].group == group,
            collect(enumerate(pages)))
        ) do (index, page)
            active = index == current_index
            DOM.li(
                DOM.a(
                page.title;
                href = page_href(page, index),
                class = active ? "lc-nav-entry is-active" : "lc-nav-entry",
                dataPageId = page.id,
                ariaCurrent = active ? "page" : "false"
            )
            )
        end
        DOM.div(
            DOM.p(group; class = "lc-nav-group-title"),
            DOM.ul(entries...; class = "lc-nav-list");
            class = "lc-nav-group",
            dataNavGroup = group
        )
    end

    return DOM.aside(
        DOM.div(
            DOM.img(; src = assets.logo, alt = "", class = "lc-brand-logo"),
            DOM.div(
                DOM.strong("LineCableModels.jl"),
                DOM.span("Live manual");
                class = "lc-brand-copy"
            );
            class = "lc-sidebar-header"
        ),
        DOM.div(
            DOM.input(;
                id = "page-search",
                type = "search",
                placeholder = "Search pages",
                autocomplete = "off",
                ariaLabel = "Search pages"
            );
            class = "lc-search"
        ),
        DOM.nav(
            group_nodes...;
            id = "page-navigation",
            class = "lc-sidebar-nav",
            ariaLabel = "Manual pages"
        ),
        DOM.div(
            DOM.span("Bonito showcase"),
            DOM.span("Bonito · WGLMakie");
            class = "lc-sidebar-footer"
        );
        id = "lc-sidebar",
        class = "lc-sidebar",
        ariaLabel = "Page navigation"
    )
end

function page_control(label; id, target, target_index, title, aria_label)
    if isnothing(target)
        return DOM.span(
            label;
            id,
            class = "lc-navbar-button is-disabled",
            title,
            ariaLabel = aria_label,
            ariaDisabled = "true"
        )
    end
    return DOM.a(
        label;
        id,
        href = page_href(target, target_index),
        class = "lc-navbar-button",
        title,
        ariaLabel = aria_label
    )
end

function navbar(pages, current_index::Integer; preparing::Bool = false)
    current = pages[current_index]
    previous = current_index == 1 ? nothing : pages[current_index - 1]
    next = current_index == length(pages) ? nothing : pages[current_index + 1]
    return DOM.header(
        DOM.button(
            "☰";
            id = "lc-sidebar-toggle",
            class = "lc-sidebar-toggle",
            type = "button",
            ariaLabel = "Toggle page navigation",
            ariaControls = "lc-sidebar",
            ariaExpanded = "true"
        ),
        DOM.nav(
            DOM.span(
                current.group;
                id = "lc-current-group",
                class = "lc-breadcrumb-parent"
            ),
            DOM.span("›"; class = "lc-breadcrumb-separator", ariaHidden = "true"),
            DOM.span(current.title; id = "lc-current-title");
            class = "lc-breadcrumb",
            ariaLabel = "Current page"
        ),
        DOM.div(
            page_control(
                "‹";
                id = "lc-previous-page",
                target = previous,
                target_index = current_index - 1,
                title = "Previous page",
                aria_label = "Previous page"
            ),
            DOM.output(
                "$(current_index) / $(length(pages))";
                id = "lc-page-position",
                class = "lc-page-position",
                ariaLive = "polite"
            ),
            page_control(
                "›";
                id = "lc-next-page",
                target = next,
                target_index = current_index + 1,
                title = "Next page",
                aria_label = "Next page"
            ),
            DOM.button(
                "⛶";
                id = "lc-fullscreen-toggle",
                class = "lc-navbar-button",
                type = "button",
                title = "Presentation mode",
                ariaLabel = "Enter presentation mode",
                ariaPressed = "false"
            ),
            DOM.button(
                "⎙";
                id = "lc-print",
                class = "lc-navbar-button",
                type = "button",
                title = "Print or save as PDF",
                ariaLabel = "Print or save as PDF"
            );
            class = "lc-navbar-actions",
            ariaLabel = "Page controls"
        ),
        DOM.span(
            preparing ? "Preparing application case" : "Interactive showcase";
            class = "lc-navbar-status"
        );
        class = "lc-navbar"
    )
end

function documenter_shell(
        assets,
        pages,
        current_index::Integer,
        manual_root;
        preparing::Bool = false
)
    main = DOM.main(
        navbar(pages, current_index; preparing),
        DOM.div(manual_root; class = "lc-page-frame");
        class = "lc-main"
    )
    return DOM.div(
        assets.lato_css,
        assets.julia_mono_css,
        DOM.style(assets.app_css),
        navigation(assets, pages, current_index),
        main;
        id = "lc-live-docs",
        class = preparing ? "lc-documenter is-preparing" : "lc-documenter",
        ariaBusy = string(preparing)
    )
end

function render_page(session::Session, page)
    result = page.build(session)
    result isa NamedTuple && hasproperty(result, :body) || throw(ArgumentError(
        "page $(repr(page.id)) build function must return a named tuple with a `body` field"
    ))
    body = result.body isa Tuple ? result.body : (result.body,)
    attributes = (;
        id = "page-$(page.id)",
        class = "lc-page is-active $(page.class)",
        dataPageId = page.id,
        dataNavGroup = page.group,
        dataNavTitle = page.title,
        ariaHidden = "false"
    )
    return DOM.section(body...; attributes...)
end

function preparation_page(session::Session, page)
    body = slide(
        page.title,
        article(md"""
        This page is waiting for its process-wide numerical model. The server remains responsive and this status is updated from the actual Julia preparation phases.
        """),
        preparation_status(session);
        lede = md"The application case is being prepared once and will be reused by every browser session."
    )
    return DOM.section(
        body...;
        id = "page-$(page.id)",
        class = "lc-page is-active $(page.class)",
        dataPageId = page.id,
        dataNavGroup = page.group,
        dataNavTitle = page.title,
        ariaHidden = "false"
    )
end

function bind_preparation_reload(session::Session)
    ready = map(
        state -> state.phase === :ready,
        session,
        preparation_state()
    )
    Bonito.onjs(
        session,
        ready,
        js"""
isReady => {
    if (isReady) {
        window.location.reload();
    }
}
"""
    )
    notify(ready)
    return ready
end

function bind_shell_behavior(session::Session, root)
    Bonito.onload(
        session, root,
        js"""
(root) => {
    if (root.__lineCableManual) {
        return;
    }

    root.__lineCableManual = true;
    const manualRoot = root.querySelector("#linecable-app");
    const search = root.querySelector("#page-search");
    const sidebarToggle = root.querySelector("#lc-sidebar-toggle");
    const previousLink = root.querySelector("#lc-previous-page");
    const nextLink = root.querySelector("#lc-next-page");
    const fullscreenToggle = root.querySelector("#lc-fullscreen-toggle");
    const printButton = root.querySelector("#lc-print");
    const navEntries = Array.from(root.querySelectorAll(".lc-nav-entry"));
    const navGroups = Array.from(root.querySelectorAll(".lc-nav-group"));
    const mobileSidebar = window.matchMedia("(max-width: 56rem)");
    const visualViewport = window.visualViewport;
    const presentationModeKey = "linecable:presentation-mode";
    let touchStart = null;
    let firstViewportFrame = 0;
    let secondViewportFrame = 0;
    let viewportTimer = 0;
    let resolutionQuery = null;

    const applyViewport = () => {
        const width = document.documentElement.clientWidth;
        const height = document.documentElement.clientHeight;
        if (width <= 0 || height <= 0) {
            return;
        }
        root.style.setProperty("--lc-viewport-width", `${width}px`);
        root.style.setProperty("--lc-viewport-height", `${height}px`);
        root.dataset.viewportSize = `${width}x${height}`;
        root.dataset.pixelRatio = String(window.devicePixelRatio || 1);
    };

    const settleViewport = () => {
        applyViewport();
        cancelAnimationFrame(firstViewportFrame);
        cancelAnimationFrame(secondViewportFrame);
        clearTimeout(viewportTimer);
        firstViewportFrame = requestAnimationFrame(() => {
            secondViewportFrame = requestAnimationFrame(applyViewport);
        });
        viewportTimer = window.setTimeout(applyViewport, 240);
    };

    const bindResolutionQuery = () => {
        if (resolutionQuery) {
            resolutionQuery.removeEventListener("change", handleResolutionChange);
        }
        resolutionQuery = window.matchMedia(
            `(resolution: ${window.devicePixelRatio || 1}dppx)`
        );
        resolutionQuery.addEventListener("change", handleResolutionChange);
    };

    function handleResolutionChange() {
        bindResolutionQuery();
        settleViewport();
    }

    const releaseViewport = () => {
        window.removeEventListener("resize", settleViewport);
        visualViewport?.removeEventListener("resize", settleViewport);
        resolutionQuery?.removeEventListener("change", handleResolutionChange);
        cancelAnimationFrame(firstViewportFrame);
        cancelAnimationFrame(secondViewportFrame);
        clearTimeout(viewportTimer);
    };

    window.addEventListener("resize", settleViewport, { passive: true });
    visualViewport?.addEventListener("resize", settleViewport, { passive: true });
    window.addEventListener("pagehide", releaseViewport, { once: true });
    bindResolutionQuery();
    settleViewport();

    const navigate = (link) => {
        if (link instanceof HTMLAnchorElement && link.href) {
            window.location.assign(link.href);
        }
    };

    search.addEventListener("input", () => {
        const query = search.value.trim().toLocaleLowerCase();
        navEntries.forEach((entry) => {
            const group = entry.closest(".lc-nav-group").dataset.navGroup;
            const searchable = `${group} ${entry.textContent}`.toLocaleLowerCase();
            entry.closest("li").hidden = query.length > 0 && !searchable.includes(query);
        });
        navGroups.forEach((group) => {
            const visibleEntry = Array.from(group.querySelectorAll("li")).some((item) => !item.hidden);
            group.hidden = !visibleEntry;
        });
    });
    search.addEventListener("keydown", (event) => {
        if (event.key === "Escape") {
            search.value = "";
            search.dispatchEvent(new Event("input"));
            search.blur();
        }
        event.stopPropagation();
    });

    const interactiveTarget = (target) => Boolean(target.closest(
        "input, button, a, select, textarea, canvas, [contenteditable='true'], [role='slider']"
    ));
    root.addEventListener("keydown", (event) => {
        if ((event.ctrlKey || event.metaKey) && event.key === "/") {
            event.preventDefault();
            search.focus();
            return;
        }
        if (interactiveTarget(event.target)) {
            return;
        }
        if (event.key === "ArrowLeft" || event.key === "PageUp") {
            event.preventDefault();
            navigate(previousLink);
        } else if (event.key === "ArrowRight" || event.key === "PageDown") {
            event.preventDefault();
            navigate(nextLink);
        }
    });

    sidebarToggle.addEventListener("click", () => {
        if (mobileSidebar.matches) {
            const open = root.classList.toggle("is-sidebar-open");
            sidebarToggle.setAttribute("aria-expanded", String(open));
        } else {
            const collapsed = root.classList.toggle("is-sidebar-collapsed");
            sidebarToggle.setAttribute("aria-expanded", String(!collapsed));
        }
    });

    manualRoot.addEventListener("pointerdown", (event) => {
        if (event.pointerType === "touch" && !interactiveTarget(event.target)) {
            touchStart = { x: event.clientX, y: event.clientY };
        }
        if (!interactiveTarget(event.target)) {
            manualRoot.focus({ preventScroll: true });
        }
    });
    manualRoot.addEventListener("pointerup", (event) => {
        if (!touchStart || event.pointerType !== "touch") {
            touchStart = null;
            return;
        }
        const dx = event.clientX - touchStart.x;
        const dy = event.clientY - touchStart.y;
        touchStart = null;
        if (Math.abs(dx) >= 60 && Math.abs(dx) > 1.25 * Math.abs(dy)) {
            navigate(dx < 0 ? nextLink : previousLink);
        }
    });
    manualRoot.addEventListener("pointercancel", () => {
        touchStart = null;
    });

    const readPresentationMode = () => {
        try {
            return window.sessionStorage.getItem(presentationModeKey) === "true";
        } catch (error) {
            console.warn("Presentation-mode state is unavailable.", error);
            return false;
        }
    };
    const rememberPresentationMode = (active) => {
        try {
            if (active) {
                window.sessionStorage.setItem(presentationModeKey, "true");
            } else {
                window.sessionStorage.removeItem(presentationModeKey);
            }
        } catch (error) {
            console.warn("Presentation-mode state could not be saved.", error);
        }
    };
    const setFocusMode = (active, persist = true) => {
        root.classList.toggle("is-focus-mode", active);
        fullscreenToggle.setAttribute("aria-pressed", String(active));
        fullscreenToggle.setAttribute(
            "aria-label",
            active ? "Exit presentation mode" : "Enter presentation mode"
        );
        fullscreenToggle.title = active ? "Exit presentation mode" : "Presentation mode";
        if (persist) {
            rememberPresentationMode(active);
        }
        settleViewport();
    };
    setFocusMode(readPresentationMode(), false);

    fullscreenToggle.addEventListener("click", async () => {
        if (root.classList.contains("is-focus-mode")) {
            setFocusMode(false);
            if (document.fullscreenElement === root) {
                await document.exitFullscreen();
            }
            return;
        }
        setFocusMode(true);
        if (!root.requestFullscreen) {
            return;
        }
        try {
            await root.requestFullscreen();
        } catch (error) {
            console.warn("Browser fullscreen is unavailable; presentation mode remains active.", error);
        }
    });
    document.addEventListener("fullscreenchange", () => {
        if (document.fullscreenElement === root) {
            setFocusMode(true);
        } else {
            settleViewport();
        }
    });

    printButton.addEventListener("click", () => window.print());
    root.dataset.manualState = "ready";
}
"""
    )
    return root
end

function manual(session::Session, assets, page; pages = RENDERED_PAGES)
    pages = select_rendered_pages(pages)
    current_index = page_index(page, pages)
    preparing = hasproperty(page, :ready) && !page.ready()
    rendered_page = preparing ? preparation_page(session, page) :
                    render_page(session, page)
    page_nodes = DOM.div(
        rendered_page;
        id = "linecable-pages",
        class = "lc-pages"
    )
    manual_root = DOM.div(
        page_nodes;
        id = "linecable-app",
        class = "lc-manual",
        tabindex = "0",
        ariaLabel = "LineCableModels live manual pages"
    )
    root = documenter_shell(assets, pages, current_index, manual_root; preparing)
    preparing && bind_preparation_reload(session)
    return bind_shell_behavior(session, root)
end

function cable_app(
        assets = app_assets();
        pages = RENDERED_PAGES,
        page = first(select_rendered_pages(pages))
)
    pages = select_rendered_pages(pages)
    return App(
        session -> manual(session, assets, page; pages);
        title = "LineCableModels live manual",
        indicator = ConnectionIndicator(; size = 8, top = "20px", right = "16px")
    )
end

function cable_routes(
        assets = app_assets();
        pages = RENDERED_PAGES
)
    pages = select_rendered_pages(pages)
    pairs = map(enumerate(pages)) do (index, page)
        page_route(page, index) => cable_app(assets; pages, page)
    end
    return Routes(pairs...)
end

function start_page_preparations!(pages = RENDERED_PAGES)
    pages = select_rendered_pages(pages)
    return map(
        page -> hasproperty(page, :prepare) ? page.prepare() : nothing,
        pages
    )
end
