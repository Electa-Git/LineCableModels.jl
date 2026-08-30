using Bonito
using Markdown

include(joinpath(@__DIR__, "PowerImpedanceDiagramExt.jl"))
include(joinpath(@__DIR__, "PageAuthoring.jl"))
using .PageAuthoring

const LATO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/lato-font/3.0.0/css/lato-font.min.css"
const JULIA_MONO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/juliamono/0.050/juliamono.min.css"
const DOCUMENTATION_LOGO = joinpath(
    @__DIR__, "..", "..", "docs", "src", "assets", "logo.svg")
const DECK_DIRECTORY = joinpath(@__DIR__, "pages")
const REQUIRED_DECK_FIELDS = (:id, :group, :title, :order, :render, :class, :setup, :pages)
const REQUIRED_DECK_PAGE_FIELDS = (:id, :title, :render, :class, :build)
const SLUG_PATTERN = r"^[a-z0-9]+(?:-[a-z0-9]+)*$"

function validate_slug(value, label::AbstractString, source::AbstractString)
    value isa AbstractString && !isempty(value) || throw(ArgumentError(
        "$label in $(basename(source)) must be a non-empty string"
    ))
    occursin(SLUG_PATTERN, value) || throw(ArgumentError(
        "$label $(repr(value)) must contain lowercase words separated by hyphens"
    ))
    return value
end

function validate_deck_page(page, deck_id::AbstractString, source::AbstractString)
    page isa NamedTuple || throw(ArgumentError(
        "deck $(repr(deck_id)) in $(basename(source)) contains a page that is not a named tuple"
    ))
    missing = filter(field -> !hasproperty(page, field), REQUIRED_DECK_PAGE_FIELDS)
    isempty(missing) || throw(ArgumentError(
        "page in deck $(repr(deck_id)) is missing fields: $(join(missing, ", "))"
    ))
    validate_slug(page.id, "page id", source)
    page.title isa AbstractString && !isempty(page.title) || throw(ArgumentError(
        "page title in deck $(repr(deck_id)) must be a non-empty string"
    ))
    page.render isa Bool || throw(ArgumentError(
        "page render flag in deck $(repr(deck_id)) must be true or false"
    ))
    page.class isa AbstractString || throw(ArgumentError(
        "page class in deck $(repr(deck_id)) must be a string"
    ))
    page.build isa Function || throw(ArgumentError(
        "page build entry in deck $(repr(deck_id)) must be a function"
    ))
    return page
end

function validate_deck_descriptor(deck, source::AbstractString)
    deck isa NamedTuple || throw(ArgumentError(
        "deck file $(basename(source)) must return a named tuple"
    ))
    missing = filter(field -> !hasproperty(deck, field), REQUIRED_DECK_FIELDS)
    isempty(missing) || throw(ArgumentError(
        "deck file $(basename(source)) is missing fields: $(join(missing, ", "))"
    ))
    validate_slug(deck.id, "deck id", source)
    deck.group isa AbstractString && !isempty(deck.group) || throw(ArgumentError(
        "deck group in $(basename(source)) must be a non-empty string"
    ))
    deck.title isa AbstractString && !isempty(deck.title) || throw(ArgumentError(
        "deck title in $(basename(source)) must be a non-empty string"
    ))
    deck.order isa Real && isfinite(deck.order) || throw(ArgumentError(
        "deck order in $(basename(source)) must be finite"
    ))
    deck.render isa Bool || throw(ArgumentError(
        "deck render flag in $(basename(source)) must be true or false"
    ))
    deck.class isa AbstractString || throw(ArgumentError(
        "deck class in $(basename(source)) must be a string"
    ))
    deck.setup isa Function || throw(ArgumentError(
        "deck setup entry in $(basename(source)) must be a function"
    ))
    deck.pages isa Tuple || deck.pages isa AbstractVector ||
        throw(ArgumentError(
            "deck pages in $(basename(source)) must be a tuple or vector"
        ))
    isempty(deck.pages) && throw(ArgumentError(
        "deck $(repr(deck.id)) must declare at least one page"
    ))
    pages = map(page -> validate_deck_page(page, deck.id, source), deck.pages)
    page_ids = getproperty.(pages, :id)
    duplicates = unique(filter(id -> count(==(id), page_ids) > 1, page_ids))
    isempty(duplicates) || throw(ArgumentError(
        "duplicate page ids in deck $(repr(deck.id)): $(join(duplicates, ", "))"
    ))
    return merge(deck, (; pages = Tuple(pages), source = abspath(source)))
end

function load_deck_descriptors(directory::AbstractString = DECK_DIRECTORY)
    isdir(directory) || throw(ArgumentError("deck directory does not exist: $directory"))
    sources = sort(filter(path -> endswith(path, ".jl"), readdir(directory; join = true)))
    isempty(sources) && throw(ArgumentError("no Julia deck files found in $directory"))
    decks = map(sources) do source
        validate_deck_descriptor(Base.include(@__MODULE__, source), source)
    end
    ids = getproperty.(decks, :id)
    duplicates = unique(filter(id -> count(==(id), ids) > 1, ids))
    isempty(duplicates) || throw(ArgumentError(
        "duplicate deck ids: $(join(duplicates, ", "))"
    ))
    return sort(decks; by = deck -> (deck.order, basename(deck.source)))
end

rendered_deck_pages(deck) = filter(page -> page.render, deck.pages)

function select_rendered_decks(decks)
    rendered = filter(deck -> deck.render, decks)
    isempty(rendered) && throw(ArgumentError(
        "no showcase decks are enabled; set `render = true` in at least one pages/*.jl file"
    ))
    empty_decks = filter(deck -> isempty(rendered_deck_pages(deck)), rendered)
    isempty(empty_decks) || throw(ArgumentError(
        "rendered decks must contain at least one rendered page: $(join(getproperty.(empty_decks, :id), ", "))"
    ))
    return rendered
end

function find_deck(id::AbstractString, decks)
    return only(filter(deck -> deck.id == id, decks))
end

function find_deck_page(id::AbstractString, deck)
    return only(filter(page -> page.id == id, deck.pages))
end

const DECK_DESCRIPTORS = load_deck_descriptors()
const RENDERED_DECKS = select_rendered_decks(DECK_DESCRIPTORS)

function app_assets()
    return (
        lato_css = Asset(LATO_STYLESHEET),
        julia_mono_css = Asset(JULIA_MONO_STYLESHEET),
        app_css = read(joinpath(@__DIR__, "assets", "theme.css"), String),
        logo = Asset(DOCUMENTATION_LOGO)
    )
end

deck_route(deck, index::Integer) = index == 1 ? "/" : "/$(deck.id)"
deck_href(deck, index::Integer) = index == 1 ? "./" : "./$(deck.id)"
function deck_page_href(deck, deck_index::Integer, page)
    "$(deck_href(deck, deck_index))#$(page.id)"
end

function deck_index(deck, decks)
    index = findfirst(candidate -> candidate.id == deck.id, decks)
    isnothing(index) && throw(ArgumentError("deck $(repr(deck.id)) is not rendered"))
    return index
end

function deck_page_refs(decks)
    refs = NamedTuple[]
    for (deck_position, deck) in enumerate(decks)
        for (page_position, page) in enumerate(rendered_deck_pages(deck))
            push!(refs, (; deck, page, deck_position, page_position))
        end
    end
    return refs
end

function deck_page_ref_index(deck, page, refs)
    index = findfirst(ref -> ref.deck.id == deck.id && ref.page.id == page.id, refs)
    isnothing(index) && throw(ArgumentError(
        "page $(repr(page.id)) in deck $(repr(deck.id)) is not rendered"
    ))
    return index
end

function navigation(assets, decks, current_deck, current_page)
    groups = unique(deck.group for deck in decks)
    group_nodes = map(groups) do group
        deck_nodes = map(
            filter(pair -> pair[2].group == group, collect(enumerate(decks)))
        ) do (deck_position, deck)
            pages = rendered_deck_pages(deck)
            page_entries = map(pages) do page
                active = deck.id == current_deck.id && page.id == current_page.id
                DOM.li(
                    DOM.a(
                    page.title;
                    href = deck_page_href(deck, deck_position, page),
                    class = active ? "lc-nav-entry is-active" : "lc-nav-entry",
                    dataDeckId = deck.id,
                    dataDeckTitle = deck.title,
                    dataPageId = page.id,
                    dataPageTitle = page.title,
                    ariaCurrent = active ? "page" : "false"
                )
                )
            end
            length(pages) == 1 && return first(page_entries)
            return DOM.li(
                DOM.a(
                    deck.title;
                    href = deck_page_href(deck, deck_position, first(pages)),
                    class = "lc-nav-deck-entry",
                    dataDeckId = deck.id
                ),
                DOM.ul(page_entries...; class = "lc-nav-list lc-nav-page-list");
                class = "lc-nav-deck"
            )
        end
        DOM.div(
            DOM.p(group; class = "lc-nav-group-title"),
            DOM.ul(deck_nodes...; class = "lc-nav-list");
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
                placeholder = "Search decks and pages",
                autocomplete = "off",
                ariaLabel = "Search decks and pages"
            );
            class = "lc-search"
        ),
        DOM.nav(
            group_nodes...;
            id = "page-navigation",
            class = "lc-sidebar-nav",
            ariaLabel = "Manual decks and pages"
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

function page_control(label; id, target, title, aria_label)
    isnothing(target) && return DOM.a(
        label;
        id,
        class = "lc-navbar-button is-disabled",
        title,
        ariaLabel = aria_label,
        ariaDisabled = "true"
    )
    return DOM.a(
        label;
        id,
        href = deck_page_href(target.deck, target.deck_position, target.page),
        class = "lc-navbar-button",
        title,
        ariaLabel = aria_label
    )
end

function navbar(decks, current_deck, current_page; preparing::Bool = false)
    refs = deck_page_refs(decks)
    current_index = deck_page_ref_index(current_deck, current_page, refs)
    previous = current_index == 1 ? nothing : refs[current_index - 1]
    next = current_index == length(refs) ? nothing : refs[current_index + 1]
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
                current_deck.group;
                id = "lc-current-group",
                class = "lc-breadcrumb-parent"
            ),
            DOM.span("›"; class = "lc-breadcrumb-separator", ariaHidden = "true"),
            DOM.span(
                current_deck.title;
                id = "lc-current-deck",
                class = "lc-breadcrumb-parent"
            ),
            DOM.span("›"; class = "lc-breadcrumb-separator", ariaHidden = "true"),
            DOM.span(current_page.title; id = "lc-current-title");
            class = "lc-breadcrumb",
            ariaLabel = "Current page"
        ),
        DOM.div(
            page_control(
                "‹";
                id = "lc-previous-page",
                target = previous,
                title = "Previous page",
                aria_label = "Previous page"
            ),
            DOM.output(
                "$(current_index) / $(length(refs))";
                id = "lc-page-position",
                class = "lc-page-position",
                ariaLive = "polite"
            ),
            page_control(
                "›";
                id = "lc-next-page",
                target = next,
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
        decks,
        current_deck,
        current_page,
        manual_root;
        preparing::Bool = false
)
    main = DOM.main(
        navbar(decks, current_deck, current_page; preparing),
        DOM.div(manual_root; class = "lc-page-frame");
        class = "lc-main"
    )
    return DOM.div(
        assets.lato_css,
        assets.julia_mono_css,
        DOM.style(assets.app_css),
        navigation(assets, decks, current_deck, current_page),
        main;
        id = "lc-live-docs",
        class = preparing ? "lc-documenter is-preparing" : "lc-documenter",
        dataCurrentDeckId = current_deck.id,
        ariaBusy = string(preparing)
    )
end

function render_deck_page(session::Session, deck, page, state, index::Integer)
    result = page.build(session, state)
    result isa NamedTuple && hasproperty(result, :body) || throw(ArgumentError(
        "page $(repr(page.id)) in deck $(repr(deck.id)) must return a named tuple with a `body` field"
    ))
    body = result.body isa Tuple ? result.body : (result.body,)
    attributes = (;
        id = "deck-$(deck.id)-page-$(page.id)",
        class = index == 1 ? "lc-page is-active $(page.class) $(deck.class)" :
                "lc-page $(page.class) $(deck.class)",
        dataDeckId = deck.id,
        dataDeckTitle = deck.title,
        dataPageId = page.id,
        dataPageTitle = page.title,
        dataPageIndex = string(index),
        ariaHidden = string(index != 1)
    )
    index != 1 && (attributes = merge(attributes, (; hidden = true)))
    return DOM.section(body...; attributes...)
end

function render_deck_pages(session::Session, deck)
    state = deck.setup(session)
    pages = rendered_deck_pages(deck)
    rendered = map(enumerate(pages)) do (index, page)
        render_deck_page(session, deck, page, state, index)
    end
    return (; nodes = rendered, state)
end

function preparation_page(session::Session, deck, page)
    body = slide(
        deck.title,
        article(md"""
        This deck is waiting for its process-wide numerical model. The server remains responsive and this status is updated from the actual Julia preparation phases.
        """),
        preparation_status(session);
        lede = md"The application case is being prepared once and will be reused by every deck session."
    )
    return DOM.section(
        body...;
        id = "deck-$(deck.id)-page-$(page.id)",
        class = "lc-page is-active $(page.class) $(deck.class)",
        dataDeckId = deck.id,
        dataDeckTitle = deck.title,
        dataPageId = page.id,
        dataPageTitle = page.title,
        dataPageIndex = "1",
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
    const currentTitle = root.querySelector("#lc-current-title");
    const pagePosition = root.querySelector("#lc-page-position");
    const navEntries = Array.from(root.querySelectorAll(".lc-nav-entry"));
    const navGroups = Array.from(root.querySelectorAll(".lc-nav-group"));
    const deckSections = Array.from(manualRoot.querySelectorAll(".lc-page"));
    const currentDeckId = root.dataset.currentDeckId;
    const currentDeckEntries = navEntries.filter(
        (entry) => entry.dataset.deckId === currentDeckId
    );
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
        window.removeEventListener("hashchange", activateDeckPage);
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
        if (
            link instanceof HTMLAnchorElement &&
            link.href &&
            link.getAttribute("aria-disabled") !== "true"
        ) {
            window.location.assign(link.href);
        }
    };

    const configurePageControl = (link, entry) => {
        if (!(link instanceof HTMLAnchorElement)) {
            return;
        }
        if (!entry) {
            link.removeAttribute("href");
            link.classList.add("is-disabled");
            link.setAttribute("aria-disabled", "true");
            return;
        }
        link.href = entry.href;
        link.classList.remove("is-disabled");
        link.setAttribute("aria-disabled", "false");
    };

    const activateDeckPage = () => {
        const requestedId = decodeURIComponent(window.location.hash.slice(1));
        const section = deckSections.find((candidate) =>
            candidate.dataset.pageId === requestedId
        ) || deckSections[0];
        if (!section) {
            return;
        }
        deckSections.forEach((candidate) => {
            const active = candidate === section;
            candidate.hidden = !active;
            candidate.classList.toggle("is-active", active);
            candidate.setAttribute("aria-hidden", String(!active));
        });
        navEntries.forEach((entry) => {
            const active = entry.dataset.deckId === currentDeckId &&
                entry.dataset.pageId === section.dataset.pageId;
            entry.classList.toggle("is-active", active);
            entry.setAttribute("aria-current", active ? "page" : "false");
        });
        currentTitle.textContent = section.dataset.pageTitle;
        root.dataset.currentPageId = section.dataset.pageId;

        const currentEntry = currentDeckEntries.find(
            (entry) => entry.dataset.pageId === section.dataset.pageId
        );
        const globalIndex = navEntries.indexOf(currentEntry);
        configurePageControl(previousLink, navEntries[globalIndex - 1]);
        configurePageControl(nextLink, navEntries[globalIndex + 1]);
        pagePosition.textContent = `${globalIndex + 1} / ${navEntries.length}`;

        if (mobileSidebar.matches) {
            root.classList.remove("is-sidebar-open");
            sidebarToggle.setAttribute("aria-expanded", "false");
        }
        settleViewport();
    };

    window.addEventListener("hashchange", activateDeckPage);
    activateDeckPage();

    search.addEventListener("input", () => {
        const query = search.value.trim().toLocaleLowerCase();
        navEntries.forEach((entry) => {
            const group = entry.closest(".lc-nav-group").dataset.navGroup;
            const searchable = `${group} ${entry.dataset.deckTitle} ${entry.textContent}`
                .toLocaleLowerCase();
            entry.closest("li").hidden = query.length > 0 && !searchable.includes(query);
        });
        root.querySelectorAll(".lc-nav-deck").forEach((deck) => {
            const visiblePage = Array.from(deck.querySelectorAll(".lc-nav-entry")).some(
                (entry) => !entry.closest("li").hidden
            );
            deck.hidden = !visiblePage;
        });
        navGroups.forEach((group) => {
            const visibleEntry = Array.from(group.querySelectorAll(".lc-nav-entry")).some(
                (entry) => !entry.closest("li").hidden
            );
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

function manual(session::Session, assets, deck; decks = RENDERED_DECKS)
    decks = select_rendered_decks(decks)
    deck_index(deck, decks)
    pages = rendered_deck_pages(deck)
    current_page = first(pages)
    preparing = hasproperty(deck, :ready) && !deck.ready()
    rendered = preparing ?
               (;
        nodes = [preparation_page(session, deck, current_page)], state = nothing) :
               render_deck_pages(session, deck)
    page_nodes = DOM.div(
        rendered.nodes...;
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
    root = documenter_shell(
        assets,
        decks,
        deck,
        current_page,
        manual_root;
        preparing
    )
    preparing && bind_preparation_reload(session)
    return bind_shell_behavior(session, root)
end

function cable_app(
        assets = app_assets();
        decks = RENDERED_DECKS,
        deck = first(select_rendered_decks(decks))
)
    decks = select_rendered_decks(decks)
    return App(
        session -> manual(session, assets, deck; decks);
        title = "LineCableModels live manual",
        indicator = ConnectionIndicator(; size = 8, top = "20px", right = "16px")
    )
end

function cable_routes(
        assets = app_assets();
        decks = RENDERED_DECKS
)
    decks = select_rendered_decks(decks)
    pairs = map(enumerate(decks)) do (index, deck)
        deck_route(deck, index) => cable_app(assets; decks, deck)
    end
    return Routes(pairs...)
end

function start_deck_preparations!(decks = RENDERED_DECKS)
    decks = select_rendered_decks(decks)
    return map(
        deck -> hasproperty(deck, :prepare) ? deck.prepare() : nothing,
        decks
    )
end
