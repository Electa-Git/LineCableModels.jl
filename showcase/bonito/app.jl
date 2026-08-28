using Bonito
using LineCableModels
using Makie

const REVEAL_BASE_URL = "https://cdn.jsdelivr.net/npm/reveal.js@6.0.1"
const LATO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/lato-font/3.0.0/css/lato-font.min.css"
const JULIA_MONO_STYLESHEET = "https://cdnjs.cloudflare.com/ajax/libs/juliamono/0.050/juliamono.min.css"
const DOCUMENTATION_LOGO = joinpath(
    @__DIR__, "..", "..", "docs", "src", "assets", "logo.svg")
const CORE_MATERIAL, INSULATION_MATERIAL = let materials = MaterialsLibrary()
    Material(materials, :copper), Material(materials, :xlpe)
end

# Frozen visual samples from the current material-preview theme. The showcase
# intentionally does not call or extend the package's PlotBuilder machinery.
const PLOT_BACKGROUND = RGBf(0.9, 0.9, 0.9)
const COPPER_COLOR = RGBf(0.9093547, 0.8964516, 0.8706453)
const XLPE_COLOR = RGBf(0.06460433, 0.1270558, 0.1206059)
const SLIDE_DESCRIPTORS = (
    (id = "overview", group = "Showcase", title = "Live technical manual"),
    (
        id = "core-and-insulation",
        group = "Cable design",
        title = "Core and insulation"
    )
)

function design(core_radius_mm::Real, insulation_thickness_mm::Real)
    isfinite(core_radius_mm) && core_radius_mm > zero(core_radius_mm) ||
        throw(DomainError(core_radius_mm, "core radius must be positive and finite"))
    isfinite(insulation_thickness_mm) &&
    insulation_thickness_mm > zero(insulation_thickness_mm) || throw(DomainError(
        insulation_thickness_mm,
        "insulation thickness must be positive and finite"
    ))

    core_radius_m = float(core_radius_mm) * 1.0e-3
    insulation_thickness_m = float(insulation_thickness_mm) * 1.0e-3
    core = Conductor.Solid(:core_conductor, CORE_MATERIAL; r = core_radius_m)
    insulation = Insulator.Shell(
        :insulation,
        INSULATION_MATERIAL;
        t = insulation_thickness_m
    )
    root = Stack(Group(:core, core), insulation)
    return build(CableDesign, "interactive-demo", root)
end

function geometry(cable::CableDesign)
    regions = cable.geometry.regions
    core = only(filter(region -> region.source.tag === :core_conductor, regions))
    insulation = only(filter(region -> region.source.tag === :insulation, regions))
    return (
        core_r_in_m = r_in(core.primitive),
        core_r_ex_m = r_ex(core.primitive),
        core_area_m2 = area(core.primitive),
        insulation_r_in_m = r_in(insulation.primitive),
        insulation_r_ex_m = r_ex(insulation.primitive)
    )
end

function cable_figure(design_state::Observable; session = nothing)
    observe(f) = isnothing(session) ? map(f, design_state) : map(f, session, design_state)
    insulation_circle = observe() do cable
        radius_mm = geometry(cable).insulation_r_ex_m * 1.0e3
        Circle(Point2f(0), Float32(radius_mm))
    end
    core_circle = observe() do cable
        radius_mm = geometry(cable).core_r_ex_m * 1.0e3
        Circle(Point2f(0), Float32(radius_mm))
    end

    figure = Figure(size = (760, 520), backgroundcolor = PLOT_BACKGROUND)
    axis = Axis(
        figure[1, 1];
        aspect = DataAspect(),
        backgroundcolor = PLOT_BACKGROUND
    )
    insulation_plot = poly!(
        axis,
        insulation_circle;
        color = XLPE_COLOR,
        strokecolor = :black,
        strokewidth = 0.5,
        label = "XLPE insulation"
    )
    core_plot = poly!(
        axis,
        core_circle;
        color = COPPER_COLOR,
        strokecolor = :black,
        strokewidth = 0.5,
        label = "Copper core"
    )
    hidedecorations!(axis)
    hidespines!(axis)
    limits!(axis, -50, 50, -50, 50)
    axislegend(
        axis;
        position = :rt,
        backgroundcolor = (:white, 0.9),
        framecolor = RGBf(0.55, 0.55, 0.55),
        framevisible = true,
        labelcolor = RGBf(0.15, 0.15, 0.15),
        padding = (8, 8, 6, 6)
    )
    return (;
        figure,
        axis,
        insulation_plot,
        core_plot,
        insulation_circle,
        core_circle
    )
end

function reveal_assets()
    return (
        reveal_js = ES6Module("$REVEAL_BASE_URL/dist/reveal.mjs"),
        notes_js = ES6Module("$REVEAL_BASE_URL/dist/plugin/notes.mjs"),
        reset_css = Asset("$REVEAL_BASE_URL/dist/reset.css"),
        reveal_css = Asset("$REVEAL_BASE_URL/dist/reveal.css"),
        lato_css = Asset(LATO_STYLESHEET),
        julia_mono_css = Asset(JULIA_MONO_STYLESHEET),
        app_css = Asset(joinpath(@__DIR__, "assets", "theme.css")),
        logo = Asset(DOCUMENTATION_LOGO)
    )
end

function navigation(assets)
    groups = unique(descriptor.group for descriptor in SLIDE_DESCRIPTORS)
    group_nodes = map(groups) do group
        entries = map(
            filter(pair -> pair[2].group == group,
            collect(enumerate(SLIDE_DESCRIPTORS)))
        ) do (index, descriptor)
            button_class = index == 1 ? "lc-nav-entry is-active" : "lc-nav-entry"
            DOM.li(
                DOM.button(
                descriptor.title;
                type = "button",
                class = button_class,
                dataSlideIndex = string(index - 1),
                dataSlideId = descriptor.id,
                ariaCurrent = index == 1 ? "page" : "false"
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
                id = "slide-search",
                type = "search",
                placeholder = "Search slides",
                autocomplete = "off",
                ariaLabel = "Search slides"
            );
            class = "lc-search"
        ),
        DOM.nav(
            group_nodes...;
            id = "slide-navigation",
            class = "lc-sidebar-nav",
            ariaLabel = "Presentation sections"
        ),
        DOM.div(
            DOM.span("Bonito showcase"),
            DOM.span("reveal.js 6.0.1");
            class = "lc-sidebar-footer"
        );
        id = "lc-sidebar",
        class = "lc-sidebar",
        ariaLabel = "Slide navigation"
    )
end

function slide(descriptor, content...; class)
    return DOM.section(
        content...;
        id = "slide-$(descriptor.id)",
        class,
        dataNavGroup = descriptor.group,
        dataNavTitle = descriptor.title
    )
end

function deck(session::Session, assets)
    core_slider = Bonito.Slider(
        collect(5.0:0.1:25.0);
        value = 12.5,
        id = "core-radius-input",
        ariaLabel = "Core radius in millimetres"
    )
    insulation_slider = Bonito.Slider(
        collect(1.0:0.1:20.0);
        value = 8.0,
        id = "insulation-thickness-input",
        ariaLabel = "Insulation thickness in millimetres"
    )
    design_state = map(design, session, core_slider.value, insulation_slider.value)
    geometry_state = map(geometry, session, design_state)
    plot = cable_figure(design_state; session)

    core_value = map(
        values -> "$(round(values.core_r_ex_m * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    insulation_value = map(
        values -> "$(round((values.insulation_r_ex_m - values.insulation_r_in_m) * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    outer_diameter = map(
        values -> "$(round(2values.insulation_r_ex_m * 1.0e3; digits = 1)) mm",
        session,
        geometry_state
    )
    core_area = map(
        values -> "$(round(values.core_area_m2 * 1.0e6; digits = 2)) mm²",
        session,
        geometry_state
    )
    design_type = map(
        cable -> "$(nameof(typeof(cable))){$(eltype(cable)), …}",
        session,
        design_state
    )

    overview_descriptor, interactive_descriptor = SLIDE_DESCRIPTORS
    overview_slide = slide(
        overview_descriptor,
        DOM.article(
            DOM.h1("Live technical manual"),
            DOM.p(
                "This showcase behaves like a documentation page while reveal.js quietly provides presentation navigation.";
                class = "lc-lede"
            ),
            DOM.div(
                DOM.strong("Interactive example"),
                DOM.p(
                    "The cable on the next page is rebuilt from ordinary LineCableModels constructors whenever either native Bonito slider changes."
                );
                class = "lc-doc-note"
            ),
            DOM.h2("Showcase boundary"),
            DOM.ul(
                DOM.li("The application shell mirrors the current Documenter dark theme."),
                DOM.li("The WGLMakie figure is a small showcase-local renderer."),
                DOM.li("No PlotBuilder integration or package API adapter is introduced.")
            ),
            DOM.pre(
                DOM.code(
                    "CableDesign\n└── Stack\n    ├── Group(:core, Conductor.Solid(...))\n    └── Insulator.Shell(...)"
                );
                class = "lc-code-block"
            ),
            DOM.p(
                "Use the sidebar, arrow keys, or swipe to move between pages.";
                class = "lc-page-hint"
            );
            class = "lc-article"
        ),
        DOM.aside(
            "Frame this as a live technical manual rather than a conventional slide deck.";
            class = "notes"
        );
        class = "lc-slide lc-overview-slide"
    )

    control_panel = DOM.div(
        DOM.label(
            DOM.div(
                DOM.span("Core radius"),
                DOM.output(core_value),
                class = "lc-control-heading"
            ),
            core_slider;
            id = "core-radius-control",
            class = "lc-control"
        ),
        DOM.label(
            DOM.div(
                DOM.span("Insulation thickness"),
                DOM.output(insulation_value),
                class = "lc-control-heading"
            ),
            insulation_slider;
            id = "insulation-thickness-control",
            class = "lc-control"
        ),
        DOM.dl(
            DOM.div(DOM.dt("Conductor radius"), DOM.dd(core_value)),
            DOM.div(DOM.dt("Insulation thickness"), DOM.dd(insulation_value)),
            DOM.div(DOM.dt("Outer cable diameter"), DOM.dd(outer_diameter)),
            DOM.div(DOM.dt("Conductor area"), DOM.dd(core_area));
            class = "lc-values"
        ),
        DOM.p(DOM.span("Object: "), DOM.code(design_type); class = "lc-diagnostic");
        class = "lc-controls"
    )

    interactive_slide = slide(
        interactive_descriptor,
        DOM.div(
            DOM.h1("Core and insulation"),
            DOM.p(
                "Adjust the dimensions to reconstruct the domain object and update the cross-section in place."
            );
            class = "lc-page-heading"
        ),
        DOM.div(
            control_panel,
            DOM.div(
                WGLMakie.WithConfig(plot.figure; resize_to = :parent);
                id = "cable-plot",
                class = "lc-plot-host"
            );
            class = "lc-interactive-grid"
        ),
        DOM.aside(
            "Every slider input reconstructs the CableDesign; only the existing Makie circle observables are updated.";
            class = "notes"
        );
        class = "lc-slide lc-interactive-slide"
    )

    slides = DOM.div(
        overview_slide,
        interactive_slide;
        id = "linecable-slides",
        class = "slides"
    )
    reveal_root = DOM.div(
        slides;
        id = "linecable-app",
        class = "reveal lc-deck",
        tabindex = "0",
        ariaLabel = "LineCableModels live manual pages"
    )
    main = DOM.main(
        DOM.header(
            DOM.button(
                "☰";
                id = "lc-sidebar-toggle",
                class = "lc-sidebar-toggle",
                type = "button",
                ariaLabel = "Toggle slide navigation",
                ariaControls = "lc-sidebar",
                ariaExpanded = "false"
            ),
            DOM.nav(
                DOM.span(
                    first(SLIDE_DESCRIPTORS).group;
                    id = "lc-current-group",
                    class = "lc-breadcrumb-parent"
                ),
                DOM.span("›"; class = "lc-breadcrumb-separator", ariaHidden = "true"),
                DOM.span(first(SLIDE_DESCRIPTORS).title; id = "lc-current-title");
                class = "lc-breadcrumb",
                ariaLabel = "Current page"
            ),
            DOM.span("Interactive showcase"; class = "lc-navbar-status");
            class = "lc-navbar"
        ),
        DOM.div(reveal_root; class = "lc-deck-frame");
        class = "lc-main"
    )
    root = DOM.div(
        assets.lato_css,
        assets.julia_mono_css,
        assets.reset_css,
        assets.reveal_css,
        assets.app_css,
        navigation(assets),
        main;
        id = "lc-live-docs",
        class = "lc-documenter"
    )

    Bonito.onload(
        session, root,
        js"""
(root) => {
    if (root.__lineCableReveal) {
        return;
    }

    const revealRoot = root.querySelector("#linecable-app");
    const main = root.querySelector(".lc-main");
    const controls = root.querySelector(".lc-controls");
    const search = root.querySelector("#slide-search");
    const sidebarToggle = root.querySelector("#lc-sidebar-toggle");
    const currentGroup = root.querySelector("#lc-current-group");
    const currentTitle = root.querySelector("#lc-current-title");
    const navEntries = Array.from(root.querySelectorAll(".lc-nav-entry"));
    const navGroups = Array.from(root.querySelectorAll(".lc-nav-group"));

    root.dataset.revealState = "loading";
    Promise.all([$(assets.reveal_js), $(assets.notes_js)]).then(([revealModule, notesModule]) => {
        if (root.__lineCableReveal) {
            return;
        }

        const Reveal = revealModule.default;
        const RevealNotes = notesModule.default;
        const deck = new Reveal(revealRoot, {
            width: 1280,
            height: 720,
            embedded: true,
            disableLayout: true,
            center: false,
            controls: false,
            progress: false,
            slideNumber: false,
            transition: "none",
            backgroundTransition: "none",
            hash: true,
            keyboardCondition: "focused",
            touch: true,
            overview: true,
            plugins: [RevealNotes],
        });
        root.__lineCableReveal = deck;

        let initialized = false;
        let dispatchingResize = false;
        let layoutFrame = 0;
        const relayout = () => {
            if (!initialized) {
                return;
            }
            window.cancelAnimationFrame(layoutFrame);
            layoutFrame = window.requestAnimationFrame(() => {
                deck.layout();
                window.requestAnimationFrame(() => {
                    dispatchingResize = true;
                    window.dispatchEvent(new Event("resize"));
                    dispatchingResize = false;
                });
            });
        };
        const onWindowResize = () => {
            if (!dispatchingResize) {
                relayout();
            }
        };

        const updateNavigation = () => {
            const indices = deck.getIndices();
            navEntries.forEach((entry) => {
                const active = Number(entry.dataset.slideIndex) === indices.h;
                entry.classList.toggle("is-active", active);
                entry.setAttribute("aria-current", active ? "page" : "false");
            });

            const activeSlide = deck.getCurrentSlide();
            if (activeSlide) {
                currentGroup.textContent = activeSlide.dataset.navGroup || "Manual";
                currentTitle.textContent = activeSlide.dataset.navTitle || "";
            }
        };

        navEntries.forEach((entry) => {
            entry.addEventListener("click", () => {
                deck.slide(Number(entry.dataset.slideIndex), 0);
                root.classList.remove("is-sidebar-open");
                sidebarToggle.setAttribute("aria-expanded", "false");
                revealRoot.focus({ preventScroll: true });
            });
        });

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

        root.addEventListener("keydown", (event) => {
            if ((event.ctrlKey || event.metaKey) && event.key === "/") {
                event.preventDefault();
                search.focus();
            }
        });
        sidebarToggle.addEventListener("click", () => {
            const open = root.classList.toggle("is-sidebar-open");
            sidebarToggle.setAttribute("aria-expanded", String(open));
        });
        revealRoot.addEventListener("pointerdown", (event) => {
            if (!event.target.closest("input, button, a, select, textarea")) {
                revealRoot.focus({ preventScroll: true });
            }
        });

        ["pointerdown", "pointermove", "touchstart", "touchmove"].forEach((name) => {
            controls.addEventListener(name, (event) => event.stopPropagation(), { passive: true });
        });
        controls.addEventListener("keydown", (event) => {
            if (event.target.matches('input[type="range"]') && event.key.startsWith("Arrow")) {
                event.stopPropagation();
            }
        });

        deck.on("ready", () => {
            initialized = true;
            root.dataset.revealState = "ready";
            updateNavigation();
            relayout();
        });
        deck.on("slidechanged", () => {
            updateNavigation();
            relayout();
        });
        deck.on("overviewhidden", relayout);
        window.addEventListener("resize", onWindowResize);
        const resizeObserver = new ResizeObserver(relayout);
        resizeObserver.observe(main);
        root.__lineCableResizeObserver = resizeObserver;
        deck.initialize().catch((error) => {
            root.dataset.revealState = "error";
            console.error("Unable to initialize reveal.js", error);
        });
    }).catch((error) => {
        root.dataset.revealState = "error";
        console.error("Unable to load reveal.js", error);
    });
}
""")
    return root
end

function cable_app(assets = reveal_assets())
    return App(
        session -> deck(session, assets);
        title = "LineCableModels live manual",
        indicator = ConnectionIndicator()
    )
end
