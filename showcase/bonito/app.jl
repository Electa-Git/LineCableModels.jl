using Bonito
using LineCableModels
using Makie

const REVEAL_BASE_URL = "https://cdn.jsdelivr.net/npm/reveal.js@6.0.1"
const CORE_MATERIAL, INSULATION_MATERIAL = let materials = MaterialsLibrary()
    Material(materials, :copper), Material(materials, :xlpe)
end

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
        core_r_in_m = r_in(core.shape),
        core_r_ex_m = r_ex(core.shape),
        core_area_m2 = area(core.shape),
        insulation_r_in_m = r_in(insulation.shape),
        insulation_r_ex_m = r_ex(insulation.shape)
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

    background = RGBf(0.035, 0.043, 0.055)
    figure = Figure(size = (760, 520), backgroundcolor = background)
    axis = Axis(
        figure[1, 1];
        aspect = DataAspect(),
        backgroundcolor = background
    )
    insulation_plot = poly!(
        axis,
        insulation_circle;
        color = RGBf(0.22, 0.55, 0.85),
        strokecolor = RGBf(0.62, 0.82, 1.0),
        strokewidth = 2,
        label = "XLPE insulation"
    )
    core_plot = poly!(
        axis,
        core_circle;
        color = RGBf(0.82, 0.43, 0.18),
        strokecolor = RGBf(1.0, 0.72, 0.43),
        strokewidth = 2,
        label = "Copper core"
    )
    hidedecorations!(axis)
    hidespines!(axis)
    limits!(axis, -50, 50, -50, 50)
    axislegend(
        axis;
        position = :rt,
        framevisible = false,
        labelcolor = :white
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
        theme_css = Asset("$REVEAL_BASE_URL/dist/theme/black.css"),
        app_css = Asset(joinpath(@__DIR__, "assets", "theme.css"))
    )
end

function deck(session::Session, assets)
    core_slider = Bonito.Slider(
        collect(5.0:0.1:25.0);
        value = 12.5,
        id = "core-radius-input",
        aria_label = "Core radius in millimetres"
    )
    insulation_slider = Bonito.Slider(
        collect(1.0:0.1:20.0);
        value = 8.0,
        id = "insulation-thickness-input",
        aria_label = "Insulation thickness in millimetres"
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

    title_slide = DOM.section(
        DOM.div(
            DOM.p("Interactive geometry prototype"; class = "lc-eyebrow"),
            DOM.h1("LineCableModels interactive cable geometry"),
            DOM.h2("Bonito + reveal.js + WGLMakie"),
            DOM.p(
                "The next slide rebuilds and renders a live CableDesign.";
                class = "lc-lede"
            );
            class = "lc-title-content"
        ),
        DOM.aside(
            "Introduce the stack and advance to the live domain-object demonstration.";
            class = "notes"
        );
        class = "lc-slide lc-title-slide"
    )

    control_panel = DOM.div(
        DOM.div(
            DOM.div(
                DOM.span("Core radius"),
                DOM.output(core_value),
                class = "lc-control-heading"
            ),
            core_slider;
            id = "core-radius-control",
            class = "lc-control"
        ),
        DOM.div(
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

    interactive_slide = DOM.section(
        DOM.h2("Interactive cable design"),
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

    root = DOM.div(
        assets.reset_css,
        assets.reveal_css,
        assets.theme_css,
        assets.app_css,
        DOM.div(title_slide, interactive_slide; id = "linecable-slides", class = "slides");
        id = "linecable-app",
        class = "reveal lc-deck"
    )

    Bonito.onload(
        session, root, js"""
 (root) => {
     if (root.__lineCableReveal) {
         return;
     }

     Promise.all([$(assets.reveal_js), $(assets.notes_js)]).then(([revealModule, notesModule]) => {
         if (root.__lineCableReveal) {
             return;
         }

         const Reveal = revealModule.default;
         const RevealNotes = notesModule.default;
         const deck = new Reveal(root, {
             width: 1280,
             height: 720,
             margin: 0.04,
             controls: true,
             progress: true,
             hash: true,
             slideNumber: "c/t",
             center: false,
             transition: "none",
             keyboard: true,
             touch: true,
             overview: true,
             plugins: [RevealNotes],
         });
         root.__lineCableReveal = deck;

         let dispatchingResize = false;
         const relayout = () => {
             deck.layout();
             window.requestAnimationFrame(() => {
                 dispatchingResize = true;
                 window.dispatchEvent(new Event("resize"));
                 dispatchingResize = false;
             });
         };
         const onWindowResize = () => {
             if (!dispatchingResize) {
                 relayout();
             }
         };

         const controls = root.querySelector(".lc-controls");
         ["pointerdown", "pointermove", "touchstart", "touchmove"].forEach((name) => {
             controls.addEventListener(name, (event) => event.stopPropagation(), { passive: true });
         });
         controls.addEventListener("keydown", (event) => {
             if (event.target.matches('input[type="range"]') && event.key.startsWith("Arrow")) {
                 event.stopPropagation();
             }
         });

         deck.on("ready", relayout);
         deck.on("slidechanged", relayout);
         window.addEventListener("resize", onWindowResize);
         deck.initialize();
     });
 }
 """)
    return root
end

function cable_app(assets = reveal_assets())
    return App(
        session -> deck(session, assets);
        title = "LineCableModels interactive cable geometry",
        indicator = ConnectionIndicator()
    )
end
