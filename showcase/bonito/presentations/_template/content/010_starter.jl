module StarterDeck

using Bonito
using Makie
using Markdown
using Observables: off
using WGLMakie
using ..PageAuthoring

const EXAMPLE_IMAGE = normpath(joinpath(@__DIR__, "..", "assets", "example-diagram.svg"))
const PREFLIGHT_STATE = preflight_state(
    label = "Starter resource not activated"
)

function activate_starter_resource!()
    phase = PREFLIGHT_STATE[].phase
    phase === :hot && return nothing
    phase ∉ (:cold, :error) && return nothing
    set_preflight!(
        PREFLIGHT_STATE,
        :preparing,
        0.5,
        "Preparing the starter resource…"
    )
    return @async begin
        yield()
        set_preflight!(PREFLIGHT_STATE, :hot, 1.0, "Starter resource hot")
    end
end

const PREFLIGHT_RESOURCE = preflight_resource(
    id = "starter-resource",
    title = "Starter process-wide resource",
    activate = activate_starter_resource!,
    state = PREFLIGHT_STATE
)

function example_figure(amplitude, session::Session)
    points = map(session, amplitude) do value
        [Point2f(x, value * sin(x)) for x in range(0, 2pi; length = 160)]
    end
    figure = Figure(size = (760, 460), backgroundcolor = RGBf(0.9, 0.9, 0.9))
    axis = Axis(
        figure[1, 1];
        xlabel = "x",
        ylabel = "a sin(x)",
        backgroundcolor = RGBf(0.9, 0.9, 0.9)
    )
    line = lines!(axis, points; color = RGBf(0.251, 0.388, 0.847), linewidth = 3)
    limits!(axis, 0, 2pi, -2.2, 2.2)
    return (; figure, axis, line, points)
end

function setup(session::Session)
    slider = Bonito.Slider(
        collect(0.1:0.1:2.0);
        value = 1.0,
        id = "starter-amplitude-input",
        ariaLabel = "Example amplitude"
    )
    choice = Bonito.Dropdown(
        ["Baseline", "Alternative", "Stress case"];
        id = "starter-choice-input",
        ariaLabel = "Example scenario"
    )
    run = Bonito.Button("Run action"; id = "starter-run-action")
    amplitude = map(value -> Float64(value), session, slider.value)
    amplitude_text = map(value -> string(round(value; digits = 1)), session, amplitude)
    choice_text = map(string, session, choice.value)
    action_count = Observable(0)
    action_text = map(
        count -> count == 0 ? "No action requested" : "Action requested $count time(s)",
        session,
        action_count
    )
    handler = on(run.value) do clicked
        clicked && (action_count[] += 1)
        return nothing
    end
    on(session.on_close) do _
        off(handler)
        return nothing
    end
    plot = example_figure(amplitude, session)
    return (;
        slider,
        choice,
        run,
        amplitude,
        amplitude_text,
        choice_text,
        action_count,
        action_text,
        plot
    )
end

panel(title, content) = webpart(prose(content); kind = :panel, title)

function prose_page(::Session, _)
    content = article(md"""
    ## Literate content

    Write ordinary Julia Markdown. It supports **emphasis**, lists, tables,
    code, links, and Bonito-rendered mathematics such as ``Z(j\omega)``.

    | Quantity | Meaning |
    |:--|:--|
    | ``Z`` | Series impedance |
    | ``Y`` | Shunt admittance |

    ```math
    \mathbf Z(\omega)=\mathbf R(\omega)+j\omega\mathbf L(\omega)
    ```

    ```julia
    value = map(transform, session, widget.value)
    ```
    """;
        hint = "Use layout_single when prose or one figure should own the canvas.")
    body = slide("Prose and mathematics", layout_single(content; height = "auto"))
    return (; body)
end

function two_columns_page(::Session, _)
    left = panel("Explanation", md"""
    Two-column pages work well for narrative beside evidence. Choose
    `:equal`, `:narrow_wide`, or `:wide_narrow` instead of writing CSS tracks.
    """)
    right = panel("Evidence", md"""
    - Reusable spacing
    - Responsive stacking
    - Stable print behavior
    - No presentation-local CSS
    """)
    body = slide(
        "Two columns",
        layout_two_columns(left, right; ratio = :narrow_wide, height = "24rem")
    )
    return (; body)
end

function top_split_page(::Session, _)
    top = panel("Summary", md"A full-width summary can introduce two detailed regions.")
    left = panel("Method", md"Describe the method, assumptions, or input data here.")
    right = panel("Result", md"Place a result, interpretation, or comparison here.")
    body = slide(
        "Top row and split detail",
        layout_top_split(top, left, right; ratio = :equal, height = "28rem")
    )
    return (; body)
end

function three_columns_page(::Session, _)
    columns = layout_three_columns(
        panel("Case A", md"Use concise comparable content."),
        panel("Case B", md"Keep corresponding concepts aligned."),
        panel("Case C", md"The layout stacks on narrow screens.");
        height = "24rem"
    )
    return (; body = slide("Three-way comparison", columns))
end

function controls_plot_page(::Session, state)
    controls = webpart(
        control(
            "Amplitude",
            state.slider;
            value = state.amplitude_text,
            id = "starter-amplitude-control"
        ),
        control(
            "Scenario",
            state.choice;
            value = state.choice_text,
            id = "starter-choice-control"
        ),
        state.run,
        status_line(state.action_text),
        value_list(
            "Amplitude" => state.amplitude_text,
            "Scenario" => state.choice_text
        );
        kind = :controls
    )
    plot = webpart(
        wgl_figure(state.plot.figure);
        kind = :plot,
        title = "Reactive Makie figure",
        meta = "persistent plot identity",
        id = "starter-live-plot"
    )
    footer = panel(
        "Pattern",
        md"Create widgets and figures once in `setup(session)`, then reuse them from every page builder that needs the same state."
    )
    body = slide(
        "Controls and live output",
        layout_controls_plot(
            controls,
            plot;
            footer,
            height = "32rem",
            class = "lc-template-controls-plot"
        );
        lede = md"Native Bonito widgets drive one persistent WGLMakie figure."
    )
    return (; body)
end

function split_bottom_page(::Session, state)
    left = panel("Selected scenario", md"The current choice is **$(state.choice_text)**.")
    right = panel("Live amplitude", md"The current amplitude is **$(state.amplitude_text)**.")
    bottom = panel(
        "Shared conclusion",
        md"A merged bottom region is useful for conclusions, provenance, or one common action."
    )
    body = slide(
        "Split content with footer",
        layout_split_bottom(left, right, bottom; height = "26rem")
    )
    return (; body)
end

function two_rows_page(::Session, _)
    top = panel("Network or overview", md"Use the shorter row for context or a system diagram.")
    bottom = panel("Detailed result", md"Use the taller row for a chart, table, or diagnostic output.")
    body = slide(
        "Stacked visual regions",
        layout_two_rows(top, bottom; ratio = :short_tall, height = "30rem")
    )
    return (; body)
end

function dashboard_page(::Session, _)
    dashboard = layout_quad(
        panel("Input", md"Primary input or assumption."),
        panel("State", md"Current operating state."),
        panel("Metric", md"A compact result or KPI."),
        panel("Interpretation", md"What the audience should conclude.");
        height = "30rem"
    )
    return (; body = slide("Four-panel dashboard", dashboard))
end

function image_page(::Session, _)
    image = webpart(
        local_image(
            EXAMPLE_IMAGE;
            alt = "Example diagram made from simple geometric elements",
            caption = "Presentation-local assets live beside the content directory."
        );
        kind = :panel,
        class = "lc-template-image-panel"
    )
    explanation = panel("Local assets", md"""
    Resolve assets from `@__DIR__` and pass their absolute path to
    `local_image`. The publisher remains unaware of presentation-specific
    images.
    """)
    body = slide(
        "Images and captions",
        layout_sidebar_main(explanation, image; height = "28rem")
    )
    return (; body)
end

const DECK = deck_descriptor(
    id = "starter-deck",
    group = "Authoring examples",
    title = "Starter deck",
    order = 10,
    render = true,
    setup = setup,
    resources = (PREFLIGHT_RESOURCE,),
    pages = (
        deck_page("Prose and mathematics"; id = "prose", build = prose_page),
        deck_page("Two columns"; id = "two-columns", build = two_columns_page),
        deck_page("Top and split"; id = "top-split", build = top_split_page),
        deck_page("Three columns"; id = "three-columns", build = three_columns_page),
        deck_page("Controls and plot"; id = "controls-plot", build = controls_plot_page),
        deck_page("Split with footer"; id = "split-bottom", build = split_bottom_page),
        deck_page("Two rows"; id = "two-rows", build = two_rows_page),
        deck_page("Dashboard"; id = "dashboard", build = dashboard_page),
        deck_page("Images"; id = "images", build = image_page)
    )
)

end

StarterDeck.DECK
