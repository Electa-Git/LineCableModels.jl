module PageAuthoring

using Bonito
using Markdown
using WGLMakie

export article
export action_row
export color_key
export control
export deck_descriptor
export deck_page
export diagnostic
export local_image
export layout_controls_plot
export layout_quad
export layout_sidebar_main
export layout_single
export layout_split_bottom
export layout_three_columns
export layout_top_split
export layout_two_columns
export layout_two_rows
export logo_row
export preparation_failed
export preparation_ready
export preparation_state
export preparation_status
export preflight_failed
export preflight_ready
export preflight_resource
export preflight_state
export preflight_summary
export prose
export reset_selectors!
export set_preparation!
export set_preflight!
export slide
export status_line
export title_canvas
export value_list
export webgrid
export webpart
export wgl_figure

const AREA_NAME = r"^[A-Za-z_][A-Za-z0-9_-]*$"
const WEBPART_KINDS = (:plain, :panel, :controls, :plot)
const COLUMN_RATIOS = Dict(
    :equal => "repeat(2, minmax(0, 1fr))",
    :narrow_wide => "minmax(14rem, 2fr) minmax(0, 5fr)",
    :wide_narrow => "minmax(0, 5fr) minmax(14rem, 2fr)"
)
const ROW_RATIOS = Dict(
    :equal => "repeat(2, minmax(0, 1fr))",
    :short_tall => "minmax(0, 2fr) minmax(0, 5fr)",
    :tall_short => "minmax(0, 5fr) minmax(0, 2fr)"
)

"""
    wgl_figure(figure; resize_to=:parent, screen_config...)

Embed a WGLMakie figure without exposing WGLMakie's per-canvas loading spinner.
The publisher's preflight UI owns preparation feedback; the hidden spinner node
is retained only as WGLMakie's documented removable initialization sentinel.
"""
function wgl_figure(figure; resize_to = :parent, screen_config...)
    haskey(screen_config, :spinner) && throw(ArgumentError(
        "wgl_figure owns the WGLMakie spinner; use the publisher preflight UI for loading feedback"
    ))
    silent_spinner = DOM.div(;
        class = "wglmakie-spinner",
        hidden = true,
        ariaHidden = "true"
    )
    return WGLMakie.WithConfig(
        figure;
        resize_to,
        spinner = silent_spinner,
        screen_config...
    )
end

classes(parts...) = join(filter(!isempty, string.(parts)), " ")
children(content) = content isa Tuple ? content : (content,)

"""Arrange related action buttons in one responsive row."""
function action_row(actions...; class::AbstractString = "")
    return DOM.div(
        actions...;
        class = classes("lc-action-row", class)
    )
end

"""
    reset_selectors!(selector => default, ...)

Restore Bonito selectors without notifying observers whose value is already at
its declared default. Existing reactive pipelines remain responsible for any
resulting visual or numerical update.
"""
function reset_selectors!(selector_defaults::Pair...)
    for (selector, default) in selector_defaults
        observable = hasproperty(selector, :value) ? selector.value : selector
        isequal(observable[], default) || (observable[] = default)
    end
    return nothing
end

function with_id(builder, id, content...; attributes...)
    if isnothing(id)
        return builder(content...; attributes...)
    end
    return builder(content...; id = id, attributes...)
end

"""
    prose(content; class="")

Mark content as Julia/Bonito-authored prose. Julia Markdown converts its LaTeX
nodes directly to Bonito KaTeX components.
"""
function prose(content; class::AbstractString = "")
    return DOM.div(
        content;
        class = classes("lc-markdown", "lc-bonito-markdown", class),
        dataMathRenderer = "bonito"
    )
end

"""
    article(content; hint=nothing, class="")

Place long-form Markdown in the narrow documentation reading column.
"""
function article(content; hint = nothing, class::AbstractString = "")
    hint_node = isnothing(hint) ? () :
                (DOM.p(hint; class = "lc-page-hint"),)
    return DOM.article(
        prose(content),
        hint_node...;
        class = classes("lc-article", class)
    )
end

"""
    slide(title, content...; lede=nothing)

Create the standard page heading and return a body tuple accepted by the
showcase page-discovery renderer.
"""
function slide(title::AbstractString, content...; lede = nothing)
    lede_node = isnothing(lede) ? () :
                (prose(lede; class = "lc-page-lede"),)
    heading = DOM.div(
        DOM.h1(title),
        lede_node...;
        class = "lc-page-heading"
    )
    return (heading, content...)
end

"""
    title_canvas(title; subtitle=nothing, event=nothing, logo=nothing,
                 footer=nothing, class="")

Create an unpanelled, full-page title composition for opening slides. Content
stays in one visual canvas while the shared theme owns its responsive scale and
spacing.
"""
function title_canvas(
        title::AbstractString;
        subtitle = nothing,
        event = nothing,
        logo = nothing,
        footer = nothing,
        class::AbstractString = ""
)
    brand = isnothing(logo) ? () :
            (DOM.div(logo; class = "lc-title-brand"),)
    event_line = isnothing(event) ? () :
                 (DOM.p(event; class = "lc-title-event"),)
    subtitle_line = isnothing(subtitle) ? () :
                    (DOM.p(subtitle; class = "lc-title-subtitle"),)
    footer_block = isnothing(footer) ? () :
                   (DOM.footer(children(footer)...; class = "lc-title-footer"),)
    return DOM.div(
        brand...,
        DOM.div(
            event_line...,
            DOM.h1(title; class = "lc-title-heading"),
            subtitle_line...;
            class = "lc-title-copy"
        ),
        footer_block...;
        class = classes("lc-title-canvas", class)
    )
end

function area_token(value)
    if isnothing(value) || value === Symbol(".")
        return "."
    elseif value isa Symbol || value isa AbstractString
        token = string(value)
        occursin(AREA_NAME, token) || throw(ArgumentError(
            "invalid grid-area name $(repr(token)); use letters, digits, underscores, or hyphens"
        ))
        return token
    end
    throw(ArgumentError(
        "grid areas must contain symbols, strings, `nothing`, or `Symbol(\".\")`"
    ))
end

function area_tokens(areas::AbstractMatrix)
    isempty(areas) && throw(ArgumentError("a web grid must contain at least one cell"))
    return map(area_token, areas)
end

function area_names(tokens::AbstractMatrix{<:AbstractString})
    names = String[]
    for row in axes(tokens, 1), column in axes(tokens, 2)

        name = tokens[row, column]
        name == "." || name in names || push!(names, name)
    end
    isempty(names) &&
        throw(ArgumentError("a web grid must contain at least one named area"))
    return names
end

function validate_rectangular_areas(tokens, names)
    for name in names
        positions = findall(==(name), tokens)
        rows = first.(Tuple.(positions))
        columns = last.(Tuple.(positions))
        row_range = minimum(rows):maximum(rows)
        column_range = minimum(columns):maximum(columns)
        all(tokens[row, column] == name for row in row_range, column in column_range) ||
            throw(ArgumentError(
                "grid area $(repr(name)) must form one rectangle when cells are merged"
            ))
    end
    return nothing
end

function css_areas(tokens)
    rows = map(axes(tokens, 1)) do row
        "'$(join(tokens[row, :], " "))'"
    end
    return join(rows, " ")
end

stacked_css_areas(names) = join(("'$name'" for name in names), " ")

"""
    webgrid(areas; cells...)

Lay out named webparts using a Julia matrix. Repeating a name merges those
cells. Named cells are supplied as keywords matching the matrix symbols.
"""
function webgrid(
        areas::AbstractMatrix;
        columns = "repeat($(size(areas, 2)), minmax(0, 1fr))",
        compact_columns = columns,
        rows = "repeat($(size(areas, 1)), minmax(0, 1fr))",
        gap = "1rem",
        width = "100%",
        height = "100%",
        stack::Bool = true,
        stack_rows = nothing,
        stack_height = "auto",
        id = nothing,
        class::AbstractString = "",
        cells...
)
    tokens = area_tokens(areas)
    names = area_names(tokens)
    validate_rectangular_areas(tokens, names)
    resolved_stack_rows = isnothing(stack_rows) ?
                          "repeat($(length(names)), auto)" : stack_rows

    expected = Set(Symbol.(names))
    provided = Set(keys(cells))
    missing = sort!(collect(setdiff(expected, provided)); by = string)
    unexpected = sort!(collect(setdiff(provided, expected)); by = string)
    isempty(missing) || throw(ArgumentError(
        "web grid is missing cells: $(join(string.(missing), ", "))"
    ))
    isempty(unexpected) || throw(ArgumentError(
        "web grid received cells absent from its area matrix: $(join(string.(unexpected), ", "))"
    ))

    cell_nodes = map(names) do name
        content = cells[Symbol(name)]
        DOM.div(
            children(content)...;
            class = classes("lc-grid-cell", "lc-grid-cell-$name"),
            style = Styles("grid-area" => name)
        )
    end
    style = Styles(
        "--lc-authored-columns" => columns,
        "--lc-compact-columns" => compact_columns,
        "--lc-authored-rows" => rows,
        "--lc-authored-areas" => css_areas(tokens),
        "--lc-authored-gap" => gap,
        "--lc-authored-height" => height,
        "--lc-stacked-areas" => stacked_css_areas(names),
        "--lc-stacked-rows" => resolved_stack_rows,
        "--lc-stacked-height" => stack_height
    )
    attributes = (
        class = classes("lc-canvas", class),
        dataStack = string(stack)
    )
    if !isnothing(id)
        attributes = merge(attributes, (; id))
    end
    return Grid(
        cell_nodes...;
        columns = "var(--lc-effective-columns)",
        rows = "var(--lc-effective-rows)",
        areas = "var(--lc-effective-areas)",
        gap = "var(--lc-effective-gap)",
        width = width,
        height = "var(--lc-effective-height)",
        align_items = "stretch",
        justify_items = "stretch",
        style,
        attributes...
    )
end

function layout_tracks(ratios, ratio::Symbol, axis::AbstractString)
    haskey(ratios, ratio) || throw(ArgumentError(
        "unsupported $axis ratio $(repr(ratio)); expected one of $(join(sort!(collect(keys(ratios)); by = string), ", "))"
    ))
    return ratios[ratio]
end

function master_class(name::AbstractString, class::AbstractString)
    classes("lc-master-layout", "lc-layout-$name", class)
end

"""Lay out one full-width authored region."""
function layout_single(
        content;
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    return webgrid(
        reshape([:main], 1, 1);
        columns = "minmax(0, 1fr)",
        rows = "minmax(0, 1fr)",
        gap,
        height,
        stack = false,
        id,
        class = master_class("single", class),
        main = content
    )
end

"""Lay out two responsive columns using a semantic width ratio."""
function layout_two_columns(
        left,
        right;
        ratio::Symbol = :equal,
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    columns = layout_tracks(COLUMN_RATIOS, ratio, "column")
    return webgrid(
        reshape([:left, :right], 1, 2);
        columns,
        rows = "minmax(0, 1fr)",
        gap,
        height,
        stack_rows = "repeat(2, auto)",
        id,
        class = master_class("two-columns", class),
        left,
        right
    )
end

"""Lay out three equal responsive columns."""
function layout_three_columns(
        left,
        center,
        right;
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    return webgrid(
        reshape([:left, :center, :right], 1, 3);
        columns = "repeat(3, minmax(0, 1fr))",
        rows = "minmax(0, 1fr)",
        gap,
        height,
        stack_rows = "repeat(3, auto)",
        id,
        class = master_class("three-columns", class),
        left,
        center,
        right
    )
end

"""Lay out a full-width top region over two responsive columns."""
function layout_top_split(
        top,
        left,
        right;
        ratio::Symbol = :equal,
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    columns = layout_tracks(COLUMN_RATIOS, ratio, "column")
    return webgrid(
        [:top :top; :left :right];
        columns,
        rows = "auto minmax(0, 1fr)",
        gap,
        height,
        stack_rows = "repeat(3, auto)",
        id,
        class = master_class("top-split", class),
        top,
        left,
        right
    )
end

"""Lay out two responsive columns over a full-width bottom region."""
function layout_split_bottom(
        left,
        right,
        bottom;
        ratio::Symbol = :equal,
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    columns = layout_tracks(COLUMN_RATIOS, ratio, "column")
    return webgrid(
        [:left :right; :bottom :bottom];
        columns,
        rows = "minmax(0, 1fr) auto",
        gap,
        height,
        stack_rows = "repeat(3, auto)",
        id,
        class = master_class("split-bottom", class),
        left,
        right,
        bottom
    )
end

"""Lay out a narrow contextual sidebar beside a primary region."""
function layout_sidebar_main(
        sidebar,
        main;
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    return webgrid(
        reshape([:sidebar, :main], 1, 2);
        columns = COLUMN_RATIOS[:narrow_wide],
        rows = "minmax(0, 1fr)",
        gap,
        height,
        stack_rows = "repeat(2, auto)",
        id,
        class = master_class("sidebar-main", class),
        sidebar,
        main
    )
end

"""
Lay out controls beside a primary plot, with an optional footer. Set
`footer_placement = :under_plot` to keep the controls at full height while the
footer occupies only the main column.
"""
function layout_controls_plot(
        controls,
        plot;
        footer = nothing,
        footer_placement::Symbol = :shared,
        footer_rows = "minmax(0, 1fr) auto",
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    isnothing(footer) && return layout_two_columns(
        controls,
        plot;
        ratio = :narrow_wide,
        height,
        gap,
        id,
        class = classes("lc-layout-controls-plot", class)
    )
    areas = if footer_placement === :shared
        [:controls :plot; :footer :footer]
    elseif footer_placement === :under_plot
        [:controls :plot; :controls :footer]
    else
        throw(ArgumentError(
            "unsupported footer placement $(repr(footer_placement)); expected :shared or :under_plot"
        ))
    end
    return webgrid(
        areas;
        columns = COLUMN_RATIOS[:narrow_wide],
        rows = footer_rows,
        gap,
        height,
        stack_rows = "repeat(3, auto)",
        id,
        class = master_class("controls-plot", class),
        controls,
        plot,
        footer
    )
end

"""Lay out two responsive rows using a semantic height ratio."""
function layout_two_rows(
        top,
        bottom;
        ratio::Symbol = :equal,
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    rows = layout_tracks(ROW_RATIOS, ratio, "row")
    return webgrid(
        reshape([:top, :bottom], 2, 1);
        columns = "minmax(0, 1fr)",
        rows,
        gap,
        height,
        stack_rows = "repeat(2, auto)",
        id,
        class = master_class("two-rows", class),
        top,
        bottom
    )
end

"""Lay out a responsive two-by-two dashboard."""
function layout_quad(
        top_left,
        top_right,
        bottom_left,
        bottom_right;
        height = "100%",
        gap = "1rem",
        id = nothing,
        class::AbstractString = ""
)
    return webgrid(
        [:top_left :top_right; :bottom_left :bottom_right];
        columns = "repeat(2, minmax(0, 1fr))",
        rows = "repeat(2, minmax(0, 1fr))",
        gap,
        height,
        stack_rows = "repeat(4, auto)",
        id,
        class = master_class("quad", class),
        top_left,
        top_right,
        bottom_left,
        bottom_right
    )
end

"""
    webpart(content...; kind=:plain, title=nothing, meta=nothing)

Wrap content in one of the presentation's stable panel forms.
"""
function webpart(
        content...;
        kind::Symbol = :plain,
        title = nothing,
        meta = nothing,
        id = nothing,
        class::AbstractString = "",
        body_class::AbstractString = ""
)
    kind in WEBPART_KINDS || throw(ArgumentError(
        "unsupported webpart kind $(repr(kind)); expected one of $(join(WEBPART_KINDS, ", "))"
    ))

    if kind === :plot
        body = with_id(
            DOM.div,
            id,
            content...;
            class = classes(
                "lc-webpart-body",
                "lc-plot-host",
                isnothing(title) ? "lc-webpart lc-webpart-plot" : "",
                body_class
            )
        )
        isnothing(title) && return body
        caption_meta = isnothing(meta) ? () : (DOM.span(meta),)
        caption = DOM.figcaption(DOM.strong(title), caption_meta...)
        return DOM.figure(
            caption,
            body;
            class = classes(
                "lc-webpart",
                "lc-webpart-plot",
                "lc-application-figure",
                class
            )
        )
    end

    heading = if isnothing(title)
        ()
    else
        heading_content = isnothing(meta) ? (DOM.strong(title),) :
                          (DOM.strong(title), DOM.span(meta))
        (DOM.div(heading_content...; class = "lc-webpart-header"),)
    end
    if kind === :panel
        body = DOM.div(
            content...;
            class = classes("lc-webpart-body", body_class)
        )
        return with_id(
            DOM.div,
            id,
            heading...,
            body;
            class = classes(
                "lc-webpart",
                "lc-webpart-panel",
                isnothing(title) ? "" : "lc-webpart-has-header",
                class
            )
        )
    end
    kind_class = kind === :controls ? "lc-webpart-controls lc-controls" :
                 "lc-webpart-plain"
    return with_id(
        DOM.div,
        id,
        heading...,
        content...;
        class = classes("lc-webpart", kind_class, class)
    )
end

"""Create a consistently labelled native Bonito control."""
function control(
        label::AbstractString,
        widget;
        value = nothing,
        id = nothing,
        class::AbstractString = ""
)
    value_node = isnothing(value) ? () : (DOM.output(value),)
    heading = DOM.div(
        DOM.span(label),
        value_node...;
        class = "lc-control-heading"
    )
    return with_id(
        DOM.label,
        id,
        heading,
        widget;
        class = classes("lc-control", class)
    )
end

"""Create a definition-list panel from `label => value` pairs."""
function value_list(entries::Pair...; class::AbstractString = "")
    rows = map(entries) do entry
        DOM.div(DOM.dt(first(entry)), DOM.dd(last(entry)))
    end
    return DOM.dl(rows...; class = classes("lc-values", class))
end

"""Create a small color key from `label => swatch_class` pairs."""
function color_key(
        entries::Pair...;
        separator = "→",
        class::AbstractString = ""
)
    nodes = Any[]
    for (index, entry) in enumerate(entries)
        index == 1 || push!(nodes,
            DOM.span(separator; class = "lc-color-key-separator", ariaHidden = "true")
        )
        push!(nodes,
            DOM.span(;
                class = classes("lc-color-key-swatch", last(entry)),
                ariaHidden = "true"
            ),
            DOM.span(first(entry); class = "lc-color-key-label")
        )
    end
    return DOM.div(nodes...; class = classes("lc-color-key", class))
end

function status_line(content; class::AbstractString = "")
    DOM.p(content; class = classes("lc-status-line", class))
end

function preflight_snapshot(phase::Symbol, progress::Real, label::AbstractString)
    isfinite(progress) || throw(ArgumentError("preparation progress must be finite"))
    return (;
        phase,
        progress = clamp(Float64(progress), 0.0, 1.0),
        label = string(label)
    )
end

"""Create independent process-wide state for one opt-in preflight resource."""
function preflight_state(;
        phase::Symbol = :cold,
        progress::Real = 0.0,
        label::AbstractString = "Not activated"
)
    return Observable(preflight_snapshot(phase, progress, label))
end

"""Declare a process-wide resource that a presentation may activate explicitly."""
function preflight_resource(;
        id,
        title,
        activate::Function,
        state = preflight_state(),
        required::Bool = true
)
    resource_id = string(id)
    occursin(r"^[a-z0-9]+(?:-[a-z0-9]+)*$", resource_id) ||
        throw(ArgumentError(
            "preflight resource ids must contain lowercase words separated by hyphens"
        ))
    resource_title = string(title)
    isempty(resource_title) && throw(ArgumentError(
        "preflight resource titles must not be empty"
    ))
    snapshot = state[]
    all(hasproperty(snapshot, field) for field in (:phase, :progress, :label)) ||
        throw(ArgumentError(
            "preflight state must expose phase, progress, and label fields"
        ))
    return (;
        id = resource_id,
        title = resource_title,
        required,
        state,
        activate
    )
end

function set_preflight!(state, phase::Symbol, progress::Real, label::AbstractString)
    state[] = preflight_snapshot(phase, progress, label)
    return state[]
end

preflight_ready(state) = state[].phase === :hot
preflight_failed(state) = state[].phase === :error

"""Summarize required resources as one cold, preparing, hot, or error state."""
function preflight_summary(resources)
    required = filter(resource -> resource.required, collect(resources))
    isempty(required) && return preflight_snapshot(
        :hot,
        1.0,
        "No preflight resources required"
    )
    states = getindex.(getproperty.(required, :state))
    progress = sum(state.progress for state in states) / length(states)
    failure = findfirst(state -> state.phase === :error, states)
    !isnothing(failure) && return preflight_snapshot(
        :error,
        progress,
        states[failure].label
    )
    all(state -> state.phase === :hot, states) && return preflight_snapshot(
        :hot,
        1.0,
        "All required resources are hot"
    )
    active = findfirst(state -> state.phase ∉ (:cold, :hot), states)
    !isnothing(active) && return preflight_snapshot(
        :preparing,
        progress,
        states[active].label
    )
    return preflight_snapshot(:cold, progress, "Activation required")
end

const PREPARATION_STATE = preflight_state(
    phase = :cold,
    label = "Waiting for numerical preparation…"
)

function set_preparation!(phase::Symbol, progress::Real, label::AbstractString)
    return set_preflight!(PREPARATION_STATE, phase, progress, label)
end

function set_preparation!(state, phase::Symbol, progress::Real, label::AbstractString)
    return set_preflight!(state, phase, progress, label)
end

preparation_ready(state = PREPARATION_STATE) = preflight_ready(state)
preparation_failed(state = PREPARATION_STATE) = preflight_failed(state)
preparation_state() = PREPARATION_STATE

function preparation_value(session, state, transform)
    return isnothing(session) ? transform(state[]) :
           map(transform, session, state)
end

"""Create the shell-owned display for a real process-wide preflight state."""
function preparation_status(
        session::Union{Nothing, Session} = nothing;
        state = PREPARATION_STATE,
        title::AbstractString = "Preparing resources"
)
    label = preparation_value(session, state, snapshot -> snapshot.label)
    progress = preparation_value(
        session,
        state,
        snapshot -> string(round(Int, 100snapshot.progress))
    )
    progress_style = preparation_value(
        session,
        state,
        snapshot -> "width: $(100snapshot.progress)%"
    )
    phase = preparation_value(
        session,
        state,
        snapshot -> string(snapshot.phase)
    )
    busy = preparation_value(
        session,
        state,
        snapshot -> string(snapshot.phase ∉ (:hot, :error))
    )
    return DOM.div(
        DOM.div(
            DOM.strong(title),
            DOM.output(
                label;
                id = "lc-preparation-label",
                ariaLive = "polite"
            );
            class = "lc-preparation-heading"
        ),
        DOM.div(
            DOM.div(;
                class = "lc-preparation-progress-bar",
                style = progress_style
            );
            id = "lc-preparation-progress",
            class = "lc-preparation-progress",
            role = "progressbar",
            ariaLabel = "Preflight preparation progress",
            ariaValueMin = "0",
            ariaValueMax = "100",
            ariaValueNow = progress
        );
        id = "lc-preparation",
        class = "lc-preparation",
        dataPhase = phase,
        ariaBusy = busy
    )
end

function diagnostic(
        value;
        prefix = "",
        suffix = "",
        class::AbstractString = ""
)
    return DOM.p(
        prefix,
        DOM.code(value),
        suffix;
        class = classes("lc-diagnostic", class)
    )
end

"""Serve a local image through Bonito, with an optional caption."""
function local_image(
        path::AbstractString;
        alt::AbstractString = "",
        caption = nothing,
        class::AbstractString = ""
)
    image = DOM.img(;
        src = Asset(abspath(path)),
        alt,
        class = classes("lc-local-image", isnothing(caption) ? class : "")
    )
    isnothing(caption) && return image
    return DOM.figure(
        image,
        DOM.figcaption(caption);
        class = classes("lc-local-image-figure", class)
    )
end

"""Arrange institutional or partner marks in one responsive logo row."""
function logo_row(logos...; class::AbstractString = "")
    return DOM.div(logos...; class = classes("lc-logo-row", class))
end

"""
    deck_page(title; id, render=true, class="", build)

Declare one addressable page inside a deck. The build function receives the
current Bonito `Session` and the state returned by the deck's `setup` function.
"""
function deck_page(
        build::Function,
        title::AbstractString;
        id,
        render::Bool = true,
        class::AbstractString = ""
)
    return (;
        id = string(id),
        title = string(title),
        render,
        class,
        build
    )
end

function deck_page(
        title::AbstractString;
        id,
        render::Bool = true,
        class::AbstractString = "",
        build
)
    return deck_page(build, title; id, render, class)
end

"""
    deck_descriptor(; id, group, title, order, pages, ...)

Construct the discoverable descriptor returned by one deck file. `setup` runs
once per browser session; every rendered page in the deck receives its result.
"""
function deck_descriptor(;
        id,
        group,
        title,
        order,
        render::Bool = true,
        class::AbstractString = "",
        setup::Function = _ -> nothing,
        resources = (),
        pages,
        extras...
)
    return (;
        id = string(id),
        group = string(group),
        title = string(title),
        order,
        render,
        class,
        setup,
        resources = Tuple(resources),
        pages = Tuple(pages),
        extras...
    )
end

end
