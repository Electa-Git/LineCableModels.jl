module PageAuthoring

using Bonito
using Markdown

export article
export color_key
export control
export diagnostic
export local_image
export page_descriptor
export preparation_failed
export preparation_ready
export preparation_state
export preparation_status
export prose
export set_preparation!
export slide
export status_line
export value_list
export webgrid
export webpart

const AREA_NAME = r"^[A-Za-z_][A-Za-z0-9_-]*$"
const WEBPART_KINDS = (:plain, :panel, :controls, :plot)
const PREPARATION_STATE = Observable((;
    phase = :idle,
    progress = 0.0,
    label = "Waiting for numerical preparation…"
))

classes(parts...) = join(filter(!isempty, string.(parts)), " ")
children(content) = content isa Tuple ? content : (content,)

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
    kind_class = kind === :controls ? "lc-webpart-controls lc-controls" :
                 kind === :panel ? "lc-webpart-panel" : "lc-webpart-plain"
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

function set_preparation!(phase::Symbol, progress::Real, label::AbstractString)
    isfinite(progress) || throw(ArgumentError("preparation progress must be finite"))
    PREPARATION_STATE[] = (;
        phase,
        progress = clamp(Float64(progress), 0.0, 1.0),
        label = string(label)
    )
    return PREPARATION_STATE[]
end

preparation_ready() = PREPARATION_STATE[].phase === :ready
preparation_failed() = PREPARATION_STATE[].phase === :error
preparation_state() = PREPARATION_STATE

function preparation_value(session, transform)
    return isnothing(session) ? transform(PREPARATION_STATE[]) :
           map(transform, session, PREPARATION_STATE)
end

"""Create the shell-owned display for real process-wide preparation phases."""
function preparation_status(session::Union{Nothing, Session} = nothing)
    label = preparation_value(session, state -> state.label)
    progress = preparation_value(session, state -> string(round(Int, 100state.progress)))
    progress_style = preparation_value(
        session,
        state -> "width: $(100state.progress)%"
    )
    phase = preparation_value(session, state -> string(state.phase))
    busy = preparation_value(session, state -> string(state.phase ∉ (:ready, :error)))
    return DOM.div(
        DOM.div(
            DOM.strong("Preparing application case"),
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
            ariaLabel = "Application case preparation progress",
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

"""Construct the discoverable descriptor returned by a page file."""
function page_descriptor(;
        id,
        group,
        title,
        order,
        render::Bool = true,
        class::AbstractString = "",
        build,
        extras...
)
    return (;
        id = string(id),
        group = string(group),
        title = string(title),
        order,
        render,
        class,
        build,
        extras...
    )
end

end
