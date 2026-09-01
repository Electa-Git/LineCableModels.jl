abstract type AbstractToolbarItem end

struct ToolbarButton <: AbstractToolbarItem
    action::Symbol
    icon::Any
    label::Union{Nothing,String}
    tooltip::String
end

struct ToolbarDropdown <: AbstractToolbarItem
    action::Symbol
    icon::Any
    label::String
    options::Vector{Pair{Symbol,String}}
    selected::Symbol
    tooltip::String
end

struct ToolbarSeparator <: AbstractToolbarItem end

struct ToolbarEvent
    namespace::Symbol
    action::Symbol
    value::Union{Nothing,Symbol}
end

struct ToolbarBinding{F}
    namespace::Symbol
    callback::F
end

struct Toolbar
    dom::Any
    events::Observable{Union{Nothing,ToolbarEvent}}
end

const TOOLBAR_NAME_PATTERN = r"^[a-z][a-z0-9_]*$"
const TOOLBAR_ORIENTATIONS = (:horizontal, :vertical)
const TOOLBAR_SIZES = (:small, :medium, :large)

function valid_toolbar_name(name::Symbol, kind::AbstractString)
    occursin(TOOLBAR_NAME_PATTERN, string(name)) || throw(ArgumentError(
        "$kind must use lowercase snake_case and start with a letter: $name"
    ))
    return name
end

function ToolbarButton(
        action::Symbol;
        icon,
        label::Union{Nothing,AbstractString}=nothing,
        tooltip::Union{Nothing,AbstractString}=nothing
    )
    valid_toolbar_name(action, "toolbar action")
    resolved_label = isnothing(label) ? nothing : string(label)
    resolved_tooltip = if isnothing(tooltip)
        isnothing(resolved_label) ? replace(string(action), '_' => ' ') : resolved_label
    else
        string(tooltip)
    end
    return ToolbarButton(action, icon, resolved_label, resolved_tooltip)
end

function ToolbarDropdown(
        action::Symbol,
        options;
        icon,
        label::AbstractString,
        selected::Union{Nothing,Symbol}=nothing,
        tooltip::Union{Nothing,AbstractString}=nothing
    )
    valid_toolbar_name(action, "toolbar action")
    converted = Pair{Symbol,String}[
        Symbol(first(option)) => string(last(option)) for option in options
    ]
    isempty(converted) && throw(ArgumentError("toolbar dropdown requires an option"))
    keys = first.(converted)
    allunique(keys) || throw(ArgumentError(
        "toolbar dropdown option names must be unique: $action"
    ))
    foreach(key -> valid_toolbar_name(key, "toolbar option"), keys)
    resolved_selected = isnothing(selected) ? first(keys) : selected
    resolved_selected in keys || throw(ArgumentError(
        "selected option $resolved_selected is not declared for $action"
    ))
    resolved_tooltip = isnothing(tooltip) ? string(label) : string(tooltip)
    return ToolbarDropdown(
        action,
        icon,
        string(label),
        converted,
        resolved_selected,
        resolved_tooltip
    )
end

function ToolbarBinding(
        callback;
        namespace::Symbol
    )
    valid_toolbar_name(namespace, "toolbar namespace")
    return ToolbarBinding(namespace, callback)
end

function ToolbarBinding(
        owner::Module;
        namespace::Symbol=Symbol(lowercase(string(nameof(owner)))),
        handler::Symbol=:handle_toolbar
    )
    isdefined(owner, handler) || throw(ArgumentError(
        "module $(nameof(owner)) does not define $handler(event)"
    ))
    callback = getfield(owner, handler)
    probe = ToolbarEvent(namespace, :action, nothing)
    applicable(callback, probe) || throw(ArgumentError(
        "$(nameof(owner)).$handler does not accept a ToolbarEvent"
    ))
    return ToolbarBinding(callback; namespace)
end

function ToolbarBinding(
        owner::Module,
        state;
        namespace::Symbol=Symbol(lowercase(string(nameof(owner)))),
        handler::Symbol=:handle_toolbar
    )
    isdefined(owner, handler) || throw(ArgumentError(
        "module $(nameof(owner)) does not define $handler(state, event)"
    ))
    callback = getfield(owner, handler)
    probe = ToolbarEvent(namespace, :action, nothing)
    applicable(callback, state, probe) || throw(ArgumentError(
        "$(nameof(owner)).$handler does not accept the supplied state and a ToolbarEvent"
    ))
    return ToolbarBinding(event -> callback(state, event); namespace)
end

(binding::ToolbarBinding)(event::ToolbarEvent) = binding.callback(event)

toolbar_event_name(event::ToolbarEvent) = "$(event.namespace).$(event.action)"

function toolbar_icon(name::Symbol)
    children = if name == :reset
        (
            SVG.path(; d="M3 12a9 9 0 1 0 3-6.7"),
            SVG.path(; d="M3 3v6h6"),
        )
    elseif name == :fit
        (
            SVG.path(; d="M8 3H3v5"),
            SVG.path(; d="M16 3h5v5"),
            SVG.path(; d="M21 16v5h-5"),
            SVG.path(; d="M3 16v5h5"),
        )
    elseif name == :zoom_in
        (
            SVG.circle(; cx="11", cy="11", r="7"),
            SVG.path(; d="m20 20-4-4"),
            SVG.path(; d="M11 8v6M8 11h6"),
        )
    elseif name == :zoom_out
        (
            SVG.circle(; cx="11", cy="11", r="7"),
            SVG.path(; d="m20 20-4-4"),
            SVG.path(; d="M8 11h6"),
        )
    elseif name == :previous
        (
            SVG.path(; d="m15 18-6-6 6-6"),
        )
    elseif name == :next
        (
            SVG.path(; d="m9 18 6-6-6-6"),
        )
    elseif name == :play
        (
            SVG.path(; d="m8 5 11 7-11 7z"),
        )
    elseif name == :layers
        (
            SVG.path(; d="m12 2 9 5-9 5-9-5z"),
            SVG.path(; d="m3 12 9 5 9-5"),
            SVG.path(; d="m3 17 9 5 9-5"),
        )
    elseif name == :download
        (
            SVG.path(; d="M12 3v12"),
            SVG.path(; d="m7 10 5 5 5-5"),
            SVG.path(; d="M5 21h14"),
        )
    else
        throw(ArgumentError("unknown toolbar icon: $name"))
    end
    attributes = Dict{Symbol,Any}(
        Symbol("stroke-width") => "1.8",
        Symbol("stroke-linecap") => "round",
        Symbol("stroke-linejoin") => "round",
        Symbol("aria-hidden") => "true",
    )
    return SVG.svg(
        children...;
        attributes...,
        viewBox="0 0 24 24",
        fill="none",
        stroke="currentColor",
        focusable="false",
        class="lc-toolbar-icon"
    )
end

toolbar_icon(node) = node

function emit_toolbar_event!(events, binding, action, value=nothing)
    event = ToolbarEvent(binding.namespace, action, value)
    events[] = event
    binding(event)
    return event
end

function toolbar_button(item::ToolbarButton, events, binding)
    clicks = Observable(false)
    on(clicks) do _
        emit_toolbar_event!(events, binding, item.action)
        return nothing
    end
    icon_only = isnothing(item.label)
    content = if icon_only
        toolbar_icon(item.icon)
    else
        DOM.span(
            toolbar_icon(item.icon),
            DOM.span(item.label; class="lc-toolbar-button-label");
            class="lc-toolbar-button-content"
        )
    end
    class = icon_only ?
        "lc-toolbar-button lc-toolbar-button-icon-only" :
        "lc-toolbar-button"
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => item.tooltip,
        Symbol("data-toolbar-action") => string(item.action),
    )
    return DOM.button(
        content;
        attributes...,
        type="button",
        class,
        title=item.tooltip,
        onclick=js"event => $(clicks).notify(true)"
    )
end

function toolbar_dropdown(item::ToolbarDropdown, events, binding)
    selection = Observable(string(item.selected))
    on(selection) do selected
        emit_toolbar_event!(events, binding, item.action, Symbol(selected))
        return nothing
    end
    options = [
        DOM.option(
            label;
            value=string(key),
            selected=(key == item.selected)
        ) for (key, label) in item.options
    ]
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => item.tooltip,
        Symbol("data-toolbar-action") => string(item.action),
    )
    select = DOM.select(
        options...;
        attributes...,
        class="lc-toolbar-select",
        title=item.tooltip,
        onchange=js"event => $(selection).notify(event.currentTarget.value)"
    )
    return DOM.div(
        toolbar_icon(item.icon),
        DOM.span(item.label; class="lc-toolbar-dropdown-label"),
        select;
        class="lc-toolbar-dropdown"
    )
end

toolbar_item(::ToolbarSeparator, _, _) = DOM.span(
    ;
    class="lc-toolbar-separator",
    role="separator"
)
toolbar_item(item::ToolbarButton, events, binding) = toolbar_button(item, events, binding)
toolbar_item(item::ToolbarDropdown, events, binding) = toolbar_dropdown(item, events, binding)

function toolbar(
        items::AbstractVector{<:AbstractToolbarItem};
        binding::ToolbarBinding,
        orientation::Symbol=:horizontal,
        size::Symbol=:medium,
        label::AbstractString="Tools"
    )
    orientation in TOOLBAR_ORIENTATIONS || throw(ArgumentError(
        "toolbar orientation must be :horizontal or :vertical"
    ))
    size in TOOLBAR_SIZES || throw(ArgumentError(
        "toolbar size must be :small, :medium, or :large"
    ))
    actions = [item.action for item in items if !(item isa ToolbarSeparator)]
    allunique(actions) || throw(ArgumentError(
        "toolbar action names must be unique within $(binding.namespace)"
    ))
    events = Observable{Union{Nothing,ToolbarEvent}}(nothing)
    children = [toolbar_item(item, events, binding) for item in items]
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => string(label),
        Symbol("data-toolbar-namespace") => string(binding.namespace),
    )
    dom = DOM.div(
        children...;
        attributes...,
        class="lc-toolbar lc-toolbar-$orientation lc-toolbar-$size",
        role="toolbar"
    )
    return Toolbar(dom, events)
end

Bonito.jsrender(session::Session, toolbar::Toolbar) = Bonito.jsrender(session, toolbar.dom)
