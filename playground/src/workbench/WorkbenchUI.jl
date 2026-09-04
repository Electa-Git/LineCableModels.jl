module WorkbenchUI

using Bonito
using UUIDs
using ..ComponentXRay
using ..Toolkit
using ..LineCableModelsPlayground: CONTROL_CONTRACT

export AbstractWorkbench,
    AbstractWorkbenchAction,
    AbstractWorkbenchView,
    ChoiceOption,
    Command,
    CommandChoice,
    Dock,
    DockTab,
    Identity,
    Inspector,
    Menu,
    MenuBar,
    NavGroup,
    NavItem,
    Sidebar,
    SplitPane,
    StatusBar,
    Toolbar,
    View,
    ViewStack,
    Workbench,
    compose,
    handle!,
    icon,
    initialize,
    workbench_app

abstract type AbstractWorkbench end
abstract type AbstractWorkbenchAction end
abstract type AbstractWorkbenchView end

const NAME_PATTERN = r"^[a-z][a-z0-9_]*$"
const BRAND_STYLES_PATH = normpath(
    joinpath(@__DIR__, "..", "..", "assets", "brand.css")
)
const WORKBENCH_STYLES_PATH = joinpath(@__DIR__, "workbench.css")
include_dependency(BRAND_STYLES_PATH)
include_dependency(WORKBENCH_STYLES_PATH)
const BRAND_STYLES = read(BRAND_STYLES_PATH, String)
const WORKBENCH_STYLES = read(WORKBENCH_STYLES_PATH, String)

function valid_name(name::Symbol, kind::AbstractString)
    occursin(NAME_PATTERN, string(name)) || throw(ArgumentError(
        "$kind must use lowercase snake_case and start with a letter: $name"
    ))
    return name
end

struct Identity{M}
    title::String
    subtitle::String
    mark::M
end

Identity(title::AbstractString, subtitle::AbstractString; mark=nothing) =
    Identity(string(title), string(subtitle), mark)

struct NavItem{A}
    id::Symbol
    label::String
    icon::Symbol
    action::A
    tooltip::String
    disabled::Bool
end

function NavItem(
        id::Symbol,
        label::AbstractString,
        action=nothing;
        icon::Symbol=:document,
        tooltip::AbstractString=label,
        disabled::Bool=false
    )
    valid_name(id, "navigation item")
    return NavItem(id, string(label), icon, action, string(tooltip), disabled)
end

struct NavGroup{T<:Tuple}
    label::String
    items::T
end

NavGroup(label::AbstractString, items::NavItem...) = NavGroup(string(label), items)

struct Sidebar{G<:Tuple,F,A}
    groups::G
    footer::F
    active::A
    initial_state::Symbol
end

function Sidebar(
        groups::NavGroup...;
        active,
        footer=nothing,
        initial_state::Symbol=:expanded
    )
    initial_state in (:expanded, :collapsed) || throw(ArgumentError(
        "sidebar initial_state must be :expanded or :collapsed"
    ))
    ids = (item.id for group in groups for item in group.items)
    allunique(ids) || throw(ArgumentError("navigation item ids must be unique"))
    return Sidebar(groups, footer, active, initial_state)
end

struct Command{A}
    action::A
    icon::Symbol
    label::Union{Nothing,String}
    tooltip::String
    disabled::Bool
end

function Command(
        action;
        icon::Symbol,
        label::Union{Nothing,AbstractString}=nothing,
        tooltip::Union{Nothing,AbstractString}=nothing,
        disabled::Bool=false
    )
    resolved_label = isnothing(label) ? nothing : string(label)
    resolved_tooltip = if isnothing(tooltip)
        isnothing(resolved_label) ? string(nameof(typeof(action))) : resolved_label
    else
        string(tooltip)
    end
    return Command(action, icon, resolved_label, resolved_tooltip, disabled)
end

struct ChoiceOption{A}
    id::Symbol
    label::String
    action::A
end

function ChoiceOption(id::Symbol, label::AbstractString, action)
    valid_name(id, "choice option")
    return ChoiceOption(id, string(label), action)
end

struct CommandChoice{O<:Tuple,S}
    label::String
    icon::Symbol
    options::O
    selected::S
    tooltip::String
    disabled::Bool
end

function CommandChoice(
        label::AbstractString,
        options::ChoiceOption...;
        selected,
        icon::Symbol=:layers,
        tooltip::AbstractString=label,
        disabled::Bool=false
    )
    isempty(options) && throw(ArgumentError("command choice requires an option"))
    allunique(option.id for option in options) || throw(ArgumentError(
        "command choice option ids must be unique"
    ))
    return CommandChoice(
        string(label),
        icon,
        options,
        selected,
        string(tooltip),
        disabled
    )
end

struct Menu{T<:Tuple}
    label::String
    items::T
end

Menu(label::AbstractString, items::Command...) = Menu(string(label), items)

struct MenuBar{T<:Tuple}
    menus::T
end

MenuBar(menus::Menu...) = MenuBar(menus)

struct Toolbar{T<:Tuple}
    items::T
    label::String
    size::Symbol
end

function Toolbar(
        items::Union{Command,CommandChoice}...;
        label::AbstractString="Workbench tools",
        size::Symbol=:medium
    )
    size in (:small, :medium, :large) || throw(ArgumentError(
        "toolbar size must be :small, :medium, or :large"
    ))
    return Toolbar(items, string(label), size)
end

struct View{C}
    id::Symbol
    title::String
    content::C
end

function View(id::Symbol, title::AbstractString, content)
    valid_name(id, "view")
    return View(id, string(title), content)
end

struct ViewStack{V<:Tuple,A}
    views::V
    active::A
end

struct SplitPane{A,B}
    first::A
    second::B
    orientation::Symbol
    ratio::Float64
    min_first::String
    min_second::String
    resizable::Bool
    responsive::Bool
end

function SplitPane(
        first,
        second;
        orientation::Symbol=:horizontal,
        ratio::Real=0.68,
        min_first::AbstractString="12rem",
        min_second::AbstractString="12rem",
        resizable::Bool=true,
        responsive::Bool=true
    )
    orientation in (:horizontal, :vertical) || throw(ArgumentError(
        "split orientation must be :horizontal or :vertical"
    ))
    0.15 <= ratio <= 0.85 || throw(ArgumentError(
        "split ratio must be between 0.15 and 0.85"
    ))
    return SplitPane(
        first,
        second,
        orientation,
        Float64(ratio),
        string(min_first),
        string(min_second),
        resizable,
        responsive
    )
end

function ViewStack(views::View...; active)
    isempty(views) && throw(ArgumentError("view stack requires a view"))
    allunique(view.id for view in views) || throw(ArgumentError(
        "view ids must be unique"
    ))
    return ViewStack(views, active)
end

struct Inspector{C}
    title::String
    subtitle::String
    content::C
end

Inspector(title::AbstractString, content; subtitle::AbstractString="Properties") =
    Inspector(string(title), string(subtitle), content)

struct DockTab{C}
    id::Symbol
    label::String
    content::C
end

function DockTab(id::Symbol, label::AbstractString, content)
    valid_name(id, "dock tab")
    return DockTab(id, string(label), content)
end

struct Dock{T<:Tuple,A}
    tabs::T
    active::A
    initial_state::Symbol
end

function Dock(tabs::DockTab...; active, initial_state::Symbol=:expanded)
    isempty(tabs) && throw(ArgumentError("dock requires a tab"))
    allunique(tab.id for tab in tabs) || throw(ArgumentError(
        "dock tab ids must be unique"
    ))
    initial_state in (:expanded, :collapsed) || throw(ArgumentError(
        "dock initial_state must be :expanded or :collapsed"
    ))
    return Dock(tabs, active, initial_state)
end

struct StatusBar{L,R}
    left::L
    right::R
end

StatusBar(left; right=nothing) = StatusBar(left, right)

struct Workbench{I,M,T,N,W,P,D,S}
    namespace::Symbol
    identity::I
    menu::M
    toolbar::T
    navigation::N
    workspace::W
    inspector::P
    output::D
    status::S
end

function Workbench(;
        namespace::Symbol,
        identity,
        navigation,
        workspace,
        menu=MenuBar(),
        toolbar=Toolbar(),
        inspector=nothing,
        output=nothing,
        status=nothing
    )
    valid_name(namespace, "workbench namespace")
    return Workbench(
        namespace,
        identity,
        menu,
        toolbar,
        navigation,
        workspace,
        inspector,
        output,
        status
    )
end

"""Create session-local state for an application."""
function initialize end

"""Describe a workbench from its application and session-local state."""
function compose end

"""Handle one typed workbench action."""
function handle! end

struct Runtime{A,S,W,X}
    application::A
    state::S
    workbench::W
    xray::X
end

function workbench_app(
        application::AbstractWorkbench;
        title::AbstractString="LineCableModels workbench",
        xray::ComponentXRay.XRayPolicy=ComponentXRay.XRayPolicy()
    )
    return App(; title=string(title)) do session
        state = initialize(application, session)
        workbench = compose(application, state)
        workbench isa Workbench || throw(ArgumentError(
            "compose must return a Workbench"
        ))
        return Runtime(application, state, workbench, xray)
    end
end

action_name(action) = action === nothing ? "none" : string(nameof(typeof(action)))

source(line) = ComponentXRay.source_reference(@__MODULE__, @__FILE__, line)

function ComponentXRay.inspection(identity::Identity)
    return ComponentXRay.ComponentInspection(
        identity;
        name="Identity",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:title, identity.title),
            ComponentXRay.PropertyInspection(:subtitle, identity.subtitle),
            ComponentXRay.PropertyInspection(:mark, isnothing(identity.mark) ? :default : :custom),
        ],
        css_scopes=[".lc-wb-identity"]
    )
end

function ComponentXRay.inspection(item::NavItem)
    return ComponentXRay.ComponentInspection(
        item;
        name="NavItem",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:id, item.id),
            ComponentXRay.PropertyInspection(:label, item.label),
            ComponentXRay.PropertyInspection(:icon, item.icon),
            ComponentXRay.PropertyInspection(:tooltip, item.tooltip),
            ComponentXRay.PropertyInspection(:disabled, item.disabled),
        ],
        actions=[ComponentXRay.ActionInspection(
            :click,
            action_name(item.action),
            item.action;
            disabled=item.disabled
        )],
        css_scopes=[".lc-wb-nav-item", ".lc-wb-nav-icon", ".lc-wb-nav-label"]
    )
end

function ComponentXRay.inspection(group::NavGroup)
    return ComponentXRay.ComponentInspection(
        group;
        name="NavGroup",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:label, group.label),
            ComponentXRay.PropertyInspection(:items, length(group.items)),
        ],
        css_scopes=[".lc-wb-nav-group", ".lc-wb-nav-heading", ".lc-wb-nav-items"]
    )
end

function ComponentXRay.inspection(sidebar::Sidebar)
    return ComponentXRay.ComponentInspection(
        sidebar;
        name="Sidebar",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:groups, length(sidebar.groups)),
            ComponentXRay.PropertyInspection(:initial_state, sidebar.initial_state),
            ComponentXRay.PropertyInspection(:footer, !isnothing(sidebar.footer)),
        ],
        bindings=[ComponentXRay.BindingInspection(
            :active,
            sidebar.active;
            notes="session-local navigation selection"
        )],
        css_scopes=[".lc-wb-sidebar", ".lc-wb-navigation"]
    )
end

function ComponentXRay.inspection(command::Command)
    label = isnothing(command.label) ? command.tooltip : command.label
    return ComponentXRay.ComponentInspection(
        command;
        name="Command",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:label, label),
            ComponentXRay.PropertyInspection(:icon, command.icon),
            ComponentXRay.PropertyInspection(:tooltip, command.tooltip),
            ComponentXRay.PropertyInspection(:disabled, command.disabled),
        ],
        actions=[ComponentXRay.ActionInspection(
            :click,
            action_name(command.action),
            command.action;
            disabled=command.disabled
        )],
        css_scopes=[".lc-wb-command", ".lc-wb-menu-command"]
    )
end

function ComponentXRay.inspection(choice::CommandChoice)
    actions = [
        ComponentXRay.ActionInspection(
            :change,
            string(option.id),
            option.action;
            disabled=choice.disabled
        ) for option in choice.options
    ]
    return ComponentXRay.ComponentInspection(
        choice;
        name="CommandChoice",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:label, choice.label),
            ComponentXRay.PropertyInspection(:icon, choice.icon),
            ComponentXRay.PropertyInspection(:options, length(choice.options)),
            ComponentXRay.PropertyInspection(:tooltip, choice.tooltip),
            ComponentXRay.PropertyInspection(:disabled, choice.disabled),
        ],
        bindings=[ComponentXRay.BindingInspection(
            :selected,
            choice.selected;
            notes="live choice selection"
        )],
        actions,
        css_scopes=[
            ".lc-wb-choice",
            ".lc-wb-choice-select",
            ".lc-wb-choice-label",
            ".lc-control-select",
        ]
    )
end

function ComponentXRay.inspection(menu::Menu)
    return ComponentXRay.ComponentInspection(
        menu;
        name="Menu",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:label, menu.label),
            ComponentXRay.PropertyInspection(:items, length(menu.items)),
        ],
        actions=[ComponentXRay.ActionInspection(:toggle, "native details", nothing)],
        css_scopes=[".lc-wb-menu", ".lc-wb-menu-popover"]
    )
end

function ComponentXRay.inspection(menubar::MenuBar)
    return ComponentXRay.ComponentInspection(
        menubar;
        name="MenuBar",
        source=source(@__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:menus, length(menubar.menus))],
        css_scopes=[
            ".lc-wb-menubar",
            ".lc-wb-menu-groups",
            ".lc-wb-menubar-tools",
            ".lc-control-select",
        ]
    )
end

function ComponentXRay.inspection(toolbar::Toolbar)
    return ComponentXRay.ComponentInspection(
        toolbar;
        name="Toolbar",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:label, toolbar.label),
            ComponentXRay.PropertyInspection(:size, toolbar.size),
            ComponentXRay.PropertyInspection(:items, length(toolbar.items)),
        ],
        css_scopes=[".lc-wb-toolbar"]
    )
end

function ComponentXRay.inspection(view::View)
    return ComponentXRay.ComponentInspection(
        view;
        name="View",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:id, view.id),
            ComponentXRay.PropertyInspection(:title, view.title),
            ComponentXRay.PropertyInspection(:content, ComponentXRay.short_type(view.content)),
        ],
        css_scopes=[".lc-wb-view"]
    )
end

function ComponentXRay.inspection(stack::ViewStack)
    return ComponentXRay.ComponentInspection(
        stack;
        name="ViewStack",
        source=source(@__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:views, length(stack.views))],
        bindings=[ComponentXRay.BindingInspection(
            :active,
            stack.active;
            notes="mounted views remain persistent"
        )],
        css_scopes=[".lc-wb-workspace", ".lc-wb-view"]
    )
end

function ComponentXRay.inspection(split::SplitPane)
    return ComponentXRay.ComponentInspection(
        split;
        name="SplitPane",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:orientation, split.orientation),
            ComponentXRay.PropertyInspection(:ratio, split.ratio),
            ComponentXRay.PropertyInspection(:min_first, split.min_first),
            ComponentXRay.PropertyInspection(:min_second, split.min_second),
            ComponentXRay.PropertyInspection(:resizable, split.resizable),
            ComponentXRay.PropertyInspection(:responsive, split.responsive),
        ],
        actions=[ComponentXRay.ActionInspection(:pointerdrag, "resize split", nothing)],
        css_scopes=[".lc-wb-split", ".lc-wb-splitter", ".lc-wb-split-region"]
    )
end

function ComponentXRay.inspection(inspector::Inspector)
    return ComponentXRay.ComponentInspection(
        inspector;
        name="Inspector",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:title, inspector.title),
            ComponentXRay.PropertyInspection(:subtitle, inspector.subtitle),
            ComponentXRay.PropertyInspection(:content, ComponentXRay.short_type(inspector.content)),
        ],
        css_scopes=[".lc-wb-inspector"]
    )
end

function ComponentXRay.inspection(tab::DockTab)
    return ComponentXRay.ComponentInspection(
        tab;
        name="DockTab",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:id, tab.id),
            ComponentXRay.PropertyInspection(:label, tab.label),
            ComponentXRay.PropertyInspection(:content, ComponentXRay.short_type(tab.content)),
        ],
        css_scopes=[".lc-wb-dock-panel"]
    )
end

function ComponentXRay.inspection(dock::Dock)
    return ComponentXRay.ComponentInspection(
        dock;
        name="Dock",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:tabs, length(dock.tabs)),
            ComponentXRay.PropertyInspection(:initial_state, dock.initial_state),
        ],
        bindings=[ComponentXRay.BindingInspection(
            :active,
            dock.active;
            notes="session-local output selection"
        )],
        actions=[ComponentXRay.ActionInspection(:pointerdrag, "resize dock", nothing)],
        css_scopes=[".lc-wb-dock", ".lc-wb-dock-header", ".lc-wb-dock-content"]
    )
end

function ComponentXRay.inspection(status::StatusBar)
    return ComponentXRay.ComponentInspection(
        status;
        name="StatusBar",
        source=source(@__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:right_content, !isnothing(status.right))],
        css_scopes=[".lc-wb-statusbar"]
    )
end

function ComponentXRay.inspection(workbench::Workbench)
    return ComponentXRay.ComponentInspection(
        workbench;
        name="Workbench",
        source=source(@__LINE__),
        parameters=[
            ComponentXRay.PropertyInspection(:namespace, workbench.namespace),
            ComponentXRay.PropertyInspection(:inspector, !isnothing(workbench.inspector)),
            ComponentXRay.PropertyInspection(:output, !isnothing(workbench.output)),
            ComponentXRay.PropertyInspection(:status, !isnothing(workbench.status)),
        ],
        css_scopes=[".lc-wb-shell", ".lc-wb-page"],
        notes=["The X-ray reports semantic workbench components, not anonymous DOM descendants."]
    )
end

function icon(name::Symbol; class::AbstractString="lc-wb-icon")
    children = if name == :workbench
        (
            SVG.rect(; x="3", y="3", width="7", height="7"),
            SVG.rect(; x="14", y="3", width="7", height="7"),
            SVG.rect(; x="3", y="14", width="7", height="7"),
            SVG.rect(; x="14", y="14", width="7", height="7"),
        )
    elseif name == :geometry || name == :cable
        (
            SVG.circle(; cx="12", cy="12", r="8.5"),
            SVG.circle(; cx="12", cy="8.5", r="2.1"),
            SVG.circle(; cx="8.7", cy="14.2", r="2.1"),
            SVG.circle(; cx="15.3", cy="14.2", r="2.1"),
        )
    elseif name == :chart
        (
            SVG.path(; d="M3 19h18M5 17V5m1 10 4-5 3 2 5-7"),
        )
    elseif name == :document
        (
            SVG.path(; d="M6 3h9l4 4v14H6zM14 3v5h5M9 13h6M9 17h6"),
        )
    elseif name == :archive
        (
            SVG.path(; d="M4 4h16v5H4zM5.5 9v11h13V9M9 13h6"),
        )
    elseif name == :terminal
        (
            SVG.path(; d="m5 7 4 5-4 5M12 17h7"),
        )
    elseif name == :cloud
        (
            SVG.path(; d="M7 18h10a4 4 0 0 0 .5-8A6 6 0 0 0 6 9a4.5 4.5 0 0 0 1 9"),
        )
    elseif name == :fit
        (
            SVG.path(; d="M8 3H3v5M16 3h5v5M21 16v5h-5M3 16v5h5"),
        )
    elseif name == :reset
        (
            SVG.path(; d="M3 12a9 9 0 1 0 3-6.7M3 3v6h6"),
        )
    elseif name == :layers
        (
            SVG.path(; d="m12 2 9 5-9 5-9-5zM3 12l9 5 9-5M3 17l9 5 9-5"),
        )
    elseif name == :download
        (
            SVG.path(; d="M12 3v12m-5-5 5 5 5-5M5 21h14"),
        )
    elseif name == :play
        (SVG.path(; d="m8 5 11 7-11 7z"),)
    elseif name == :warning
        (
            SVG.path(; d="M12 3 2.5 20h19zM12 9v5M12 17h.01"),
        )
    elseif name == :settings
        (
            SVG.circle(; cx="12", cy="12", r="3"),
            SVG.path(; d="M19 12a7 7 0 0 0-.1-1l2-1.5-2-3.4-2.3 1A7 7 0 0 0 15 6l-.3-2.5h-4L10.4 6a7 7 0 0 0-1.7 1.1l-2.3-1-2 3.4 2 1.5a7 7 0 0 0 0 2l-2 1.5 2 3.4 2.3-1A7 7 0 0 0 10.4 18l.3 2.5h4L15 18a7 7 0 0 0 1.7-1.1l2.3 1 2-3.4-2-1.5a7 7 0 0 0 0-1z"),
        )
    elseif name == :menu
        (SVG.path(; d="M4 7h16M4 12h16M4 17h16"),)
    elseif name == :chevron_left
        (SVG.path(; d="m15 18-6-6 6-6"),)
    elseif name == :chevron_up
        (SVG.path(; d="m6 15 6-6 6 6"),)
    else
        throw(ArgumentError("unknown workbench icon: $name"))
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
        class=string(class)
    )
end

active_class(session, observable, id, base) = map(session, observable) do current
    current == id ? "$base is-active" : base
end

function action_button(session, command::Command, dispatch; class="lc-wb-command")
    trigger = Observable(false)
    on(session, trigger) do _
        command.disabled || dispatch(command.action)
        return nothing
    end
    label = isnothing(command.label) ? nothing : DOM.span(
        command.label;
        class="lc-wb-command-label"
    )
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => command.tooltip,
        Symbol("data-tooltip") => command.tooltip,
    )
    node = DOM.button(
        icon(command.icon),
        label;
        attributes...,
        class,
        type="button",
        title=command.tooltip,
        disabled=command.disabled,
        onclick=js"event => $(trigger).notify(true)"
    )
    return ComponentXRay.instrument(session, node, command)
end

function choice_control(session, choice::CommandChoice, dispatch)
    changed = Observable(string(choice.selected[]))
    on(session, changed) do selected
        option = findfirst(candidate -> string(candidate.id) == selected, choice.options)
        isnothing(option) || dispatch(choice.options[option].action)
        return nothing
    end
    options = [
        DOM.option(
            option.label;
            value=string(option.id),
            selected=(option.id == choice.selected[])
        ) for option in choice.options
    ]
    select_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => choice.tooltip,
    )
    select = DOM.select(
        options...;
        select_attributes...,
        class="lc-control-select lc-wb-choice-select",
        disabled=choice.disabled,
        onchange=js"event => $(changed).notify(event.currentTarget.value)"
    )
    synchronize = js"""
    (() => {
        const select = $(select);
        const selected = $(choice.selected);
        selected.on(value => { select.value = String(value); });
    })();
    """
    node = DOM.div(
        icon(choice.icon),
        DOM.span(choice.label; class="lc-wb-choice-label"),
        select,
        DOM.script(synchronize);
        class="lc-wb-choice",
        title=choice.tooltip
    )
    return ComponentXRay.instrument(session, node, choice)
end

function render_menu(session, menu::Menu, dispatch)
    entries = map(menu.items) do command
        button = action_button(
            session,
            command,
            dispatch;
            class="lc-wb-menu-command"
        )
        return button
    end
    node = DOM.details(
        DOM.summary(menu.label),
        DOM.div(entries...; class="lc-wb-menu-popover");
        class="lc-wb-menu"
    )
    return ComponentXRay.instrument(session, node, menu)
end

function render_menubar(session, menubar::MenuBar, dispatch, namespace)
    theme_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Workbench color theme",
        Symbol("data-lc-wb-theme-selector") => "",
    )
    node = DOM.header(
        DOM.div(
            (render_menu(session, menu, dispatch) for menu in menubar.menus)...;
            class="lc-wb-menu-groups"
        ),
        DOM.div(
            DOM.span(namespace; class="lc-wb-menu-context"),
            DOM.label(
                DOM.span("Theme"),
                DOM.select(
                    DOM.option("System"; value="system"),
                    DOM.option("Dark"; value="dark", selected=true),
                    DOM.option("Light"; value="light");
                    theme_attributes...,
                    class="lc-control-select lc-wb-theme-select"
                );
                class="lc-wb-theme-control"
            );
            class="lc-wb-menubar-tools"
        );
        class="lc-wb-menubar"
    )
    return ComponentXRay.instrument(session, node, menubar)
end

function render_toolbar(session, toolbar::Toolbar, dispatch)
    items = map(toolbar.items) do item
        item isa Command ?
            action_button(session, item, dispatch) :
            choice_control(session, item, dispatch)
    end
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => toolbar.label,
    )
    node = DOM.div(
        items...;
        attributes...,
        class="lc-wb-toolbar lc-wb-toolbar-$(toolbar.size)",
        role="toolbar"
    )
    return ComponentXRay.instrument(session, node, toolbar)
end

function render_nav_item(session, item::NavItem, active, dispatch)
    trigger = Observable(false)
    on(session, trigger) do _
        item.disabled || isnothing(item.action) || dispatch(item.action)
        return nothing
    end
    class = active_class(session, active, item.id, "lc-wb-nav-item")
    attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => item.tooltip,
        Symbol("aria-disabled") => string(item.disabled),
        Symbol("data-tooltip") => item.tooltip,
    )
    node = DOM.button(
        icon(item.icon; class="lc-wb-nav-icon"),
        DOM.span(item.label; class="lc-wb-nav-label"),
        DOM.span(; class="lc-wb-active-mark");
        attributes...,
        type="button",
        class,
        disabled=item.disabled,
        onclick=js"event => $(trigger).notify(true)"
    )
    return ComponentXRay.instrument(session, node, item)
end

function render_sidebar(session, identity::Identity, sidebar::Sidebar, dispatch)
    mark = isnothing(identity.mark) ? icon(:workbench; class="lc-wb-identity-mark") :
        identity.mark
    groups = map(sidebar.groups) do group
        node = DOM.section(
            DOM.h2(group.label; class="lc-wb-nav-heading"),
            DOM.div(
                (
                    render_nav_item(session, item, sidebar.active, dispatch)
                    for item in group.items
                )...;
                class="lc-wb-nav-items"
            );
            class="lc-wb-nav-group"
        )
        return ComponentXRay.instrument(session, node, group)
    end
    footer = isnothing(sidebar.footer) ? nothing : DOM.footer(
        sidebar.footer;
        class="lc-wb-sidebar-footer"
    )
    toggle_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Collapse navigation",
        Symbol("data-tooltip") => "Collapse navigation",
        Symbol("data-lc-wb-sidebar-toggle") => "",
    )
    navigation_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Workbench",
    )
    identity_node = DOM.div(
        mark,
        DOM.span(
            DOM.strong(identity.title),
            DOM.small(identity.subtitle);
            class="lc-wb-identity-copy"
        );
        class="lc-wb-identity"
    )
    identity_node = ComponentXRay.instrument(session, identity_node, identity)
    node = DOM.aside(
        DOM.header(
            identity_node,
            DOM.button(
                icon(:chevron_left);
                toggle_attributes...,
                type="button",
                class="lc-wb-sidebar-toggle"
            );
            class="lc-wb-sidebar-header"
        ),
        DOM.nav(
            groups...;
            navigation_attributes...,
            class="lc-wb-navigation"
        ),
        footer;
        class="lc-wb-sidebar"
    )
    return ComponentXRay.instrument(session, node, sidebar)
end

function render_viewstack(session, stack::ViewStack)
    views = map(stack.views) do view
        attributes = Dict{Symbol,Any}(
            Symbol("aria-label") => view.title,
            Symbol("data-view") => string(view.id),
        )
        node = DOM.section(
            view.content;
            attributes...,
            class=active_class(session, stack.active, view.id, "lc-wb-view"),
            role="region"
        )
        return ComponentXRay.instrument(session, node, view)
    end
    node = DOM.main(views...; class="lc-wb-workspace")
    return ComponentXRay.instrument(session, node, stack)
end

function split_script(root, orientation)
    return js"""
    (() => {
        const root = $(root);
        const splitter = root.querySelector(':scope > [data-lc-wb-splitter]');
        if (!splitter) return;

        let dragging = false;
        let firstSize = 0;
        let secondSize = 0;
        let start = 0;

        splitter.addEventListener('pointerdown', event => {
            dragging = true;
            const first = root.querySelector(':scope > .lc-wb-split-first');
            const second = root.querySelector(':scope > .lc-wb-split-second');
            const firstBounds = first.getBoundingClientRect();
            const secondBounds = second.getBoundingClientRect();
            firstSize = $(orientation == :horizontal) ? firstBounds.width : firstBounds.height;
            secondSize = $(orientation == :horizontal) ? secondBounds.width : secondBounds.height;
            start = $(orientation == :horizontal) ? event.clientX : event.clientY;
            splitter.setPointerCapture(event.pointerId);
            root.classList.add('is-resizing');
            event.preventDefault();
        });

        splitter.addEventListener('pointermove', event => {
            if (!dragging) return;
            const position = $(orientation == :horizontal) ? event.clientX : event.clientY;
            const delta = position - start;
            const total = firstSize + secondSize;
            const minimum = Math.min(140, total * 0.35);
            const nextFirst = Math.max(minimum, Math.min(total - minimum, firstSize + delta));
            const nextSecond = total - nextFirst;
            if ($(orientation == :horizontal)) {
                root.style.gridTemplateColumns = `${nextFirst}px 5px ${nextSecond}px`;
            } else {
                root.style.gridTemplateRows = `${nextFirst}px 5px ${nextSecond}px`;
            }
        });

        const stop = event => {
            if (!dragging) return;
            dragging = false;
            root.classList.remove('is-resizing');
            if (splitter.hasPointerCapture(event.pointerId)) {
                splitter.releasePointerCapture(event.pointerId);
            }
        };
        splitter.addEventListener('pointerup', stop);
        splitter.addEventListener('pointercancel', stop);
    })();
    """
end

function Bonito.jsrender(session::Session, split::SplitPane)
    divider = split.resizable ? "5px" : "1px"
    template = if split.orientation == :horizontal
        "minmax($(split.min_first), $(split.ratio)fr) $divider " *
        "minmax($(split.min_second), $(1 - split.ratio)fr)"
    else
        "minmax($(split.min_first), $(split.ratio)fr) $divider " *
        "minmax($(split.min_second), $(1 - split.ratio)fr)"
    end
    style = split.orientation == :horizontal ?
        "grid-template-columns: $template;" :
        "grid-template-rows: $template;"
    attributes = Dict{Symbol,Any}(
        Symbol("data-orientation") => string(split.orientation),
        Symbol("data-responsive") => string(split.responsive),
    )
    splitter_attributes = Dict{Symbol,Any}(
        Symbol("data-lc-wb-splitter") => "",
        Symbol("aria-label") => "Drag to resize split pane",
        Symbol("aria-orientation") => split.orientation == :horizontal ?
            "vertical" : "horizontal",
        :role => "separator",
    )
    class = "lc-wb-split lc-wb-split-$(split.orientation)" *
        (split.resizable ? "" : " is-fixed")
    splitter = DOM.div(
        DOM.span(split.orientation == :horizontal ? "↔" : "↕";
            Symbol("aria-hidden") => "true",
            class="lc-wb-splitter-handle"
        );
        splitter_attributes...,
        class="lc-wb-splitter"
    )
    root = DOM.div(
        DOM.div(split.first; class="lc-wb-split-region lc-wb-split-first"),
        splitter,
        DOM.div(split.second; class="lc-wb-split-region lc-wb-split-second");
        attributes...,
        class,
        style
    )
    script = split.resizable ? DOM.script(split_script(root, split.orientation)) : nothing
    instrumented = ComponentXRay.instrument(session, root, split)
    return Bonito.jsrender(
        session,
        DOM.div(instrumented, script; class="lc-wb-split-mount")
    )
end

function render_inspector(session, inspector::Inspector)
    node = DOM.aside(
        DOM.header(
            DOM.div(
                DOM.span(inspector.subtitle; class="lc-wb-panel-kicker"),
                DOM.h2(inspector.title)
            );
            class="lc-wb-inspector-header"
        ),
        DOM.div(inspector.content; class="lc-wb-inspector-content");
        class="lc-wb-inspector"
    )
    return ComponentXRay.instrument(session, node, inspector)
end

function render_dock(session, dock::Dock)
    tab_buttons = map(dock.tabs) do tab
        trigger = Observable(false)
        on(session, trigger) do _
            dock.active[] = tab.id
            return nothing
        end
        return DOM.button(
            tab.label;
            class=active_class(session, dock.active, tab.id, "lc-wb-dock-tab"),
            type="button",
            onclick=js"event => $(trigger).notify(true)"
        )
    end
    panels = map(dock.tabs) do tab
        attributes = Dict{Symbol,Any}(
            Symbol("aria-label") => tab.label,
        )
        node = DOM.section(
            tab.content;
            attributes...,
            class=active_class(session, dock.active, tab.id, "lc-wb-dock-panel"),
        )
        return ComponentXRay.instrument(session, node, tab)
    end
    splitter_attributes = Dict{Symbol,Any}(
        Symbol("data-lc-wb-dock-splitter") => "",
    )
    toggle_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => "Collapse output panel",
        Symbol("data-tooltip") => "Collapse output panel",
        Symbol("data-lc-wb-dock-toggle") => "",
    )
    node = DOM.section(
        DOM.div(; splitter_attributes..., class="lc-wb-dock-splitter"),
        DOM.header(
            DOM.div(tab_buttons...; class="lc-wb-dock-tabs", role="tablist"),
            DOM.button(
                icon(:chevron_up);
                toggle_attributes...,
                type="button",
                class="lc-wb-dock-toggle"
            );
            class="lc-wb-dock-header"
        ),
        DOM.div(panels...; class="lc-wb-dock-content");
        class="lc-wb-dock"
    )
    return ComponentXRay.instrument(session, node, dock)
end

function render_status(session, status::StatusBar)
    right = isnothing(status.right) ? nothing : DOM.div(
        status.right;
        class="lc-wb-status-right"
    )
    node = DOM.footer(
        DOM.div(status.left; class="lc-wb-status-left"),
        right;
        class="lc-wb-statusbar"
    )
    return ComponentXRay.instrument(session, node, status)
end

function intrinsic_script(root)
    return js"""
    (() => {
        const root = $(root);
        const sidebarToggle = root.querySelector('[data-lc-wb-sidebar-toggle]');
        const dockToggle = root.querySelector('[data-lc-wb-dock-toggle]');
        const splitter = root.querySelector('[data-lc-wb-dock-splitter]');
        const page = root.closest('.lc-wb-page');
        const themeSelector = root.querySelector('[data-lc-wb-theme-selector]');
        const systemTheme = window.matchMedia('(prefers-color-scheme: light)');
        const themeStorageKey = 'lcm.playground.theme';
        const legacyThemeStorageKey = 'lcm.workbench.theme';
        const menus = Array.from(root.querySelectorAll('.lc-wb-menu'));

        const validTheme = value => ['system', 'dark', 'light'].includes(value);
        const applyTheme = preference => {
            const selected = validTheme(preference) ? preference : 'dark';
            const resolved = selected === 'system'
                ? (systemTheme.matches ? 'light' : 'dark')
                : selected;
            page.dataset.theme = selected;
            page.dataset.resolvedTheme = resolved;
            document.documentElement.dataset.lcmTheme = selected;
            document.documentElement.dataset.lcmResolvedTheme = resolved;
            themeSelector.value = selected;
        };

        let initialTheme = 'dark';
        try {
            const storedTheme = window.localStorage.getItem(themeStorageKey)
                || window.localStorage.getItem(legacyThemeStorageKey);
            if (validTheme(storedTheme)) initialTheme = storedTheme;
        } catch (_) {
            // Storage may be disabled; the selector remains session-local.
        }
        applyTheme(initialTheme);
        themeSelector.addEventListener('change', event => {
            const preference = event.currentTarget.value;
            applyTheme(preference);
            try {
                window.localStorage.setItem(themeStorageKey, preference);
                window.localStorage.removeItem(legacyThemeStorageKey);
            } catch (_) {
                // Applying the selected theme does not depend on persistence.
            }
        });
        systemTheme.addEventListener('change', () => {
            if (page.dataset.theme === 'system') applyTheme('system');
        });
        window.addEventListener('storage', event => {
            if (![themeStorageKey, legacyThemeStorageKey].includes(event.key)) return;
            let storedTheme = null;
            try {
                storedTheme = window.localStorage.getItem(themeStorageKey)
                    || window.localStorage.getItem(legacyThemeStorageKey);
            } catch (_) {
                // Fall through to the deterministic dark default.
            }
            applyTheme(storedTheme);
        });

        menus.forEach(menu => {
            menu.addEventListener('toggle', () => {
                if (!menu.open) return;
                menus.forEach(candidate => {
                    if (candidate !== menu) candidate.removeAttribute('open');
                });
            });
            menu.addEventListener('click', event => {
                if (event.target.closest('.lc-wb-menu-command')) {
                    menu.removeAttribute('open');
                }
            });
        });

        root.addEventListener('pointerdown', event => {
            menus.forEach(menu => {
                if (!menu.contains(event.target)) menu.removeAttribute('open');
            });
        });

        root.addEventListener('keydown', event => {
            if (event.key !== 'Escape') return;
            menus.forEach(menu => menu.removeAttribute('open'));
        });

        const synchronizeSidebar = () => {
            const collapsed = root.dataset.sidebarState === 'collapsed';
            sidebarToggle.setAttribute('aria-label', collapsed
                ? 'Expand navigation'
                : 'Collapse navigation');
            sidebarToggle.dataset.tooltip = collapsed
                ? 'Expand navigation'
                : 'Collapse navigation';
            sidebarToggle.querySelector('svg').style.transform = collapsed
                ? 'rotate(180deg)'
                : 'none';
        };

        const synchronizeDock = () => {
            const collapsed = root.dataset.dockState === 'collapsed';
            dockToggle.setAttribute('aria-label', collapsed
                ? 'Restore output panel'
                : 'Collapse output panel');
            dockToggle.dataset.tooltip = collapsed
                ? 'Restore output panel'
                : 'Collapse output panel';
            dockToggle.querySelector('svg').style.transform = collapsed
                ? 'rotate(180deg)'
                : 'none';
        };

        sidebarToggle.addEventListener('click', () => {
            root.dataset.sidebarState = root.dataset.sidebarState === 'collapsed'
                ? 'expanded'
                : 'collapsed';
            synchronizeSidebar();
        });

        if (dockToggle) {
            dockToggle.addEventListener('click', () => {
                root.dataset.dockState = root.dataset.dockState === 'collapsed'
                    ? 'expanded'
                    : 'collapsed';
                synchronizeDock();
            });
        }

        if (splitter) {
            let dragging = false;
            splitter.addEventListener('pointerdown', event => {
                if (root.dataset.dockState === 'collapsed') return;
                dragging = true;
                splitter.setPointerCapture(event.pointerId);
                root.classList.add('is-splitting');
                event.preventDefault();
            });
            splitter.addEventListener('pointermove', event => {
                if (!dragging) return;
                const bounds = root.getBoundingClientRect();
                const height = Math.min(
                    bounds.height * 0.55,
                    Math.max(112, bounds.bottom - event.clientY - 26)
                );
                root.style.setProperty('--lc-wb-dock-height', `${height}px`);
            });
            const stop = event => {
                if (!dragging) return;
                dragging = false;
                root.classList.remove('is-splitting');
                if (splitter.hasPointerCapture(event.pointerId)) {
                    splitter.releasePointerCapture(event.pointerId);
                }
            };
            splitter.addEventListener('pointerup', stop);
            splitter.addEventListener('pointercancel', stop);
        }

        synchronizeSidebar();
        synchronizeDock();
    })();
    """
end

function Bonito.jsrender(session::Session, runtime::Runtime)
    ComponentXRay.set_policy!(session, runtime.xray)
    workbench = runtime.workbench
    dispatch = action -> handle!(runtime.application, runtime.state, action)
    scope = "$(workbench.namespace)-$(uuid4())"
    inspector = isnothing(workbench.inspector) ? nothing :
        render_inspector(session, workbench.inspector)
    output = isnothing(workbench.output) ? nothing :
        render_dock(session, workbench.output)
    status = isnothing(workbench.status) ? nothing : render_status(session, workbench.status)
    attributes = Dict{Symbol,Any}(
        Symbol("data-lc-workbench") => scope,
        Symbol("data-sidebar-state") => string(workbench.navigation.initial_state),
        Symbol("data-dock-state") => isnothing(workbench.output) ?
            "collapsed" : string(workbench.output.initial_state),
    )
    page_attributes = Dict{Symbol,Any}(
        Symbol("data-theme") => "dark",
        Symbol("data-resolved-theme") => "dark",
    )
    root = DOM.div(
        render_sidebar(
            session,
            workbench.identity,
            workbench.navigation,
            dispatch
        ),
        render_menubar(session, workbench.menu, dispatch, workbench.namespace),
        render_toolbar(session, workbench.toolbar, dispatch),
        render_viewstack(session, workbench.workspace),
        inspector,
        output,
        status;
        attributes...,
        id=scope,
        class=isnothing(inspector) ? "lc-wb-shell without-inspector" : "lc-wb-shell"
    )
    instrumented = ComponentXRay.instrument(session, root, workbench)
    xray = ComponentXRay.install(session, root)
    return Bonito.jsrender(
        session,
        DOM.div(
            DOM.style(BRAND_STYLES),
            DOM.style(CONTROL_CONTRACT),
            DOM.style(Toolkit.TOOLKIT_STYLES),
            DOM.style(WORKBENCH_STYLES),
            instrumented,
            DOM.script(intrinsic_script(root)),
            xray;
            page_attributes...,
            class="lc-wb-page"
        )
    )
end

end
