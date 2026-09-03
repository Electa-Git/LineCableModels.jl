const RIBBON_STYLES_PATH = normpath(joinpath(@__DIR__, "..", "assets", "ribbon.css"))
include_dependency(RIBBON_STYLES_PATH)
const RIBBON_STYLES = read(RIBBON_STYLES_PATH, String)

struct RibbonGroup{T<:Tuple}
    label::String
    items::T
    size::Symbol
end

function RibbonGroup(
        label::AbstractString,
        items::AbstractToolbarItem...;
        size::Symbol=:medium
    )
    size in TOOLBAR_SIZES || throw(ArgumentError(
        "ribbon group size must be :small, :medium, or :large"
    ))
    isempty(items) && throw(ArgumentError("ribbon group requires an item"))
    actions = [item.action for item in items if !(item isa ToolbarSeparator)]
    allunique(actions) || throw(ArgumentError(
        "ribbon action names must be unique within group `$label`"
    ))
    return RibbonGroup(string(label), items, size)
end

struct RibbonTab{G<:Tuple}
    id::Symbol
    label::String
    groups::G
    disabled::Bool
end

function RibbonTab(
        id::Symbol,
        label::AbstractString,
        groups::RibbonGroup...;
        disabled::Bool=false
    )
    valid_toolbar_name(id, "ribbon tab")
    isempty(groups) && throw(ArgumentError("ribbon tab requires a group"))
    return RibbonTab(id, string(label), groups, disabled)
end

struct Ribbon
    dom::Any
    events::Observable{Union{Nothing,ToolbarEvent}}
    behavior::Any
end

function ribbon_group(group::RibbonGroup, events, binding)
    controls = [toolbar_item(item, events, binding) for item in group.items]
    return DOM.section(
        DOM.div(controls...; class="lc-ribbon-group-controls"),
        DOM.h3(group.label; class="lc-ribbon-group-label");
        class="lc-ribbon-group lc-ribbon-group-$(group.size)",
        var"data-ribbon-group"=""
    )
end

function Ribbon(
        tabs::RibbonTab...;
        binding::ToolbarBinding,
        active::Union{Nothing,Symbol}=nothing,
        collapsed::Bool=false,
        quick_access::Tuple=(),
        label::AbstractString="Application ribbon"
    )
    isempty(tabs) && throw(ArgumentError("ribbon requires a tab"))
    tab_ids = [tab.id for tab in tabs]
    allunique(tab_ids) || throw(ArgumentError("ribbon tab ids must be unique"))
    resolved_active = isnothing(active) ? first(tab_ids) : active
    resolved_active in tab_ids || throw(ArgumentError(
        "active ribbon tab is not declared: $resolved_active"
    ))
    active_tab = only(tab for tab in tabs if tab.id == resolved_active)
    active_tab.disabled && throw(ArgumentError("active ribbon tab cannot be disabled"))
    all(item isa ToolbarButton for item in quick_access) || throw(ArgumentError(
        "ribbon quick access accepts ToolbarButton items only"
    ))
    allunique(item.action for item in quick_access) || throw(ArgumentError(
        "ribbon quick-access action names must be unique"
    ))

    events = Observable{Union{Nothing,ToolbarEvent}}(nothing)
    ribbon_id = "lc-ribbon-$(uuid4())"
    tab_buttons = map(tabs) do tab
        panel_id = "$ribbon_id-panel-$(tab.id)"
        attributes = Dict{Symbol,Any}(
            Symbol("aria-controls") => panel_id,
            Symbol("aria-selected") => string(tab.id == resolved_active),
            Symbol("data-ribbon-tab") => string(tab.id),
            :id => "$ribbon_id-tab-$(tab.id)",
        )
        return DOM.button(
            tab.label;
            attributes...,
            type="button",
            role="tab",
            tabindex=tab.id == resolved_active ? "0" : "-1",
            disabled=tab.disabled,
            class="lc-ribbon-tab"
        )
    end
    tab_panels = map(tabs) do tab
        strip = DOM.div(
            (ribbon_group(group, events, binding) for group in tab.groups)...;
            class="lc-ribbon-group-strip",
            var"data-ribbon-group-strip"=""
        )
        overflow_menu = DOM.div(; class="lc-ribbon-overflow-menu")
        overflow = DOM.details(
            DOM.summary(
                DOM.span("More"),
                toolbar_icon(:next);
                var"aria-label"="More ribbon groups"
            ),
            overflow_menu;
            class="lc-ribbon-overflow",
            var"data-ribbon-overflow"="",
            hidden=true
        )
        attributes = Dict{Symbol,Any}(
            Symbol("aria-labelledby") => "$ribbon_id-tab-$(tab.id)",
            Symbol("data-ribbon-panel") => string(tab.id),
            :id => "$ribbon_id-panel-$(tab.id)",
        )
        return DOM.div(
            strip,
            overflow;
            attributes...,
            role="tabpanel",
            hidden=tab.id != resolved_active,
            class="lc-ribbon-panel"
        )
    end
    quick = DOM.div(
        (toolbar_button(item, events, binding) for item in quick_access)...;
        class="lc-ribbon-quick-access",
        var"aria-label"="Quick access",
        role="toolbar"
    )
    collapse = DOM.button(
        toolbar_icon(:collapse),
        DOM.span("Collapse ribbon"; class="lc-visually-hidden");
        type="button",
        class="lc-ribbon-collapse",
        var"data-ribbon-collapse"="",
        var"aria-expanded"=string(!collapsed),
        title=collapsed ? "Expand ribbon" : "Collapse ribbon"
    )
    root_attributes = Dict{Symbol,Any}(
        Symbol("aria-label") => string(label),
        Symbol("data-ribbon-root") => "",
        Symbol("data-collapsed") => string(collapsed),
        :id => ribbon_id,
    )
    root = DOM.section(
        DOM.style(RIBBON_STYLES),
        DOM.header(
            quick,
            DOM.div(
                tab_buttons...;
                class="lc-ribbon-tabs",
                role="tablist",
                var"aria-label"=string(label)
            ),
            collapse;
            class="lc-ribbon-tab-row"
        ),
        DOM.div(tab_panels...; class="lc-ribbon-panels");
        root_attributes...,
        class="lc-ribbon",
        role="region"
    )

    behavior = js"""
        (element) => {
            const tabs = Array.from(element.querySelectorAll('[data-ribbon-tab]'));
            const panels = Array.from(element.querySelectorAll('[data-ribbon-panel]'));
            const collapse = element.querySelector('[data-ribbon-collapse]');

            const closeMenus = (except = null) => {
                element.querySelectorAll('details[open]').forEach(menu => {
                    if (menu !== except) menu.removeAttribute('open');
                });
            };

            const rebalance = (panel) => {
                if (!panel || panel.hidden || element.dataset.collapsed === 'true') return;
                const strip = panel.querySelector('[data-ribbon-group-strip]');
                const overflow = panel.querySelector('[data-ribbon-overflow]');
                const menu = overflow.querySelector('.lc-ribbon-overflow-menu');

                Array.from(menu.children).forEach(group => strip.append(group));
                overflow.hidden = false;
                strip.style.maxWidth = `calc(100% - ${overflow.offsetWidth + 8}px)`;
                while (strip.scrollWidth > strip.clientWidth + 1 && strip.children.length > 1) {
                    menu.prepend(strip.lastElementChild);
                }
                const hasOverflow = menu.children.length > 0;
                overflow.hidden = !hasOverflow;
                if (!hasOverflow) strip.style.maxWidth = '100%';
                if (!hasOverflow) overflow.removeAttribute('open');
            };

            const rebalanceActive = () => {
                requestAnimationFrame(() => {
                    panels.filter(panel => !panel.hidden).forEach(rebalance);
                });
            };

            const activate = (id, focus = false) => {
                const target = tabs.find(tab => tab.dataset.ribbonTab === id && !tab.disabled);
                if (!target) return;
                tabs.forEach(tab => {
                    const selected = tab === target;
                    tab.setAttribute('aria-selected', String(selected));
                    tab.tabIndex = selected ? 0 : -1;
                });
                panels.forEach(panel => {
                    panel.hidden = panel.dataset.ribbonPanel !== id;
                });
                closeMenus();
                if (focus) target.focus();
                rebalanceActive();
            };

            tabs.forEach(tab => {
                tab.addEventListener('click', () => activate(tab.dataset.ribbonTab));
                tab.addEventListener('keydown', event => {
                    if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
                    event.preventDefault();
                    const enabled = tabs.filter(candidate => !candidate.disabled);
                    const index = enabled.indexOf(tab);
                    let target = index;
                    if (event.key === 'ArrowLeft') target = (index - 1 + enabled.length) % enabled.length;
                    if (event.key === 'ArrowRight') target = (index + 1) % enabled.length;
                    if (event.key === 'Home') target = 0;
                    if (event.key === 'End') target = enabled.length - 1;
                    activate(enabled[target].dataset.ribbonTab, true);
                });
            });

            collapse.addEventListener('click', () => {
                const collapsed = element.dataset.collapsed !== 'true';
                element.dataset.collapsed = String(collapsed);
                collapse.setAttribute('aria-expanded', String(!collapsed));
                collapse.title = collapsed ? 'Expand ribbon' : 'Collapse ribbon';
                closeMenus();
                if (!collapsed) rebalanceActive();
            });

            element.addEventListener('pointerdown', event => {
                const current = event.target.closest('details');
                closeMenus(current);
            });
            element.addEventListener('click', event => {
                if (event.target.closest('[data-toolbar-action]')) closeMenus();
            });
            const outsideClick = event => {
                if (!element.contains(event.target)) closeMenus();
            };
            document.addEventListener('pointerdown', outsideClick);

            const observer = new ResizeObserver(rebalanceActive);
            observer.observe(element);
            rebalanceActive();
            window.addEventListener('pagehide', () => {
                observer.disconnect();
                document.removeEventListener('pointerdown', outsideClick);
            }, {once: true});
        }
    """
    return Ribbon(root, events, behavior)
end

function Bonito.jsrender(session::Session, ribbon::Ribbon)
    Bonito.onload(session, ribbon.dom, ribbon.behavior)
    return Bonito.jsrender(session, ribbon.dom)
end
