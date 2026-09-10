"""
    NavigationButton(label; href, icon=nothing, target="_top")

Navigate to an owned, root-relative page using the shared secondary-button
appearance. `icon` is optional rendered icon content; `label` remains the
accessible name when a sidebar contracts to an icon rail. `target` is `_top`
(leave an embedded workbench) or `_self` (stay in its browsing context).

Rendering installs no callback and starts no operation. Native link semantics
preserve keyboard navigation and opening in another tab without a live session.
Script URLs, network-relative URLs, backslashes and whitespace are rejected.
"""
struct NavigationButton{I}
    "Visible label and accessible name."
    label::String
    "Root-relative destination, optionally including a query or fragment."
    href::String
    "Rendered decorative icon, or nothing."
    icon::I
    "Navigation context: _top or _self."
    target::String

    function NavigationButton(label::AbstractString; href::AbstractString,
            icon=nothing, target::AbstractString="_top")
        isempty(strip(label)) && throw(ArgumentError("navigation requires a label"))
        startswith(href, "/") && !startswith(href, "//") &&
            !occursin(r"[\s\x00-\x1f\x7f\\]", href) ||
            throw(ArgumentError("navigation requires an owned root-relative destination"))
        target in ("_top", "_self") || throw(ArgumentError("unsupported navigation context"))
        return new{typeof(icon)}(string(label), string(href), icon, string(target))
    end
end

function Bonito.jsrender(session::Session, button::NavigationButton)
    glyph = isnothing(button.icon) ? nothing : DOM.span(button.icon;
        class="lc-navigation-icon", var"aria-hidden"="true")
    node = DOM.a(glyph, DOM.span(button.label; class="lc-navigation-label");
        class="lc-button lc-button-secondary lc-navigation-button",
        href=button.href, target=button.target, title=button.label,
        var"aria-label"=button.label)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, button))
end

function ComponentXRay.inspection(button::NavigationButton)
    return ComponentXRay.ComponentInspection(button; name="NavigationButton",
        source=toolkit_source(@__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:label, button.label),
            ComponentXRay.PropertyInspection(:href, button.href),
            ComponentXRay.PropertyInspection(:target, button.target)],
        css_scopes=[".lc-button", ".lc-navigation-button", ".lc-navigation-icon", ".lc-navigation-label"])
end
