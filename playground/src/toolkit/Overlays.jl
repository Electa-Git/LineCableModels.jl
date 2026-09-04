const OVERLAY_KINDS = (:info, :success, :warning, :danger)

function overlay_kind(kind)
    resolved = Symbol(kind)
    resolved in OVERLAY_KINDS || throw(ArgumentError(
        "overlay kind must be :info, :success, :warning, or :danger"))
    return resolved
end

"""Describe one typed action exposed by a modal [`Dialog`](@ref)."""
struct DialogAction{C,V}
    "Stable action identity."
    id::Symbol
    "Visible button label."
    label::String
    "Optional Julia callback invoked on activation."
    callback::C
    "Value carried by the emitted event."
    value::V
    "Visual action treatment."
    kind::Symbol
    "Whether activation closes the dialog."
    closes::Bool
    "Whether activation is disabled."
    disabled::Bool
end

function DialogAction(id, label::AbstractString; callback=nothing, value=id,
        kind::Symbol=:secondary, closes::Bool=true, disabled::Bool=false)
    kind in (:primary, :secondary, :danger) || throw(ArgumentError(
        "dialog action kind must be :primary, :secondary, or :danger"))
    return DialogAction(form_name(id, "dialog action"), string(label), callback,
        value, kind, closes, disabled)
end

"""Record a dialog action without coupling the component to page logic."""
struct DialogEvent{T}
    "Dialog identity."
    dialog::Symbol
    "Activated action identity."
    action::Symbol
    "Action payload."
    value::T
end

"""
    Dialog(name, title, content; actions=(), open=false, dismissible=true)

Create a reusable native modal dialog with focus trapping and restoration.
Only the `open` observable and typed action events cross the browser boundary.
"""
struct Dialog{C,A}
    "Stable dialog identity."
    name::Symbol
    "Visible heading."
    title::String
    "Renderable dialog body."
    content::C
    "Ordered action tuple."
    actions::A
    "Whether the modal is requested open."
    open::Observable{Bool}
    "Latest dialog event."
    event::Observable{Any}
    "Whether Escape and backdrop clicks dismiss the dialog."
    dismissible::Bool
    "Accessible dialog description."
    description::Union{Nothing,String}
end

function Dialog(name, title::AbstractString, content; actions=(), open::Bool=false,
        dismissible::Bool=true, description::Union{Nothing,AbstractString}=nothing)
    resolved_actions = Tuple(actions)
    all(action -> action isa DialogAction, resolved_actions) || throw(ArgumentError(
        "every dialog action must be a DialogAction"))
    allunique(action.id for action in resolved_actions) || throw(ArgumentError(
        "dialog action identities must be unique"))
    return Dialog(form_name(name, "dialog name"), string(title), content,
        resolved_actions, Observable(open), Observable{Any}(nothing), dismissible,
        isnothing(description) ? nothing : string(description))
end

"""Request that `dialog` open in its current browser session."""
open_dialog!(dialog::Dialog) = (dialog.open[] = true; dialog)

"""Request that `dialog` close without invoking an action callback."""
function close_dialog!(dialog::Dialog; action::Symbol=:dismiss, value=nothing)
    dialog.event[] = DialogEvent(dialog.name, action, value)
    dialog.open[] = false
    return dialog
end

function activate_dialog!(dialog::Dialog, id::Symbol)
    index = findfirst(action -> action.id == id, dialog.actions)
    isnothing(index) && throw(KeyError(id))
    action = dialog.actions[index]
    action.disabled && return dialog
    if !isnothing(action.callback)
        if applicable(action.callback, action.value)
            action.callback(action.value)
        elseif applicable(action.callback)
            action.callback()
        else
            throw(ArgumentError("dialog callback accepts neither zero arguments nor its action value"))
        end
    end
    dialog.event[] = DialogEvent(dialog.name, action.id, action.value)
    action.closes && (dialog.open[] = false)
    return dialog
end

function dialog_close_icon()
    return SVG.svg(SVG.path(; d="M6 6l12 12M18 6 6 18"); viewBox="0 0 24 24",
        fill="none", stroke="currentColor", var"stroke-width"="1.8",
        var"stroke-linecap"="round", var"aria-hidden"="true")
end

function Bonito.jsrender(session::Session, dialog::Dialog)
    action_wire = Observable("")
    dismiss_wire = Observable(0)
    on(session, action_wire) do id
        isempty(id) || activate_dialog!(dialog, Symbol(id))
        return nothing
    end
    on(session, dismiss_wire) do _
        dialog.dismissible && close_dialog!(dialog)
        return nothing
    end
    actions = map(dialog.actions) do action
        DOM.button(action.label; type="button",
            class="lc-button lc-dialog-action is-$(action.kind)",
            disabled=action.disabled,
            onclick=js"event => $(action_wire).notify($(string(action.id)))")
    end
    close_button = dialog.dismissible ? DOM.button(dialog_close_icon(); type="button",
        class="lc-dialog-close", title="Close", var"aria-label"="Close",
        onclick=js"event => $(dismiss_wire).notify($(dismiss_wire).value + 1)") : nothing
    description_id = "$(dialog.name)-description"
    node = DOM.dialog(
        DOM.div(DOM.h2(dialog.title; id="$(dialog.name)-title"), close_button;
            class="lc-dialog-header"),
        isnothing(dialog.description) ? nothing :
            DOM.p(dialog.description; id=description_id, class="lc-dialog-description"),
        DOM.div(dialog.content; class="lc-dialog-body"),
        isempty(actions) ? nothing : DOM.div(actions...; class="lc-dialog-actions");
        id=string(dialog.name), class="lc-dialog", var"aria-labelledby"="$(dialog.name)-title",
        var"aria-describedby"=isnothing(dialog.description) ? nothing : description_id,
        oncancel=dialog.dismissible ? js"event => { event.preventDefault(); $(dismiss_wire).notify($(dismiss_wire).value + 1); }" :
            js"event => event.preventDefault()")
    behavior = js"""
    (() => {
        const dialog = $(node);
        const state = $(dialog.open);
        let restoreFocus = null;
        const synchronize = requested => {
            if (requested && !dialog.open) {
                restoreFocus = document.activeElement;
                dialog.showModal();
                const first = dialog.querySelector('button:not(:disabled), input:not(:disabled), select:not(:disabled), textarea:not(:disabled)');
                if (first) first.focus({preventScroll: true});
            } else if (!requested && dialog.open) {
                dialog.close();
                if (restoreFocus && restoreFocus.isConnected) restoreFocus.focus({preventScroll: true});
            }
        };
        state.on(synchronize);
        synchronize(state.value);
        dialog.addEventListener('click', event => {
            if (!$(dialog.dismissible) || event.target !== dialog) return;
            const bounds = dialog.getBoundingClientRect();
            const inside = event.clientX >= bounds.left && event.clientX <= bounds.right
                && event.clientY >= bounds.top && event.clientY <= bounds.bottom;
            if (!inside) $(dismiss_wire).notify($(dismiss_wire).value + 1);
        });
    })();
    """
    rendered = ComponentXRay.instrument(session,
        DOM.div(node, DOM.script(behavior); style="display: contents;"), dialog)
    return Bonito.jsrender(session, rendered)
end

"""Construct an informational modal with one acknowledgement action."""
function MessageDialog(name, title, message; acknowledge_label="Close", kwargs...)
    action = DialogAction(:acknowledge, acknowledge_label; kind=:primary)
    return Dialog(name, title, DOM.p(message; class="lc-dialog-copy");
        actions=(action,), kwargs...)
end

"""Construct a confirmation modal with cancel and confirm actions."""
function ConfirmDialog(name, title, message; onconfirm=nothing,
        confirm_label="Confirm", cancel_label="Cancel", danger::Bool=false, kwargs...)
    actions = (
        DialogAction(:cancel, cancel_label),
        DialogAction(:confirm, confirm_label; callback=onconfirm,
            value=true, kind=danger ? :danger : :primary),
    )
    return Dialog(name, title, DOM.p(message; class="lc-dialog-copy"); actions, kwargs...)
end

"""Construct a modal whose body is an existing reusable [`Form`](@ref)."""
function FormDialog(name, title, form::Form; dismissible::Bool=true, description=nothing)
    dialog = Dialog(name, title, form; dismissible, description)
    on(form.submitted) do payload
        # The form remains the payload owner. Dialog diagnostics report only
        # the lifecycle event and never duplicate submitted secret fields.
        isnothing(payload) || close_dialog!(dialog; action=:submit)
        return nothing
    end
    on(form.cancelled) do count
        count > 0 && close_dialog!(dialog; action=:cancel)
        return nothing
    end
    return dialog
end

"""Create an inline semantic notice with optional dismissal."""
struct InlineNotice{C}
    "Semantic notice kind."
    kind::Symbol
    "Visible heading."
    title::String
    "Renderable notice body."
    content::C
    "Whether the notice is currently visible."
    visible::Observable{Bool}
    "Whether a dismiss action is offered."
    dismissible::Bool
end

function InlineNotice(kind, title::AbstractString, content; dismissible::Bool=false,
        visible::Bool=true)
    return InlineNotice(overlay_kind(kind), string(title), content,
        Observable(visible), dismissible)
end

"""Hide a dismissible inline notice."""
dismiss_notice!(notice::InlineNotice) = (notice.visible[] = false; notice)

function Bonito.jsrender(session::Session, notice::InlineNotice)
    dismiss = Observable(0)
    on(session, dismiss) do _
        notice.dismissible && dismiss_notice!(notice)
        return nothing
    end
    class = map(session, notice.visible) do visible
        visible ? "lc-notice is-$(notice.kind)" : "lc-notice is-$(notice.kind) is-hidden"
    end
    button = notice.dismissible ? DOM.button(dialog_close_icon(); type="button",
        class="lc-notice-close", title="Dismiss", var"aria-label"="Dismiss notice",
        onclick=js"event => $(dismiss).notify($(dismiss).value + 1)") : nothing
    node = DOM.aside(DOM.div(DOM.strong(notice.title), notice.content;
        class="lc-notice-copy"), button; class, role=notice.kind == :danger ? "alert" : "status")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, notice))
end

"""Represent one immutable notification managed by a [`ToastCenter`](@ref)."""
struct ToastEntry
    "Stable notification identity."
    id::UUID
    "Semantic notification kind."
    kind::Symbol
    "Visible heading."
    title::String
    "Visible message."
    message::String
    "Automatic dismissal delay in milliseconds, or zero to persist."
    timeout_ms::Int
end

function ToastEntry(kind, title::AbstractString, message::AbstractString;
        id::UUID=uuid4(), timeout_ms::Integer=5000)
    timeout_ms >= 0 || throw(ArgumentError("toast timeout must be nonnegative"))
    return ToastEntry(id, overlay_kind(kind), string(title), string(message), Int(timeout_ms))
end

"""Manage a bounded, session-local stack of transient notifications."""
struct ToastCenter
    "Ordered notification entries."
    entries::Observable{Vector{ToastEntry}}
    "Maximum retained notifications."
    capacity::Int
    "Accessible region label."
    label::String
end

function ToastCenter(; capacity::Integer=4, label::AbstractString="Notifications")
    capacity > 0 || throw(ArgumentError("toast capacity must be positive"))
    return ToastCenter(Observable(ToastEntry[]), Int(capacity), string(label))
end

"""Append a toast and evict the oldest entry when capacity is exceeded."""
function push_toast!(center::ToastCenter, entry::ToastEntry)
    entries = [center.entries[]; entry]
    center.entries[] = entries[max(1, length(entries) - center.capacity + 1):end]
    return entry
end

push_toast!(center::ToastCenter, kind, title, message; kwargs...) =
    push_toast!(center, ToastEntry(kind, title, message; kwargs...))

"""Dismiss a toast by stable identity."""
function dismiss_toast!(center::ToastCenter, id::UUID)
    center.entries[] = filter(entry -> entry.id != id, center.entries[])
    return center
end

"""Remove all current notifications."""
clear_toasts!(center::ToastCenter) = (center.entries[] = ToastEntry[]; center)

struct RenderedToast
    center::ToastCenter
    entry::ToastEntry
end

function Bonito.jsrender(session::Session, rendered::RenderedToast)
    entry = rendered.entry
    dismiss = Observable(0)
    on(session, dismiss) do _
        dismiss_toast!(rendered.center, entry.id)
        return nothing
    end
    node = DOM.div(
        DOM.div(DOM.strong(entry.title), DOM.span(entry.message); class="lc-toast-copy"),
        DOM.button(dialog_close_icon(); type="button", class="lc-toast-close",
            title="Dismiss", var"aria-label"="Dismiss notification",
            onclick=js"event => $(dismiss).notify($(dismiss).value + 1)");
        class="lc-toast is-$(entry.kind)", role=entry.kind == :danger ? "alert" : "status")
    timer = entry.timeout_ms == 0 ? nothing : DOM.script(js"""
    (() => {
        const timer = window.setTimeout(
            () => $(dismiss).notify($(dismiss).value + 1),
            $(entry.timeout_ms)
        );
        window.addEventListener('pagehide', () => window.clearTimeout(timer), {once: true});
    })();
    """)
    return Bonito.jsrender(session, DOM.div(node, timer; style="display: contents;"))
end

function Bonito.jsrender(session::Session, center::ToastCenter)
    rendered = Observable([RenderedToast(center, entry) for entry in center.entries[]])
    on(session, center.entries) do entries
        rendered[] = [RenderedToast(center, entry) for entry in entries]
        return nothing
    end
    node = DOM.div(KeyedList(rendered; key=item -> item.entry.id);
        class="lc-toast-center", role="region", var"aria-label"=center.label)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, center))
end

function ComponentXRay.inspection(dialog::Dialog)
    return ComponentXRay.ComponentInspection(dialog; name="Dialog",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:name, dialog.name),
            ComponentXRay.PropertyInspection(:title, dialog.title),
            ComponentXRay.PropertyInspection(:actions, length(dialog.actions)),
            ComponentXRay.PropertyInspection(:dismissible, dialog.dismissible)],
        bindings=[ComponentXRay.BindingInspection(:open, dialog.open),
            ComponentXRay.BindingInspection(:event, dialog.event)],
        css_scopes=[".lc-dialog", ".lc-dialog-header", ".lc-dialog-actions"])
end

function ComponentXRay.inspection(notice::InlineNotice)
    return ComponentXRay.ComponentInspection(notice; name="InlineNotice",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:kind, notice.kind),
            ComponentXRay.PropertyInspection(:title, notice.title),
            ComponentXRay.PropertyInspection(:dismissible, notice.dismissible)],
        bindings=[ComponentXRay.BindingInspection(:visible, notice.visible)],
        css_scopes=[".lc-notice", ".lc-notice-copy"])
end

function ComponentXRay.inspection(center::ToastCenter)
    return ComponentXRay.ComponentInspection(center; name="ToastCenter",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:capacity, center.capacity),
            ComponentXRay.PropertyInspection(:label, center.label)],
        bindings=[ComponentXRay.BindingInspection(:entries, center.entries)],
        css_scopes=[".lc-toast-center", ".lc-toast"])
end
