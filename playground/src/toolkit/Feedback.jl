"""
    StatusIndicator(label; tone=:neutral, busy=false)

Display an observable status label using the shared semantic colours and bold
type. `tone` is `:neutral`, `:info`, `:success`, `:warning`, or `:danger`.
`busy` displays indeterminate activity; neither colour nor animation establishes
readiness. The owner supplies the evidence and may change the three observables.
"""
struct StatusIndicator
    "Visible, non-sensitive status text."
    label::Observable{String}
    "Semantic colour role."
    tone::Observable{Symbol}
    "Whether the owner has reported ongoing activity."
    busy::Observable{Bool}
end

function feedback_tone(tone)
    tone in (:neutral, :info, :success, :warning, :danger) ||
        throw(ArgumentError("unsupported status tone"))
    return string(tone)
end

function StatusIndicator(label::AbstractString; tone::Symbol=:neutral, busy::Bool=false)
    feedback_tone(tone)
    return StatusIndicator(Observable(string(label)), Observable(tone), Observable(busy))
end

function Bonito.jsrender(session::Session, status::StatusIndicator)
    node = DOM.span(status.label; class="lc-status-indicator lc-activity-status",
        role="status", var"aria-atomic"="true",
        var"data-tone"=feedback_tone(status.tone[]),
        var"data-busy"=string(status.busy[]))
    # Bonito's generic binding updates properties; CSS and ARIA consume attributes.
    tone = map(feedback_tone, session, status.tone)
    onjs(session, tone, js"value => $(node).setAttribute('data-tone', value)")
    onjs(session, status.busy, js"value => $(node).setAttribute('data-busy', String(value))")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, status))
end

"""
    ActionButton(label; busy_label="Working…", busy=false, disabled=false)

Emit explicit activations through `clicks`. The owner sets `busy` while its
action is pending and clears it on completion or failure. Busy buttons show the
shared activity glyph and label, expose `aria-busy`, and reject further clicks.
Construction and rendering perform no action and start no background task.
"""
struct ActionButton
    "Label while no action is pending."
    label::String
    "Label during the owner's action."
    busy_label::String
    "Owner-reported activity."
    busy::Observable{Bool}
    "Owner-reported disabled state."
    disabled::Observable{Bool}
    "Monotonic count of accepted activations."
    clicks::Observable{Int}
end

ActionButton(label::AbstractString; busy_label::AbstractString="Working…", busy::Bool=false,
    disabled::Bool=false) = ActionButton(string(label), string(busy_label),
        Observable(busy), Observable(disabled), Observable(0))

function Bonito.jsrender(session::Session, button::ActionButton)
    incoming = Observable(0)
    on(session, incoming) do _
        button.busy[] || button.disabled[] || (button.clicks[] += 1)
        return nothing
    end
    node = DOM.button(map(b -> b ? button.busy_label : button.label, session, button.busy);
        type="button", class="lc-button lc-button-secondary",
        disabled=map((b, d) -> b || d, session, button.busy, button.disabled),
        var"data-busy"=string(button.busy[]),
        var"aria-busy"=string(button.busy[]),
        onclick=js"event => $(incoming).notify($(incoming).value + 1)")
    onjs(session, button.busy, js"""value => {
        $(node).setAttribute('data-busy', String(value));
        $(node).setAttribute('aria-busy', String(value));
    }""")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, button))
end

function ComponentXRay.inspection(status::StatusIndicator)
    return ComponentXRay.ComponentInspection(status; name="StatusIndicator",
        source=toolkit_source(@__FILE__, @__LINE__),
        bindings=[ComponentXRay.BindingInspection(:label, status.label),
            ComponentXRay.BindingInspection(:tone, status.tone),
            ComponentXRay.BindingInspection(:busy, status.busy)],
        css_scopes=[".lc-status-indicator", ".lc-activity-status"])
end

function ComponentXRay.inspection(button::ActionButton)
    return ComponentXRay.ComponentInspection(button; name="ActionButton",
        source=toolkit_source(@__FILE__, @__LINE__),
        bindings=[ComponentXRay.BindingInspection(:busy, button.busy),
            ComponentXRay.BindingInspection(:disabled, button.disabled),
            ComponentXRay.BindingInspection(:clicks, button.clicks)],
        css_scopes=[".lc-button"])
end
