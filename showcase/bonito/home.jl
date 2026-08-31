module PlaygroundHome

using Bonito
using Markdown
using Observables: off
using ..PageAuthoring

const MAX_ACTIVITY_LINES = 80

function first_target(presentation)
    deck = first(presentation.decks)
    page = first(filter(candidate -> candidate.render, deck.pages))
    return "./presentations/$(presentation.slug)/$(deck.id)#$(page.id)"
end

function unique_resources(presentations)
    resources = Any[]
    for presentation in presentations, resource in presentation.resources

        any(candidate -> candidate.state === resource.state, resources) ||
            push!(resources, resource)
    end
    return resources
end

function unique_resource_states(presentations)
    return getproperty.(unique_resources(presentations), :state)
end

function activation_observable(session::Session, selector, presentations)
    states = unique_resource_states(presentations)
    return map(
        (presentation, snapshots...) -> preflight_summary(presentation.resources),
        session,
        selector.value,
        states...
    )
end

function activate!(presentation)
    for resource in presentation.resources
        resource.required || continue
        phase = resource.state[].phase
        phase === :hot && continue
        phase ∉ (:cold, :error) && continue
        try
            resource.activate()
        catch error
            set_preflight!(
                resource.state,
                :error,
                resource.state[].progress,
                "Activation failed: $(sprint(showerror, error))"
            )
        end
    end
    return preflight_summary(presentation.resources)
end

function activity_writer(terminal::TerminalOutput)
    lines = String[]
    activity_lock = ReentrantLock()
    function replace!(messages)
        lock(activity_lock) do
            empty!(lines)
            Base.append!(lines, string.(messages))
            length(lines) > MAX_ACTIVITY_LINES &&
                deleteat!(lines, 1:(length(lines) - MAX_ACTIVITY_LINES))
            terminal.content[] = isempty(lines) ? "" : "$(join(lines, '\n'))\n"
        end
        return nothing
    end
    function push_line!(message)
        lock(activity_lock) do
            push!(lines, string(message))
            length(lines) > MAX_ACTIVITY_LINES && popfirst!(lines)
            terminal.content[] = "$(join(lines, '\n'))\n"
        end
        return nothing
    end
    return (; replace!, push_line!)
end

function phase_marker(phase)
    phase === :hot ? "HOT" :
    phase === :error ? "ERROR" :
    phase === :cold ? "COLD" : uppercase(string(phase))
end

function activation_label(phase)
    phase === :hot ? "Activated" :
    phase === :preparing ? "Activating…" :
    phase === :error ? "Retry activation" : "Activate"
end

function presentation_report(presentation)
    messages = String[
    "LineCableModels playground preflight",
    "Presentation: $(presentation.title)",
    ""
]
    for deck in presentation.decks
        if isempty(deck.resources)
            push!(messages, "[LAZY] $(deck.title) · no preflight required")
            continue
        end
        for resource in deck.resources
            snapshot = resource.state[]
            push!(
                messages,
                "[$(phase_marker(snapshot.phase))] $(resource.title) · $(snapshot.label)"
            )
        end
    end
    summary = preflight_summary(presentation.resources)
    push!(messages, "", "[$(phase_marker(summary.phase))] $(summary.label)")
    return messages
end

function build(session::Session, presentations, assets)
    isempty(presentations) && throw(ArgumentError(
        "the playground home requires at least one discovered presentation"
    ))

    selector = Bonito.Dropdown(
        collect(presentations);
        option_to_string = presentation -> presentation.title,
        id = "playground-presentation-select",
        ariaLabel = "Available presentations"
    )
    activation = activation_observable(session, selector, presentations)
    phase = map(snapshot -> string(snapshot.phase), session, activation)
    ready = map(snapshot -> snapshot.phase === :hot, session, activation)
    initial_phase = activation[].phase
    activate_button = Bonito.Button(
        activation_label(initial_phase);
        id = "playground-presentation-activate",
        class = "lc-home-activate-button",
        style = nothing,
        disabled = initial_phase ∈ (:hot, :preparing),
        ariaDisabled = string(initial_phase ∈ (:hot, :preparing))
    )
    open_button = Bonito.Button(
        "Open";
        id = "playground-presentation-open",
        class = "lc-home-open-button",
        style = nothing,
        disabled = initial_phase !== :hot,
        ariaDisabled = string(initial_phase !== :hot)
    )
    action_buttons = DOM.div(
        activate_button,
        open_button;
        class = "lc-home-actions"
    )
    targets = first_target.(presentations)

    terminal_style = Styles(
        "width" => "100%",
        "height" => "auto",
        "min-height" => "100%",
        "padding" => "0.8rem 1rem",
        "overflow" => "visible",
        "background" => "#171b1b",
        "color" => "#dbdee0",
        "font-family" => "JuliaMono, monospace",
        "font-size" => "0.8rem",
        "line-height" => "1.45",
        "white-space" => "pre-wrap",
        "overflow-wrap" => "anywhere"
    )
    terminal = TerminalOutput(; style = terminal_style)
    activity = activity_writer(terminal)
    activity.replace!(presentation_report(selector.value[]))

    selection_observer = on(selector.value) do presentation
        activity.replace!(presentation_report(presentation))
        return nothing
    end
    resource_observers = map(unique_resources(presentations)) do resource
        on(resource.state) do snapshot
            selected = selector.value[]
            any(candidate -> candidate.state === resource.state, selected.resources) ||
                return nothing
            percent = round(Int, 100snapshot.progress)
            activity.push_line!(
                "[$(phase_marker(snapshot.phase))] $(resource.title) · $(snapshot.label) · $percent%"
            )
            return nothing
        end
    end
    activate_observer = on(activate_button.value) do clicked
        clicked || return nothing
        presentation = selector.value[]
        activity.push_line!("[ACTIVATE] Activating $(presentation.title)")
        activate!(presentation)
        return nothing
    end

    Bonito.onjs(
        session,
        open_button.value,
        js"""
clicked => {
    if (!clicked || !$(ready).value) {
        return;
    }
    const targets = $(targets);
    const index = $(selector.option_index).value - 1;
    const target = targets[index];
    if (target) {
        window.location.assign(target);
    }
}
"""
    )
    Bonito.onjs(
        session,
        phase,
        js"""
phase => {
    const actions = $(action_buttons);
    const activateButton = actions.querySelector("#playground-presentation-activate");
    const openButton = actions.querySelector("#playground-presentation-open");
    const hot = phase === "hot";
    const preparing = phase === "preparing";
    activateButton.disabled = hot || preparing;
    activateButton.setAttribute("aria-disabled", String(hot || preparing));
    activateButton.textContent = hot ? "Activated" :
        preparing ? "Activating…" :
        phase === "error" ? "Retry activation" : "Activate";
    openButton.disabled = !hot;
    openButton.setAttribute("aria-disabled", String(!hot));
}
"""
    )

    introduction = webpart(
        prose(md"""
        ## One publisher, many purposes

        The playground serves independent collections of stateful Bonito decks. A collection can be a live manual, a focused tutorial, an application study, or a presentation assembled for one audience.

        Each presentation owns its modules and runtime state. Activate its declared process-wide resources here; deck sessions, widgets, figures, and REPL workers remain lazy until opened.
        """);
        kind = :panel,
        class = "lc-home-introduction"
    )
    chooser = webpart(
        DOM.img(;
            src = assets.logo,
            alt = "LineCableModels logo",
            class = "lc-home-logo"
        ),
        control(
            "Presentation",
            selector;
            value = "$(length(presentations)) available",
            id = "playground-presentation-control"
        ),
        preparation_status(
            session;
            state = activation,
            title = "Presentation preflight"
        ),
        action_buttons;
        kind = :controls,
        class = "lc-home-chooser"
    )
    terminal_host = DOM.div(
        terminal;
        id = "playground-preflight-terminal",
        class = "lc-home-terminal"
    )
    activity_panel = webpart(
        terminal_host;
        kind = :panel,
        title = "Preflight activity",
        meta = "explicit publisher events · bounded",
        id = "playground-preflight-console",
        class = "lc-home-activity",
        body_class = "lc-webpart-body-flush"
    )
    canvas = layout_split_bottom(
        introduction,
        chooser,
        activity_panel;
        ratio = :equal,
        height = "100%",
        class = "lc-home-layout"
    )

    Bonito.onjs(
        session,
        terminal.content,
        js"""
() => {
    const terminalHost = $(terminal_host);
    requestAnimationFrame(() => {
        terminalHost.scrollTop = terminalHost.scrollHeight;
    });
}
"""
    )

    on(session.on_close) do _
        off(selection_observer)
        off(activate_observer)
        foreach(off, resource_observers)
        return nothing
    end

    return (;
        body = slide(
            "LineCableModels playground",
            canvas;
            lede = md"Choose, activate, and open what you want to publish, teach, or explore."
        ),
        selector,
        activate_button,
        open_button,
        activation,
        terminal
    )
end

end
