const WIDGET_THEME = raw"""
:root {
  color-scheme: dark;
  --lc-widget-bg: #1f2424;
  --lc-widget-panel: #252b2b;
  --lc-widget-panel-soft: #2b3232;
  --lc-widget-text: #dbdee0;
  --lc-widget-heading: #f2f2f2;
  --lc-widget-muted: #a8b0b0;
  --lc-widget-border: #5e6d6f;
  --lc-widget-link: #4eb5de;
  --lc-widget-focus: #1abc9c;
  --lc-kicker-heading-gap: 0.5rem;
  --bonito-widget-bg: #252b2b;
  --bonito-widget-fg: #dbdee0;
  --bonito-widget-border: #5e6d6f;
  --bonito-widget-hover-bg: #323a3a;
  --bonito-widget-muted-bg: #394142;
  --bonito-widget-accent: #1abc9c;
}

html,
body {
  width: 100%;
  min-width: 0;
  min-height: 100%;
  margin: 0;
  background: var(--lc-widget-bg);
  caret-color: transparent;
  cursor: default;
}

body,
button,
input,
select,
textarea {
  font-family: Lato, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}

input:not([type]):not([readonly]):not([disabled]),
input[type="text"]:not([readonly]):not([disabled]),
input[type="textfield"]:not([readonly]):not([disabled]),
input[type="search"]:not([readonly]):not([disabled]),
input[type="email"]:not([readonly]):not([disabled]),
input[type="url"]:not([readonly]):not([disabled]),
input[type="tel"]:not([readonly]):not([disabled]),
input[type="password"]:not([readonly]):not([disabled]),
input[type="number"]:not([readonly]):not([disabled]),
textarea:not([readonly]):not([disabled]),
[contenteditable="true"] {
  caret-color: auto;
  cursor: text;
  user-select: text;
  -webkit-user-select: text;
}

a[href],
button:not([disabled]),
select:not([disabled]),
summary,
label[for],
[role="button"]:not([aria-disabled="true"]),
input[type="button"]:not([disabled]),
input[type="submit"]:not([disabled]),
input[type="reset"]:not([disabled]),
input[type="checkbox"]:not([disabled]),
input[type="radio"]:not([disabled]),
input[type="range"]:not([disabled]) {
  cursor: pointer;
}

*,
*::before,
*::after {
  box-sizing: border-box;
}

.lc-widget-shell {
  display: grid;
  width: 100%;
  min-width: 0;
  min-height: 100vh;
  padding: clamp(1rem, 3vw, 1.6rem);
  grid-template-rows: auto minmax(0, 1fr);
  gap: 1rem;
  color: var(--lc-widget-text);
  background: var(--lc-widget-panel);
}

.lc-widget-header {
  display: flex;
  min-width: 0;
  align-items: end;
  justify-content: space-between;
  gap: 1rem;
}

.lc-widget-kicker {
  color: var(--lc-widget-focus);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.68rem;
  font-weight: 700;
  letter-spacing: 0.08em;
  line-height: 1;
}

.lc-widget-header h1 {
  margin: var(--lc-kicker-heading-gap) 0 0;
  color: var(--lc-widget-heading);
  font-size: clamp(1.15rem, 3vw, 1.45rem);
  line-height: 1.15;
}

.lc-widget-session {
  color: var(--lc-widget-muted);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.68rem;
  white-space: nowrap;
}

.lc-widget-content,
.lc-control-stack {
  display: grid;
  min-width: 0;
  align-content: start;
  gap: 0.9rem;
}

.lc-control-row {
  display: grid;
  min-width: 0;
  grid-template-columns: minmax(8rem, 0.75fr) minmax(0, 1.25fr);
  align-items: center;
  gap: 1rem;
}

.lc-control-copy {
  min-width: 0;
}

.lc-control-label {
  display: block;
  margin: 0 0 0.25rem;
  color: var(--lc-widget-heading);
  font-size: 0.88rem;
  font-weight: 700;
}

.lc-control-hint {
  margin: 0;
  color: var(--lc-widget-muted);
  font-size: 0.76rem;
  line-height: 1.4;
}

.lc-widget-output {
  display: flex;
  min-width: 0;
  min-height: 3rem;
  padding: 0.65rem 0.8rem;
  align-items: center;
  justify-content: space-between;
  gap: 1rem;
  background: #1d2222;
  border: 1px solid #465355;
  border-radius: 3px;
}

.lc-widget-output span:first-child {
  color: var(--lc-widget-muted);
  font-size: 0.78rem;
}

.lc-widget-value {
  overflow: hidden;
  color: var(--lc-widget-link);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.88rem;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.lc-widget-actions {
  display: flex;
  flex-wrap: wrap;
  gap: 0.6rem;
}

.lc-widget-actions button,
.lc-native-field,
.lc-native-number,
.lc-native-select {
  width: 100%;
  min-width: 0 !important;
  min-height: 2.55rem;
  margin: 0 !important;
  padding: 0.55rem 0.8rem !important;
  color: var(--lc-widget-text) !important;
  background: #202626 !important;
  border: 1px solid var(--lc-widget-border) !important;
  border-radius: 3px !important;
  box-shadow: none !important;
}

.lc-native-number {
  color-scheme: dark;
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-variant-numeric: tabular-nums;
}

.lc-native-range {
  width: 100%;
  min-width: 0;
  height: 1.5rem;
  margin: 0;
  accent-color: var(--lc-widget-focus);
  cursor: pointer;
}

.lc-widget-actions button {
  width: auto;
  flex: 1 1 8rem;
}

.lc-widget-actions button:first-child {
  color: #102523 !important;
  background: var(--lc-widget-focus) !important;
  border-color: var(--lc-widget-focus) !important;
}

.lc-widget-actions button:hover,
.lc-widget-actions button:focus-visible,
.lc-native-field:focus,
.lc-native-select:focus {
  outline: 2px solid rgb(78 181 222 / 45%) !important;
  outline-offset: 1px;
}

.lc-toolbar-showcase {
  display: grid;
  min-width: 0;
  grid-template-columns: minmax(0, 1fr) minmax(10rem, 0.42fr);
  gap: 0.8rem;
}

.lc-toolbar-example {
  display: grid;
  min-width: 0;
  padding: 0.75rem;
  align-content: start;
  justify-items: start;
  gap: 0.55rem;
  background: #1d2222;
  border: 1px solid #465355;
  border-radius: 3px;
}

.lc-toolbar-example-wide {
  grid-column: 1 / -1;
}

.lc-toolbar-example-label {
  color: var(--lc-widget-muted);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.68rem;
  letter-spacing: 0.05em;
  line-height: 1;
}

.lc-toolbar {
  --lc-toolbar-control: 2.45rem;
  --lc-toolbar-icon: 1.05rem;
  --lc-toolbar-pad-inline: 0.65rem;
  display: flex;
  width: max-content;
  max-width: 100%;
  padding: 0.28rem;
  align-items: stretch;
  gap: 0.25rem;
  color: var(--lc-widget-text);
  background: #171c1c;
  border: 1px solid var(--lc-widget-border);
  border-radius: 3px;
  user-select: none;
  -webkit-user-select: none;
}

.lc-toolbar * {
  user-select: inherit;
  -webkit-user-select: inherit;
}

.lc-toolbar-horizontal {
  flex-flow: row wrap;
}

.lc-toolbar-vertical {
  flex-direction: column;
}

.lc-toolbar-small {
  --lc-toolbar-control: 2rem;
  --lc-toolbar-icon: 0.88rem;
  --lc-toolbar-pad-inline: 0.5rem;
  font-size: 0.74rem;
}

.lc-toolbar-medium {
  font-size: 0.82rem;
}

.lc-toolbar-large {
  --lc-toolbar-control: 2.9rem;
  --lc-toolbar-icon: 1.25rem;
  --lc-toolbar-pad-inline: 0.82rem;
  font-size: 0.9rem;
}

.lc-toolbar-button,
.lc-toolbar-dropdown {
  min-width: 0;
  height: var(--lc-toolbar-control);
  color: var(--lc-widget-text);
  background: #222929;
  border: 1px solid #465355;
  border-radius: 2px;
}

.lc-toolbar-button {
  display: inline-flex;
  width: auto;
  margin: 0;
  padding: 0 var(--lc-toolbar-pad-inline);
  align-items: center;
  justify-content: center;
  font: inherit;
  cursor: pointer;
}

.lc-toolbar-button-icon-only {
  width: var(--lc-toolbar-control);
  padding: 0;
}

.lc-toolbar-button-content,
.lc-toolbar-dropdown {
  display: inline-flex;
  align-items: center;
  gap: 0.45rem;
}

.lc-toolbar-button-label,
.lc-toolbar-dropdown-label {
  line-height: 1;
  white-space: nowrap;
}

.lc-toolbar-dropdown {
  padding-left: var(--lc-toolbar-pad-inline);
  overflow: hidden;
}

.lc-toolbar-dropdown-label {
  color: var(--lc-widget-muted);
}

.lc-toolbar-icon {
  width: var(--lc-toolbar-icon);
  height: var(--lc-toolbar-icon);
  flex: 0 0 auto;
}

.lc-toolbar-select {
  width: auto;
  min-width: 5.5rem;
  height: 100%;
  margin: 0;
  padding: 0 1.9rem 0 0.15rem;
  color: var(--lc-widget-heading);
  background: transparent;
  border: 0;
  border-radius: 0;
  font: inherit;
  cursor: pointer;
}

.lc-toolbar-separator {
  width: 1px;
  margin: 0.3rem 0.12rem;
  align-self: stretch;
  background: #465355;
}

.lc-toolbar-vertical .lc-toolbar-separator {
  width: auto;
  height: 1px;
  margin: 0.12rem 0.3rem;
}

.lc-toolbar-button:hover,
.lc-toolbar-dropdown:hover {
  color: #fff;
  background: #303838;
  border-color: #708083;
}

.lc-toolbar-button:focus-visible {
  outline: 2px solid rgb(78 181 222 / 55%);
  outline-offset: 1px;
}

.lc-toolbar-select:focus-visible {
  outline: none;
}

.lc-toolbar-dropdown:focus-within {
  border-color: var(--lc-widget-focus);
  box-shadow: 0 0 0 1px rgb(78 181 222 / 45%);
}

.lc-toggle-line {
  display: flex;
  min-height: 3rem;
  padding: 0.65rem 0.8rem;
  align-items: center;
  justify-content: space-between;
  gap: 1rem;
  background: #1d2222;
  border: 1px solid #465355;
  border-radius: 3px;
}

.lc-toggle-line input {
  accent-color: var(--lc-widget-focus);
}

.lc-status-chip {
  display: inline-flex;
  width: max-content;
  padding: 0.3rem 0.55rem;
  align-items: center;
  color: #c9fff4;
  background: rgb(26 188 156 / 14%);
  border: 1px solid rgb(26 188 156 / 55%);
  border-radius: 999px;
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.75rem;
}

.lc-progress {
  width: 100%;
  height: 0.72rem;
  overflow: hidden;
  accent-color: var(--lc-widget-focus);
  background: #141919;
  border: 0;
  border-radius: 999px;
}

.lc-progress::-webkit-progress-bar {
  background: #141919;
  border-radius: 999px;
}

.lc-progress::-webkit-progress-value {
  background: var(--lc-widget-focus);
  border-radius: 999px;
}

.lc-summary-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 0.6rem;
}

.lc-summary-grid .lc-widget-output {
  display: grid;
  align-content: center;
  justify-content: stretch;
  gap: 0.2rem;
}

.lc-console {
  display: grid;
  min-width: 0;
  height: 14rem;
  min-height: 10rem;
  grid-template-rows: 2.35rem minmax(0, 1fr);
  overflow: hidden;
  background: #141919;
  border: 1px solid var(--lc-widget-border);
  border-radius: 3px;
}

.lc-console-header {
  display: flex;
  min-width: 0;
  padding: 0 0.75rem;
  align-items: center;
  justify-content: space-between;
  gap: 1rem;
  background: #222929;
  border-bottom: 1px solid #465355;
}

.lc-console-title {
  overflow: hidden;
  color: var(--lc-widget-heading);
  font-size: 0.76rem;
  font-weight: 700;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.lc-console-status {
  color: var(--lc-widget-focus);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.68rem;
  white-space: nowrap;
}

.lc-console-viewport {
  min-width: 0;
  min-height: 0;
  padding: 0.65rem 0.75rem 0.8rem;
  overflow: auto;
  caret-color: transparent;
  color: var(--lc-widget-text);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: clamp(0.72rem, 1.6vw, 0.82rem);
  font-variant-numeric: tabular-nums;
  line-height: 1.5;
  scrollbar-color: #657477 #141919;
  scrollbar-width: thin;
}

.lc-console-viewport:focus-visible {
  outline: 1px solid rgb(78 181 222 / 55%);
  outline-offset: -1px;
}

.lc-console-lines {
  display: grid;
  align-content: start;
  gap: 0.08rem;
}

.lc-console-line {
  display: grid;
  min-width: 0;
  grid-template-columns: 4.25rem minmax(0, 1fr);
  column-gap: 0.7rem;
}

.lc-console-prefix {
  color: var(--lc-widget-muted);
  text-align: right;
  white-space: pre;
}

.lc-console-message {
  min-width: 0;
  overflow-wrap: anywhere;
  white-space: pre-wrap;
}

.lc-console-input .lc-console-prefix,
.lc-console-broker .lc-console-prefix {
  color: var(--lc-widget-link);
}

.lc-console-info .lc-console-prefix {
  color: var(--lc-widget-focus);
}

.lc-console-warning .lc-console-prefix {
  color: #f3bf61;
}

.lc-console-error .lc-console-prefix {
  color: #f07178;
}

.lc-console-empty {
  color: var(--lc-widget-muted);
  font-style: italic;
}

.lc-tabs-demo {
  min-width: 0;
  height: 13rem;
  overflow: hidden;
  border: 1px solid var(--lc-widget-border);
  border-radius: 3px;
}

.lc-tab-sample {
  display: grid;
  height: 100%;
  min-width: 0;
  padding: clamp(0.9rem, 2.5vw, 1.3rem);
  align-content: start;
  gap: 0.65rem;
  overflow: auto;
}

.lc-tab-sample h2 {
  margin: 0;
  color: var(--lc-widget-heading);
  font-size: 1.05rem;
}

.lc-tab-sample p {
  max-width: 58ch;
  margin: 0;
  color: var(--lc-widget-muted);
  font-size: 0.82rem;
  line-height: 1.5;
}

.lc-tab-sample .lc-widget-output {
  margin-top: 0.25rem;
}

@media (max-width: 480px) {
  .lc-widget-header {
    align-items: start;
    flex-direction: column;
    gap: 0.35rem;
  }

  .lc-control-row,
  .lc-summary-grid {
    grid-template-columns: 1fr;
  }

  .lc-toolbar-showcase {
    grid-template-columns: 1fr;
  }

  .lc-toolbar-example-wide {
    grid-column: auto;
  }
}
"""

function widget_header(kicker, title)
    return DOM.header(
        DOM.div(
            DOM.div(kicker; class="lc-widget-kicker"),
            DOM.h1(title)
        ),
        DOM.span("session-local"; class="lc-widget-session");
        class="lc-widget-header"
    )
end

function widget_shell(kicker, title, content)
    return DOM.div(
        DOM.style(WIDGET_THEME),
        widget_header(kicker, title),
        DOM.div(content; class="lc-widget-content");
        class="lc-widget-shell"
    )
end

function widget_output(label, value)
    return DOM.div(
        DOM.span(label),
        DOM.output(value; class="lc-widget-value");
        class="lc-widget-output"
    )
end

function control_copy(label, hint)
    return DOM.div(
        DOM.span(label; class="lc-control-label"),
        DOM.p(hint; class="lc-control-hint");
        class="lc-control-copy"
    )
end

function slider_widget()
    return App(; title="Slider · LineCableModels playground") do session
        slider = Slider(
            collect(0.0:0.5:100.0);
            value=35.0,
            class="lc-native-range"
        )
        value = map(session, slider.value) do current
            return "$(round(current; digits=1)) %"
        end
        content = DOM.div(
            DOM.div(
                control_copy(
                    "Underground-cable share",
                    "A continuous selector with a reactive Julia value."
                ),
                slider;
                class="lc-control-row"
            ),
            widget_output("Current value", value);
            class="lc-control-stack"
        )
        return widget_shell("SELECTOR", "Native range slider", content)
    end
end

function toggle_widget()
    return App(; title="Toggle · LineCableModels playground") do session
        toggle = Checkbox(true)
        state = map(session, toggle.value) do enabled
            return enabled ? "enabled" : "disabled"
        end
        content = DOM.div(
            DOM.div(
                DOM.div(
                    DOM.span("Include earth return"; class="lc-control-label"),
                    DOM.p(
                        "A Boolean selector backed by an Observable{Bool}.";
                        class="lc-control-hint"
                    )
                ),
                toggle;
                class="lc-toggle-line"
            ),
            DOM.span(state; class="lc-status-chip");
            class="lc-control-stack"
        )
        return widget_shell("SELECTOR", "Checkbox and status", content)
    end
end

function dropdown_widget()
    return App(; title="Dropdown · LineCableModels playground") do session
        options = ["Solid conductor", "Concentric strands", "Tubular sheath"]
        dropdown = Dropdown(options; index=2, class="lc-native-select")
        selection = map(session, dropdown.value) do current
            return string(current)
        end
        content = DOM.div(
            DOM.div(
                control_copy(
                    "Conductor construction",
                    "Select one declared domain-model alternative."
                ),
                dropdown;
                class="lc-control-row"
            ),
            widget_output("Selected", selection);
            class="lc-control-stack"
        )
        return widget_shell("CHOICE", "Dropdown selector", content)
    end
end

function text_widget()
    return App(; title="Text input · LineCableModels playground") do session
        field = TextField(
            "B4–B5 cable corridor";
            class="lc-native-field",
            aria_label="Scenario name"
        )
        summary = map(session, field.value) do current
            isempty(strip(current)) ? "untitled scenario" : current
        end
        content = DOM.div(
            DOM.div(
                control_copy(
                    "Scenario name",
                    "Text is committed on change and remains local to this session."
                ),
                field;
                class="lc-control-row"
            ),
            widget_output("Resolved label", summary);
            class="lc-control-stack"
        )
        return widget_shell("INPUT", "Text field", content)
    end
end

function number_spinner_widget()
    return App(; title="Number spinner · LineCableModels playground") do session
        spacing = NumberInput(
            1.0;
            style=nothing,
            class="lc-native-number",
            min=0.25,
            max=10.0,
            step=0.25,
            var"aria-label"="Cable spacing in metres"
        )
        reset = Button("Redefine selector")
        on(reset.value) do _
            spacing.value[] = 1.0
            return nothing
        end
        resolved = map(session, spacing.value) do value
            return "$(round(value; digits=2)) m"
        end
        content = DOM.div(
            DOM.div(
                control_copy(
                    "Lateral spacing",
                    "Type a value or use the native increment and decrement controls."
                ),
                spacing;
                class="lc-control-row"
            ),
            widget_output("Resolved spacing", resolved),
            DOM.div(reset; class="lc-widget-actions");
            class="lc-control-stack"
        )
        return widget_shell("NUMERIC INPUT", "Number spinner", content)
    end
end

function actions_widget()
    return App(; title="Actions · LineCableModels playground") do session
        count = Observable(0)
        advance = Button("Add sample")
        reset = Button("Reset")
        on(advance.value) do _
            count[] += 1
            return nothing
        end
        on(reset.value) do _
            count[] = 0
            return nothing
        end
        label = map(session, count) do current
            return "$current captured"
        end
        content = DOM.div(
            DOM.div(advance, reset; class="lc-widget-actions"),
            widget_output("Session counter", label);
            class="lc-control-stack"
        )
        return widget_shell("ACTION", "Buttons and state", content)
    end
end

function progress_widget()
    return App(; title="Progress · LineCableModels playground") do session
        progress = Observable(25)
        advance = Button("Advance stage")
        reset = Button("Reset")
        on(advance.value) do _
            progress[] = min(progress[] + 25, 100)
            return nothing
        end
        on(reset.value) do _
            progress[] = 0
            return nothing
        end
        stage = map(session, progress) do current
            current == 0 && return "queued · 0 %"
            current < 50 && return "constructing · $current %"
            current < 100 && return "computing · $current %"
            return "ready · 100 %"
        end
        content = DOM.div(
            DOM.progress(; value=progress, max=100, class="lc-progress"),
            DOM.span(stage; class="lc-status-chip"),
            DOM.div(advance, reset; class="lc-widget-actions");
            class="lc-control-stack"
        )
        return widget_shell("STATUS", "Progress and lifecycle", content)
    end
end

function console_widget()
    return App(; title="Console · LineCableModels playground") do _
        console = ConsoleView(
            entries=[
                ConsoleEntry(:input, "using LineCableModels"),
                ConsoleEntry(:info, "Publisher session connected"),
                ConsoleEntry(:broker, "worker.status → ready"),
                ConsoleEntry(:stdout, "LineCableModels playground ready"),
            ];
            title="Julia and worker messages",
            status="live",
            max_entries=80
        )
        sample_index = Ref(0)
        samples = [
            (:input, "dispatch(:line_parameters; frequencies=0:50:2500)"),
            (:broker, "line_parameters.queued · request 7f2a"),
            (:info, "Worker accepted frequency sweep"),
            (:stdout, "result.status = :ready"),
            (:warning, "Static gallery: no numerical job was executed"),
        ]
        append_message = Button("Append message")
        clear_messages = Button("Clear")
        on(append_message.value) do _
            sample_index[] = mod1(sample_index[] + 1, length(samples))
            channel, message = samples[sample_index[]]
            append_console!(console, message; channel)
            set_console_status!(console, "$(length(console.entries[])) lines")
            return nothing
        end
        on(clear_messages.value) do _
            clear_console!(console)
            set_console_status!(console, "empty")
            return nothing
        end
        content = DOM.div(
            console,
            DOM.div(append_message, clear_messages; class="lc-widget-actions");
            class="lc-control-stack"
        )
        return widget_shell("OUTPUT", "Console viewport", content)
    end
end

function playground_widgets_theme()
    return BonitoWidgets.Theme(
        ;
        scheme=:dark,
        bg="#252b2b",
        bg_panel="#1d2222",
        bg_bar="#222929",
        bg_hover="#303838",
        text="#f2f2f2",
        text_muted="#a8b0b0",
        accent="#1abc9c",
        accent_bg="rgba(26, 188, 156, 0.14)",
        border="#5e6d6f",
        radius="3px",
        radius_sm="2px",
        font="Lato, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
        font_size="0.84rem",
        font_size_sm="0.76rem",
        bar_size="2.55rem"
    )
end

function tabs_widget()
    return App(; title="Tabs · LineCableModels playground") do session
        labels = ["Summary", "Parameters", "Messages"]
        tabs = BonitoWidgets.Tabs(
            "Summary" => DOM.div(
                DOM.h2("Persistent content"),
                DOM.p(
                    "All panels stay mounted. Switching tabs changes visibility " *
                    "without rebuilding figures, controls, or session state."
                ),
                widget_output("Calculation state", "not dispatched");
                class="lc-tab-sample"
            ),
            "Parameters" => DOM.div(
                DOM.h2("Declared inputs"),
                DOM.p(
                    "Geometry, material, and sweep controls can share this panel " *
                    "without becoming part of the document layout."
                ),
                widget_output("Earth resistivity", "100 Ω·m");
                class="lc-tab-sample"
            ),
            "Messages" => DOM.div(
                DOM.h2("Worker boundary"),
                DOM.p(
                    "A console or broker status view can remain alive here while " *
                    "another panel is visible."
                ),
                widget_output("Transport", "NATS · idle");
                class="lc-tab-sample"
            );
            active=1,
            closable=false,
            style=Styles("height" => "100%")
        )
        active_label = map(session, tabs.active) do index
            return labels[index]
        end
        content = DOM.div(
            playground_widgets_theme(),
            DOM.div(tabs; class="lc-tabs-demo"),
            widget_output("Active panel", active_label);
            class="lc-control-stack"
        )
        return widget_shell("LAYOUT", "Persistent tabs", content)
    end
end

function control_panel_widget()
    return App(; title="Control panel · LineCableModels playground") do session
        frequency = Slider(
            collect(50:50:5000);
            value=1000,
            class="lc-native-range"
        )
        model = Dropdown(
            ["Default earth", "Carson", "Pollaczek"];
            index=1,
            class="lc-native-select"
        )
        include_losses = Checkbox(true)
        reset = Button("Redefine selectors")
        on(reset.value) do _
            frequency[] = 1000
            model.option_index[] = 1
            include_losses.value[] = true
            return nothing
        end

        frequency_label = map(session, frequency.value) do current
            current >= 1000 ? "$(round(current / 1000; digits=2)) kHz" : "$current Hz"
        end
        model_label = map(session, model.value) do current
            return string(current)
        end
        loss_label = map(session, include_losses.value) do enabled
            return enabled ? "included" : "ignored"
        end

        content = DOM.div(
            DOM.div(
                control_copy("Frequency", "Immediate visual selector."),
                frequency;
                class="lc-control-row"
            ),
            DOM.div(
                control_copy("Earth model", "Discrete calculation strategy."),
                model;
                class="lc-control-row"
            ),
            DOM.div(
                DOM.span("Conductor losses"; class="lc-control-label"),
                include_losses;
                class="lc-toggle-line"
            ),
            DOM.div(
                widget_output("Frequency", frequency_label),
                widget_output("Model", model_label),
                widget_output("Losses", loss_label),
                widget_output("Execution", "not dispatched");
                class="lc-summary-grid"
            ),
            DOM.div(reset; class="lc-widget-actions");
            class="lc-control-stack"
        )
        return widget_shell("COMPOSITION", "Engineering control panel", content)
    end
end

module ToolbarDemoCallbacks

function handle_toolbar(status, event)
    value = isnothing(event.value) ? "" : " = $(event.value)"
    status[] = "$(event.namespace).$(event.action)$value"
    return nothing
end

end

function toolbar_widget()
    return App(; title="Toolbar · LineCableModels playground") do session
        status = Observable("No action dispatched")
        figure_binding = ToolbarBinding(
            ToolbarDemoCallbacks,
            status;
            namespace=:figure_tools
        )
        navigation_binding = ToolbarBinding(
            ToolbarDemoCallbacks,
            status;
            namespace=:page_navigation
        )

        figure_tools = toolbar(
            [
                ToolbarButton(:reset_view; icon=:reset, tooltip="Reset view"),
                ToolbarButton(:fit_view; icon=:fit, label="Fit view"),
                ToolbarButton(:zoom_in; icon=:zoom_in, tooltip="Zoom in"),
                ToolbarButton(:zoom_out; icon=:zoom_out, tooltip="Zoom out"),
                ToolbarSeparator(),
                ToolbarDropdown(
                    :display_layer,
                    [
                        :nominal => "Nominal",
                        :envelope => "Envelope",
                        :results => "Results",
                    ];
                    icon=:layers,
                    label="Layer",
                    selected=:nominal
                ),
                ToolbarButton(:export_view; icon=:download, label="Export"),
            ];
            binding=figure_binding,
            orientation=:horizontal,
            size=:medium,
            label="Figure tools"
        )

        navigation_tools = toolbar(
            [
                ToolbarButton(:previous_page; icon=:previous, tooltip="Previous page"),
                ToolbarButton(:run_page; icon=:play, tooltip="Run page action"),
                ToolbarButton(:next_page; icon=:next, tooltip="Next page"),
            ];
            binding=navigation_binding,
            orientation=:vertical,
            size=:small,
            label="Page navigation"
        )

        presentation_tools = toolbar(
            [
                ToolbarButton(:previous_page; icon=:previous, label="Previous"),
                ToolbarButton(:run_page; icon=:play, label="Run"),
                ToolbarButton(:next_page; icon=:next, label="Next"),
            ];
            binding=navigation_binding,
            orientation=:horizontal,
            size=:large,
            label="Presentation actions"
        )

        content = DOM.div(
            DOM.div(
                DOM.span("HORIZONTAL · MEDIUM"; class="lc-toolbar-example-label"),
                figure_tools;
                class="lc-toolbar-example"
            ),
            DOM.div(
                DOM.span("VERTICAL · SMALL"; class="lc-toolbar-example-label"),
                navigation_tools;
                class="lc-toolbar-example"
            ),
            DOM.div(
                DOM.span("HORIZONTAL · LARGE"; class="lc-toolbar-example-label"),
                presentation_tools;
                class="lc-toolbar-example lc-toolbar-example-wide"
            ),
            widget_output("Last callback", status);
            class="lc-toolbar-showcase"
        )
        return widget_shell("ACTION SURFACE", "Namespaced toolbar", content)
    end
end

const WIDGET_ROUTES = (
    "/widgets/slider" => slider_widget,
    "/widgets/toggle" => toggle_widget,
    "/widgets/dropdown" => dropdown_widget,
    "/widgets/text-input" => text_widget,
    "/widgets/number-spinner" => number_spinner_widget,
    "/widgets/actions" => actions_widget,
    "/widgets/progress" => progress_widget,
    "/widgets/console" => console_widget,
    "/widgets/tabs" => tabs_widget,
    "/widgets/toolbar" => toolbar_widget,
    "/widgets/control-panel" => control_panel_widget,
)

function register_widget_routes!(server)
    for (route, factory) in WIDGET_ROUTES
        Bonito.route!(server, route => factory())
    end
    return server
end
