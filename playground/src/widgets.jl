const WIDGET_THEME = raw"""
html,
body {
  width: 100%;
  min-width: 0;
  min-height: 100%;
  margin: 0;
  background: var(--lc-widget-bg);
  caret-color: transparent;
  cursor: default;
  font-family: Lato, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
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
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-border-soft);
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
  background: var(--lc-option-bg) !important;
  border: 1px solid var(--lc-widget-border) !important;
  border-radius: 3px !important;
  box-shadow: none !important;
}

.lc-native-number {
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
  color: var(--lc-accent-ink) !important;
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

.lc-widget-actions button:disabled {
  color: var(--lc-widget-muted) !important;
  background: var(--lc-code-bg) !important;
  border-color: var(--lc-border-soft) !important;
  cursor: not-allowed;
  opacity: 0.62;
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
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-border-soft);
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

.lc-ribbon-showcase {
  display: grid;
  min-width: 0;
  gap: 0.8rem;
}

.lc-ribbon-demo-body {
  display: grid;
  min-height: 11rem;
  grid-template-columns: minmax(0, 1fr) minmax(12rem, 0.28fr);
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-border-soft);
  border-radius: 3px;
}

.lc-ribbon-demo-workspace,
.lc-ribbon-demo-states {
  min-width: 0;
  padding: 1rem;
}

.lc-ribbon-demo-workspace {
  display: grid;
  align-content: center;
}

.lc-ribbon-demo-workspace h2 {
  margin: 0.45rem 0 0.5rem;
  color: var(--lc-widget-heading);
  font-size: 1.15rem;
}

.lc-ribbon-demo-workspace p {
  max-width: 48rem;
  margin: 0;
  color: var(--lc-widget-muted);
  font-size: 0.82rem;
  line-height: 1.5;
}

.lc-ribbon-demo-states {
  border-left: 1px solid var(--lc-border-soft);
}

.lc-ribbon-demo-states dl {
  display: grid;
  margin: 0.5rem 0 0;
  grid-template-columns: auto minmax(0, 1fr);
  font-size: 0.72rem;
}

.lc-ribbon-demo-states dt,
.lc-ribbon-demo-states dd {
  margin: 0;
  padding: 0.32rem 0;
  border-top: 1px solid var(--lc-border-soft);
}

.lc-ribbon-demo-states dt {
  padding-right: 0.7rem;
  color: var(--lc-widget-heading);
  font-weight: 700;
}

.lc-ribbon-demo-states dd {
  color: var(--lc-widget-muted);
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
  background: var(--lc-console-bg);
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
.lc-toolbar-dropdown,
.lc-toolbar-toggle,
.lc-toolbar-number {
  min-width: 0;
  height: var(--lc-toolbar-control);
  color: var(--lc-widget-text);
  background: var(--lc-toolbar-bg);
  border: 1px solid var(--lc-border-soft);
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
.lc-toolbar-dropdown,
.lc-toolbar-toggle {
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

.lc-toolbar-toggle {
  display: inline-flex;
  padding: 0 var(--lc-toolbar-pad-inline);
  align-items: center;
  gap: 0.45rem;
  cursor: pointer;
}

.lc-toolbar-toggle input {
  width: 0.95rem;
  height: 0.95rem;
  margin: 0 0 0 0.15rem;
  accent-color: var(--lc-widget-focus);
}

.lc-toolbar-toggle-icon {
  display: inline-flex;
}

.lc-toolbar-number {
  display: grid;
  height: auto;
  min-height: var(--lc-toolbar-control);
  padding: 0.24rem var(--lc-toolbar-pad-inline);
  align-content: center;
  gap: 0.18rem;
}

.lc-toolbar-number-label {
  color: var(--lc-widget-muted);
  font-size: 0.68em;
  line-height: 1;
}

.lc-toolbar-number-field {
  display: flex;
  min-width: 0;
  align-items: center;
}

.lc-toolbar-number-input {
  width: 5rem;
  min-width: 0;
  padding: 0;
  color: var(--lc-widget-heading);
  background: transparent;
  border: 0;
  outline: none;
  font: inherit;
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
}

.lc-toolbar-number-unit {
  color: var(--lc-widget-muted);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.76em;
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
  background: var(--lc-border-soft);
}

.lc-toolbar-vertical .lc-toolbar-separator {
  width: auto;
  height: 1px;
  margin: 0.12rem 0.3rem;
}

.lc-toolbar-button:hover,
.lc-toolbar-dropdown:hover,
.lc-toolbar-toggle:hover,
.lc-toolbar-number:hover {
  color: var(--lc-strong-text);
  background: var(--lc-hover-bg-strong);
  border-color: var(--lc-border);
}

.lc-toolbar-button-active,
.lc-toolbar-button-active:hover {
  color: var(--lc-accent-ink);
  background: var(--lc-widget-focus);
  border-color: var(--lc-widget-focus);
}

.lc-toolbar-button-busy {
  position: relative;
  padding-right: calc(var(--lc-toolbar-pad-inline) + 1rem);
}

.lc-toolbar-button-busy::after {
  position: absolute;
  width: 0.7rem;
  height: 0.7rem;
  right: 0.45rem;
  border: 1px solid currentColor;
  border-right-color: transparent;
  border-radius: 50%;
  content: "";
  animation: lc-toolbar-spin 0.8s linear infinite;
}

.lc-toolbar-control-disabled,
.lc-toolbar-button:disabled {
  cursor: not-allowed;
  opacity: 0.46;
}

@keyframes lc-toolbar-spin {
  to { transform: rotate(360deg); }
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
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-border-soft);
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
  color: var(--lc-success-text);
  background: var(--lc-success-bg);
  border: 1px solid var(--lc-success-border);
  border-radius: 999px;
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.75rem;
}

.lc-progress {
  width: 100%;
  height: 0.72rem;
  overflow: hidden;
  accent-color: var(--lc-widget-focus);
  background: var(--lc-console-bg);
  border: 0;
  border-radius: 999px;
}

.lc-progress::-webkit-progress-bar {
  background: var(--lc-console-bg);
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
  background: var(--lc-console-bg);
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
  background: var(--lc-toolbar-bg);
  border-bottom: 1px solid var(--lc-border-soft);
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
  scrollbar-color: var(--lc-scrollbar) var(--lc-console-bg);
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
  color: var(--lc-warning);
}

.lc-console-error .lc-console-prefix {
  color: var(--lc-danger);
}

.lc-console-empty {
  color: var(--lc-widget-muted);
  font-style: italic;
}

.lc-job-panel {
  display: grid;
  min-width: 0;
  overflow: hidden;
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-widget-border);
  border-radius: 3px;
}

.lc-job-header,
.lc-job-status-line {
  display: flex;
  min-width: 0;
  align-items: center;
  justify-content: space-between;
  gap: 1rem;
}

.lc-job-header {
  padding: 0.8rem 0.9rem;
  background: var(--lc-toolbar-bg);
  border-bottom: 1px solid var(--lc-border-soft);
}

.lc-job-heading h2 {
  margin: var(--lc-kicker-heading-gap) 0 0;
  color: var(--lc-widget-heading);
  font-size: 1rem;
  line-height: 1.2;
}

.lc-worker-state,
.lc-job-state,
.lc-job-stage,
.lc-job-result-status {
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.72rem;
}

.lc-worker-state-online {
  color: var(--lc-widget-focus);
}

.lc-worker-state-offline,
.lc-worker-state-connecting,
.lc-worker-state-degraded {
  color: var(--lc-warning);
}

.lc-job-body {
  display: grid;
  min-width: 0;
  padding: 0.9rem;
  gap: 0.75rem;
}

.lc-job-state {
  color: var(--lc-widget-link);
  font-weight: 700;
}

.lc-job-stage,
.lc-job-result-status {
  overflow: hidden;
  color: var(--lc-widget-muted);
  text-overflow: ellipsis;
  white-space: nowrap;
}

.lc-job-progress {
  position: relative;
  height: 0.8rem;
  overflow: hidden;
  background: var(--lc-console-bg);
  border-radius: 999px;
}

.lc-job-progress-fill {
  position: absolute;
  inset: 0 auto 0 0;
  width: var(--lc-job-progress);
  background: var(--lc-widget-focus);
  border-radius: inherit;
  transition: width 120ms linear;
}

.lc-job-progress-value {
  position: absolute;
  inset: 50% 0 auto;
  color: var(--lc-strong-text);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.58rem;
  line-height: 1;
  text-align: center;
  transform: translateY(-50%);
}

.lc-job-actions button:first-child {
  color: var(--lc-accent-ink) !important;
  background: var(--lc-widget-focus) !important;
}

.lc-job-result {
  display: grid;
  min-width: 0;
  gap: 0.45rem;
}

.lc-job-result-preview {
  max-height: 9rem;
  margin: 0;
  padding: 0.65rem 0.75rem;
  overflow: auto;
  color: var(--lc-widget-text);
  background: var(--lc-console-bg);
  border: 1px solid var(--lc-border-soft);
  border-radius: 3px;
  font-size: 0.72rem;
  line-height: 1.45;
  white-space: pre-wrap;
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

.lc-repeater-demo {
  display: grid;
  min-width: 0;
  align-content: start;
  gap: 0.75rem;
}

.lc-repeater-demo-fields {
  display: grid;
  min-width: 0;
  grid-template-columns:
    minmax(9rem, 1.15fr)
    minmax(7rem, 0.72fr)
    minmax(8rem, 0.9fr)
    minmax(9rem, 1fr)
    auto;
  align-items: end;
  gap: 0.6rem;
}

.lc-repeater-demo-field {
  display: grid;
  min-width: 0;
  align-content: start;
  gap: 0.3rem;
}

.lc-repeater-demo-field > span {
  color: var(--lc-widget-muted);
  font-size: 0.69rem;
  line-height: 1.2;
}

.lc-repeater-demo-field > input:not([type="checkbox"]),
.lc-repeater-demo-field > select {
  width: 100%;
  min-width: 0 !important;
  min-height: 2.25rem;
  margin: 0 !important;
  padding: 0.42rem 0.58rem !important;
  color: var(--lc-widget-text) !important;
  background: var(--lc-option-bg) !important;
  border: 1px solid var(--lc-widget-border) !important;
  border-radius: 2px !important;
  box-shadow: none !important;
}

.lc-repeater-demo-enabled {
  display: flex;
  min-height: 2.25rem;
  padding: 0 0.55rem;
  align-items: center;
  gap: 0.45rem;
  color: var(--lc-widget-text);
  background: var(--lc-option-bg);
  border: 1px solid var(--lc-widget-border);
  border-radius: 2px;
  font-size: 0.76rem;
  white-space: nowrap;
}

.lc-repeater-demo-enabled input {
  margin: 0;
  accent-color: var(--lc-widget-focus);
}

.lc-repeater-snapshot {
  max-height: 7rem;
  margin: 0;
  padding: 0.65rem 0.75rem;
  overflow: auto;
  color: var(--lc-widget-link);
  background: var(--lc-console-bg);
  border: 1px solid var(--lc-border-soft);
  border-radius: 3px;
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.7rem;
  line-height: 1.45;
  white-space: pre-wrap;
}

.lc-toolkit-split,
.lc-feedback-specimen,
.lc-data-specimen {
  display: grid;
  min-width: 0;
  align-content: start;
  gap: 0.8rem;
}

.lc-toolkit-split {
  grid-template-columns: minmax(0, 1.35fr) minmax(13rem, 0.65fr);
}

.lc-toolkit-summary,
.lc-specimen-panel {
  display: grid;
  min-width: 0;
  padding: 0.75rem;
  align-content: start;
  gap: 0.65rem;
  background: var(--lc-sunken-bg);
  border: 1px solid var(--lc-border-soft);
  border-radius: var(--lc-radius);
}

.lc-specimen-panel > h2,
.lc-toolkit-summary > h2 {
  margin: 0;
  color: var(--lc-heading);
  font-size: 0.92rem;
}

.lc-toolkit-summary pre {
  min-height: 8rem;
  margin: 0;
  padding: 0.65rem;
  overflow: auto;
  color: var(--lc-link);
  background: var(--lc-console-bg);
  border: 1px solid var(--lc-border-soft);
  font-family: "JuliaMono", "SFMono-Regular", Consolas, monospace;
  font-size: 0.67rem;
  line-height: 1.45;
  white-space: pre-wrap;
}

.lc-data-specimen-grid {
  display: grid;
  min-width: 0;
  grid-template-columns: minmax(0, 1.25fr) minmax(14rem, 0.75fr);
  gap: 0.8rem;
}

.lc-data-specimen .lc-viewport-frame {
  min-height: 19rem;
}

@media (max-width: 480px) {
  .lc-widget-header {
    align-items: start;
    flex-direction: column;
    gap: 0.35rem;
  }

  .lc-job-header,
  .lc-job-status-line {
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

  .lc-ribbon-demo-body {
    grid-template-columns: 1fr;
  }

  .lc-ribbon-demo-states {
    border-top: 1px solid var(--lc-border-soft);
    border-left: 0;
  }

  .lc-repeater-demo-fields {
    grid-template-columns: 1fr;
  }

  .lc-toolkit-split,
  .lc-data-specimen-grid {
    grid-template-columns: 1fr;
  }
}

@media (min-width: 481px) and (max-width: 900px) {
  .lc-repeater-demo-fields {
    grid-template-columns: repeat(2, minmax(0, 1fr));
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

function widget_theme_script()
    return DOM.script(js"""
    (() => {
        const storageKey = 'lcm.playground.theme';
        const legacyStorageKey = 'lcm.workbench.theme';
        const validTheme = value => ['system', 'dark', 'light'].includes(value);
        const systemTheme = window.matchMedia('(prefers-color-scheme: light)');
        let preference = 'dark';

        const applyTheme = selected => {
            preference = validTheme(selected) ? selected : 'dark';
            const resolved = preference === 'system'
                ? (systemTheme.matches ? 'light' : 'dark')
                : preference;
            document.documentElement.dataset.lcmTheme = preference;
            document.documentElement.dataset.lcmResolvedTheme = resolved;
        };

        try {
            const stored = window.localStorage.getItem(storageKey)
                || window.localStorage.getItem(legacyStorageKey);
            if (validTheme(stored)) preference = stored;
        } catch (_) {
            // Storage is optional for a standalone widget route.
        }
        applyTheme(preference);

        systemTheme.addEventListener('change', () => {
            if (preference === 'system') applyTheme('system');
        });
        window.addEventListener('storage', event => {
            if (![storageKey, legacyStorageKey].includes(event.key)) return;
            let stored = null;
            try {
                stored = window.localStorage.getItem(storageKey)
                    || window.localStorage.getItem(legacyStorageKey);
            } catch (_) {
                // Fall through to the deterministic dark default.
            }
            applyTheme(stored);
        });
    })();
    """)
end

function widget_shell(kicker, title, content)
    return DOM.div(
        DOM.style(BRAND_THEME),
        DOM.style(CONTROL_CONTRACT),
        DOM.style(Toolkit.TOOLKIT_STYLES),
        DOM.style(WIDGET_THEME),
        widget_theme_script(),
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
        dropdown = Dropdown(
            options;
            index=2,
            class="lc-control-select lc-native-select"
        )
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
            class="lc-control-input lc-native-field",
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
            class="lc-control-input lc-native-number",
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
        bg="var(--lc-panel-bg)",
        bg_panel="var(--lc-sunken-bg)",
        bg_bar="var(--lc-toolbar-bg)",
        bg_hover="var(--lc-hover-bg-strong)",
        text="var(--lc-heading)",
        text_muted="var(--lc-muted)",
        accent="var(--lc-focus)",
        accent_bg="var(--lc-success-bg)",
        border="var(--lc-border)",
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
            class="lc-control-select lc-native-select"
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

function ribbon_widget()
    return App(; title="Ribbon · LineCableModels playground") do session
        status = Observable("No command dispatched")
        binding = ToolbarBinding(
            ToolbarDemoCallbacks,
            status;
            namespace=:engineering_ribbon
        )

        project = RibbonGroup(
            "Project",
            ToolbarButton(:new_project; icon=:new_file, label="New"),
            ToolbarButton(:open_project; icon=:open, label="Open"),
            ToolbarButton(:save_project; icon=:save, label="Save");
            size=:medium
        )
        model = RibbonGroup(
            "Model",
            ToolbarButton(
                :select_bus;
                icon=:busbar,
                label="Bus",
                active=true,
                tooltip="Active placement tool"
            ),
            ToolbarButton(:select_cable; icon=:cable, label="Cable"),
            ToolbarToggle(
                :snap_to_grid;
                icon=:grid,
                label="Snap",
                checked=true
            );
            size=:medium
        )
        parameters = RibbonGroup(
            "Parameters",
            ToolbarDropdown(
                :voltage_level,
                [:kv_11 => "11 kV", :kv_33 => "33 kV", :kv_220 => "220 kV"];
                icon=:layers,
                label="Voltage",
                selected=:kv_33
            ),
            ToolbarNumber(
                :line_length;
                label="Length",
                value=12.0,
                minimum=0.1,
                maximum=500.0,
                step=0.1,
                unit="km"
            );
            size=:small
        )
        command_states = RibbonGroup(
            "Command states",
            ToolbarButton(
                :inspect_selection;
                icon=:settings,
                label="Inspect"
            ),
            ToolbarButton(
                :compile_topology;
                icon=:chart,
                label="Compiling",
                busy=true
            ),
            ToolbarButton(
                :remote_study;
                icon=:play,
                label="Remote",
                disabled=true
            );
            size=:small
        )

        home = RibbonTab(:home, "Home", project, model, parameters, command_states)
        insert = RibbonTab(
            :insert,
            "Insert",
            RibbonGroup(
                "Network",
                ToolbarButton(:insert_bus; icon=:busbar, label="Busbar"),
                ToolbarButton(:insert_cable; icon=:cable, label="Cable");
                size=:large
            ),
            RibbonGroup(
                "References",
                ToolbarButton(:insert_grid; icon=:grid, label="Grid reference"),
                ToolbarButton(
                    :insert_annotation;
                    icon=:new_file,
                    label="Annotation"
                );
                size=:medium
            )
        )
        view = RibbonTab(
            :view,
            "View",
            RibbonGroup(
                "Viewport",
                ToolbarButton(:fit_view; icon=:fit, label="Fit"),
                ToolbarButton(:reset_view; icon=:reset, label="Reset"),
                ToolbarDropdown(
                    :display_layer,
                    [:nominal => "Nominal", :results => "Results"];
                    icon=:layers,
                    label="Layer"
                );
                size=:medium
            )
        )
        analysis = RibbonTab(
            :analysis,
            "Analysis",
            RibbonGroup(
                "Studies",
                ToolbarButton(:frequency_scan; icon=:chart, label="Frequency scan"),
                ToolbarButton(:run_power_flow; icon=:play, label="Power flow");
                size=:large
            )
        )
        remote = RibbonTab(
            :remote,
            "Remote",
            RibbonGroup(
                "Worker",
                ToolbarButton(:connect_worker; icon=:settings, label="Connect");
                size=:medium
            );
            disabled=true
        )

        application_ribbon = Ribbon(
            home,
            insert,
            view,
            analysis,
            remote;
            binding,
            quick_access=(
                ToolbarButton(:quick_save; icon=:save, tooltip="Save project"),
                ToolbarButton(:quick_fit; icon=:fit, tooltip="Fit workspace"),
            ),
            label="Engineering workbench commands"
        )
        content = DOM.div(
            application_ribbon,
            DOM.div(
                DOM.div(
                    DOM.span("WORK AREA"; class="lc-toolbar-example-label"),
                    DOM.h2("Cable corridor study"),
                    DOM.p(
                        "The ribbon changes presentation and emits namespaced commands; this mock working area remains mounted."
                    );
                    class="lc-ribbon-demo-workspace"
                ),
                DOM.aside(
                    DOM.span("STATE COVERAGE"; class="lc-toolbar-example-label"),
                    DOM.dl(
                        DOM.dt("Normal"), DOM.dd("New, Open, Cable"),
                        DOM.dt("Hovered"), DOM.dd("Pointer-driven"),
                        DOM.dt("Active"), DOM.dd("Bus"),
                        DOM.dt("Busy"), DOM.dd("Compiling"),
                        DOM.dt("Disabled"), DOM.dd("Remote")
                    );
                    class="lc-ribbon-demo-states"
                );
                class="lc-ribbon-demo-body"
            ),
            widget_output("Last callback", status);
            class="lc-ribbon-showcase"
        )
        return widget_shell("COMMAND COMPOSITION", "Compact application ribbon", content)
    end
end

function repeater_demo_entry(values=(
        name="Cable 01",
        spacing=1.0,
        construction="Solid core",
        note="baseline",
        enabled=true,
    ))
    constructions = ["Solid core", "Concentric strands", "Tubular conductor"]
    selected = something(findfirst(==(string(values.construction)), constructions), 1)
    name = TextField(
        string(values.name);
        class="lc-control-input lc-native-field",
        var"aria-label"="Cable name"
    )
    spacing = NumberInput(
        Float64(values.spacing);
        style=nothing,
        class="lc-control-input lc-native-number",
        min=0.25,
        max=20.0,
        step=0.25,
        var"aria-label"="Lateral spacing in metres"
    )
    construction = Dropdown(
        constructions;
        index=selected,
        class="lc-control-select lc-native-select"
    )
    note = TextField(
        string(values.note);
        class="lc-control-input lc-native-field",
        var"aria-label"="Row note"
    )
    enabled = Checkbox(Bool(values.enabled))
    fields = (
        name=name.value,
        spacing=spacing.value,
        construction=construction.value,
        note=note.value,
        enabled=enabled.value,
    )
    content = DOM.div(
        DOM.label(
            DOM.span("Name"),
            name;
            class="lc-repeater-demo-field"
        ),
        DOM.label(
            DOM.span("Spacing [m]"),
            spacing;
            class="lc-repeater-demo-field"
        ),
        DOM.label(
            DOM.span("Construction"),
            construction;
            class="lc-repeater-demo-field"
        ),
        DOM.label(
            DOM.span("Note"),
            note;
            class="lc-repeater-demo-field"
        ),
        DOM.label(
            enabled,
            DOM.span("Enabled");
            class="lc-repeater-demo-enabled"
        );
        class="lc-repeater-demo-fields"
    )
    return RepeatEntry(fields, content)
end

function repeater_widget()
    return App(; title="Repeater · LineCableModels playground") do session
        initial = [
            repeater_demo_entry(),
            repeater_demo_entry((
                name="Cable 02",
                spacing=6.0,
                construction="Concentric strands",
                note="comparison",
                enabled=true,
            )),
        ]
        repeater = Repeater(
            repeater_demo_entry;
            initial,
            duplicate=repeater_demo_entry,
            min_entries=1,
            max_entries=8,
            label="Cable definitions",
            add_label="Add cable"
        )
        captured = Observable(repr(snapshot(repeater)))
        capture = Button("Read row values")
        on(session, capture.value) do _
            captured[] = repr(snapshot(repeater))
            return nothing
        end
        content = DOM.div(
            repeater,
            DOM.div(capture; class="lc-widget-actions"),
            DOM.pre(captured; class="lc-repeater-snapshot");
            class="lc-repeater-demo"
        )
        return widget_shell("COLLECTION EDITOR", "Repeatable mixed-control rows", content)
    end
end

function file_upload_widget()
    return App(; title="File upload · LineCableModels playground") do session
        directory = mktempdir(; prefix="lcm-upload-demo-", cleanup=true)
        on(session.on_close) do _
            isdir(directory) && rm(directory; recursive=true, force=true)
            return nothing
        end
        target = UploadTarget(
            DirectoryUploadStore(directory),
            "workbench-input";
            extensions=[".csv", ".json", ".txt"],
            max_bytes=2 * 1024 * 1024,
            retention=:session
        )
        upload = UploadField(
            target;
            label="Input data",
            hint="CSV, JSON, or plain text · maximum 2 MiB"
        )
        committed = map(session, upload.value) do value
            isnothing(value) ? "No committed input" :
                "$(value.original_name) · $(value.bytes) bytes · $(first(value.sha256, 12))…"
        end
        content = DOM.div(
            upload,
            widget_output("Committed value", committed);
            class="lc-control-stack"
        )
        return widget_shell("FILE INGRESS", "Transactional upload slot", content)
    end
end

function geographic_map_widget()
    return App(; title="Geographic map · LineCableModels playground") do session
        directory = mktempdir(; prefix="lcm-map-demo-", cleanup=true)
        on(session.on_close) do _
            isdir(directory) && rm(directory; recursive=true, force=true)
            return nothing
        end
        target = UploadTarget(
            DirectoryUploadStore(directory),
            "route-source";
            extensions=[".kml", ".kmz"],
            max_bytes=16 * 1024 * 1024,
            retention=:session,
            validate=validate_kml_upload
        )
        upload = UploadField(
            target;
            label="Route geometry",
            hint="KML or KMZ · Point and LineString placemarks · maximum 16 MiB"
        )
        viewport = GeographicMap(
            KMLUploadSource(upload);
            tags=[
                "terminal" => "Terminal",
                "joint" => "Cable joint",
                "survey" => "Survey reference",
            ],
            basemap=:osm
        )
        return widget_shell(
            "GEOGRAPHIC INPUT",
            "Line-snapped map references",
            viewport
        )
    end
end

function power_system_specimen_diagram()
    series_device = (
        id="custom:transmission-line",
        name="Transmission line",
        category="network",
        description="Two-terminal overhead transmission line",
        viewBox="-16 -32 32 64",
        width=32,
        height=64,
        svg="""
        <line x1="0" y1="-30" x2="0" y2="30" stroke="black" stroke-width="2"/>
        <path d="M -7 -10 L 0 -16 L 7 -10 M -7 10 L 0 4 L 7 10" fill="none" stroke="black" stroke-width="1"/>
        """,
        terminals=[
            (id="t1", x=0, y=-30, orientation="n"),
            (id="t2", x=0, y=30, orientation="s"),
        ],
        params=[
            (name="length", label="Length", type="number", unit="km", showOnCanvas=true),
            (name="Un", label="Voltage", type="number", unit="kV", showOnCanvas=true),
        ],
        label=(x=11, y=0, anchor="start"),
        source=(kind="inline",),
    )
    cable = (
        id="custom:cable",
        name="Underground cable",
        category="network",
        description="Two-terminal underground power cable",
        viewBox="-16 -32 32 64",
        width=32,
        height=64,
        svg="""
        <line x1="0" y1="-30" x2="0" y2="30" stroke="black" stroke-width="3"/>
        <line x1="-6" y1="-22" x2="-6" y2="22" stroke="black" stroke-width="1"/>
        <line x1="6" y1="-22" x2="6" y2="22" stroke="black" stroke-width="1"/>
        """,
        terminals=[
            (id="t1", x=0, y=-30, orientation="n"),
            (id="t2", x=0, y=30, orientation="s"),
        ],
        params=[
            (name="length", label="Length", type="number", unit="km", showOnCanvas=true),
            (name="Un", label="Voltage", type="number", unit="kV", showOnCanvas=true),
        ],
        label=(x=11, y=0, anchor="start"),
        source=(kind="inline",),
    )

    return (
        version="1",
        meta=(
            title="Industrial feeder specimen",
            description="Editable topology specimen; no numerical model attached",
            labelMode="all",
            labelFontSize=9,
            symbolStandard="iec",
        ),
        elements=[
            (id="GRID", kind="grid-source", name="Grid connection"),
            (
                id="T1",
                kind="transformer-2w",
                name="Main transformer",
                params=(S=100, ratio="220/33 kV"),
            ),
            (
                id="LINE-01",
                kind="custom:transmission-line",
                name="Overhead line",
                params=(length=12.0, Un=33.0),
                color="blue",
            ),
            (
                id="CABLE-01",
                kind="custom:cable",
                name="Underground cable",
                params=(length=1.8, Un=33.0),
                color="green",
            ),
            (
                id="LOAD",
                kind="load",
                name="Industrial load",
                params=(P=24.0, cosphi=0.95),
            ),
        ],
        buses=[
            (id="B-HV", name="220 kV bus", layout=(at=[0, -170], span=220)),
            (id="B-MV-A", name="33 kV bus A", layout=(at=[0, 50], span=220)),
            (id="B-MV-B", name="33 kV bus B", layout=(at=[0, 190], span=220)),
        ],
        junctions=[],
        wires=[
            (id="wire-grid-hv", ends=["GRID.t_bottom", "B-HV"]),
            (id="wire-hv-t1", ends=["B-HV", "T1.t1"]),
            (id="wire-t1-mva", ends=["T1.t2", "B-MV-A"]),
            (id="wire-mva-line", ends=["B-MV-A", "LINE-01.t1"]),
            (id="wire-line-mvb", ends=["LINE-01.t2", "B-MV-B"]),
            (id="wire-mvb-cable", ends=["B-MV-B", "CABLE-01.t1"]),
            (id="wire-cable-load", ends=["CABLE-01.t2", "LOAD.t_top"]),
        ],
        layout=Dict(
            "GRID" => (at=[0, -220],),
            "T1" => (at=[30, 10],),
            "LINE-01" => (at=[0, 120],),
            "CABLE-01" => (at=[0, 260],),
            "LOAD" => (at=[0, 360],),
        ),
        annotations=[],
        customKinds=[series_device, cable],
    )
end

function power_system_canvas_widget()
    return App(; title="Power-system canvas · LineCableModels playground") do _
        canvas = PowerSystemCanvas(power_system_specimen_diagram())
        return widget_shell(
            "DOMAIN CANVAS",
            "Editable single-line topology",
            canvas
        )
    end
end

function form_toolkit_widget()
    return App(; title="Forms · LineCableModels playground") do session
        scenario = TextInput(:scenario; value="North corridor", required=true)
        password = SecretInput(:credential; placeholder="Session credential",
            autocomplete="new-password", required=true)
        notes = TextAreaInput(:notes; value="Nominal study case", rows=3)
        mode = RadioGroup(:mode,
            [:planning => "Planning", :operations => "Operations"];
            selected=:planning, orientation=:horizontal)
        phase = SegmentedControl(:phase,
            [:a => "Phase A", :b => "Phase B", :c => "Phase C"];
            selected=:a)
        model = ComboBox(:earth_model,
            [:default => "Default earth", :carson => "Carson", :pollaczek => "Pollaczek"];
            selected=:default)
        outputs = MultiSelect(:outputs,
            [:r => "Resistance", :l => "Inductance", :c => "Capacitance", :g => "Conductance"];
            selected=[:r, :l], size=4)
        length = UnitNumberInput(:length; value=1.0, minimum=0.1,
            maximum=5000, step=0.1, unit="km", required=true)
        harmonics = RangeInput(:harmonics; lower=1, upper=50,
            minimum=1, maximum=100, step=1, unit="order")
        validator = payload -> begin
            errors = Dict{String,String}()
            isempty(strip(string(get(payload, "scenario", "")))) &&
                (errors["scenario"] = "Provide a scenario name.")
            return errors
        end
        form = Form(
            Field("Scenario", scenario; hint="Project-local display name."),
            Field("Credential", password;
                hint="Never mirrored into an Observable or X-ray metadata."),
            Field("Notes", notes),
            Field("Study mode", mode),
            Field("Reference phase", phase),
            Field("Earth model", model),
            Field("Requested outputs", outputs),
            Field("Line length", length),
            Field("Harmonic interval", harmonics);
            name=:toolkit_form, validator, submit_label="Validate specimen",
            reset_label="Restore defaults", cancel_label="Cancel"
        )
        summary = map(session, form.status, form.submitted) do status, payload
            if isnothing(payload)
                return "status = :$status\nNo payload submitted."
            end
            visible = Dict(key => key == "credential" ? "[redacted]" : value
                for (key, value) in payload)
            return "status = :$status\n" * String(JSON3.write(visible))
        end
        return widget_shell(
            "FORM CONTRACT",
            "Typed fields and explicit submission",
            DOM.div(
                form,
                DOM.aside(DOM.h2("Submission boundary"), DOM.pre(summary);
                    class="lc-toolkit-summary");
                class="lc-toolkit-split"
            )
        )
    end
end

function overlay_toolkit_widget()
    return App(; title="Dialogs and notices · LineCableModels playground") do session
        activity = Observable("No overlay action yet")
        message = MessageDialog(:message_specimen, "Preparation complete",
            "Static inputs have been validated. No numerical work was dispatched.";
            description="Informational modal")
        confirmation = ConfirmDialog(:confirmation_specimen, "Replace project file",
            "The previous upload will be replaced atomically.";
            onconfirm=() -> (activity[] = "Replacement confirmed"), danger=true)
        modal_form = Form(
            Field("Project label", TextInput(:project_label; value="Cable study", required=true)),
            Field("Access token", SecretInput(:access_token; autocomplete="off", required=true));
            name=:modal_form, submit_label="Accept", reset_label="Reset",
            cancel_label="Cancel"
        )
        form_dialog = FormDialog(:form_specimen, "Project details", modal_form;
            description="The same Field and Form primitives compose inside a modal.")
        for dialog in (message, confirmation, form_dialog)
            on(session, dialog.event) do event
                isnothing(event) || (activity[] = "$(event.dialog).$(event.action)")
                return nothing
            end
        end
        notice = InlineNotice(:warning, "Unsaved browser state",
            DOM.span("This specimen has not been attached to a persistent project.");
            dismissible=true)
        toasts = ToastCenter(; capacity=3)
        open_message = Button("Message dialog")
        open_confirm = Button("Confirmation dialog")
        open_form = Button("Form dialog")
        restore_notice = Button("Restore notice")
        add_toast = Button("Push toast")
        on(open_message.value) do _; open_dialog!(message); nothing; end
        on(open_confirm.value) do _; open_dialog!(confirmation); nothing; end
        on(open_form.value) do _; open_dialog!(form_dialog); nothing; end
        on(restore_notice.value) do _; notice.visible[] = true; nothing; end
        on(add_toast.value) do _
            push_toast!(toasts, :success, "Snapshot retained",
                "The toast stack is bounded and session-local."; timeout_ms=4500)
            activity[] = "Toast pushed"
            return nothing
        end
        content = DOM.div(
            notice,
            DOM.div(open_message, open_confirm, open_form, restore_notice, add_toast;
                class="lc-widget-actions"),
            widget_output("Latest event", activity),
            message, confirmation, form_dialog, toasts;
            class="lc-feedback-specimen"
        )
        return widget_shell("FEEDBACK", "Dialogs, notices, and toast lifecycle", content)
    end
end

function data_view_toolkit_widget()
    return App(; title="Data views · LineCableModels playground") do session
        table = DataTable(
            TableColumn(:quantity, "Quantity"),
            TableColumn(:value, "Value"; align=:right,
                format=value -> round(value; digits=4)),
            TableColumn(:unit, "Unit"),
            TableColumn(:state, "State");
            rows=[
                (id="r", quantity="Resistance", value=0.0643, unit="Ω/km", state="ready"),
                (id="l", quantity="Inductance", value=0.0616, unit="mH/km", state="ready"),
                (id="c", quantity="Capacitance", value=0.4762, unit="μF/km", state="ready"),
                (id="g", quantity="Conductance", value=0.0002, unit="μS/km", state="estimated"),
            ],
            row_key=(row, _) -> row.id,
            label="Cable constants"
        )
        properties = PropertyGrid(
            PropertyItem("Frequency", "1.259"; unit="kHz"),
            PropertyItem("Worker", "idle"; state=:muted),
            PropertyItem("Geometry", "validated"; state=:success),
            PropertyItem("Warnings", 1; state=:warning);
            columns=2, label="Result metadata"
        )
        disclosure = Disclosure("Calculation metadata", properties; open=true)
        viewport = ViewportFrame("One-frequency package result", table;
            footer="Persistent table · overlays do not remount content")
        ready = Button("Ready")
        loading = Button("Loading")
        empty = Button("Empty")
        failure = Button("Error")
        on(ready.value) do _; set_viewport_state!(viewport, :ready; message=""); nothing; end
        on(loading.value) do _; set_viewport_state!(viewport, :loading;
            message="Waiting for worker result…"); nothing; end
        on(empty.value) do _; set_viewport_state!(viewport, :empty;
            message="No rows match this request."); nothing; end
        on(failure.value) do _; set_viewport_state!(viewport, :error;
            message="Result artifact is unavailable."); nothing; end
        content = DOM.div(
            DOM.div(viewport, DOM.div(
                DOM.aside(DOM.h2("Viewport states"),
                    DOM.div(ready, loading, empty, failure; class="lc-widget-actions");
                    class="lc-specimen-panel"),
                disclosure;
                class="lc-data-specimen-grid"));
            class="lc-data-specimen"
        )
        return widget_shell("DATA SURFACES", "Persistent viewport and compact records", content)
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
    "/widgets/ribbon" => ribbon_widget,
    "/widgets/control-panel" => control_panel_widget,
    "/widgets/form-toolkit" => form_toolkit_widget,
    "/widgets/overlay-toolkit" => overlay_toolkit_widget,
    "/widgets/data-view-toolkit" => data_view_toolkit_widget,
    "/widgets/repeater" => repeater_widget,
    "/widgets/file-upload" => file_upload_widget,
    "/widgets/geographic-map" => geographic_map_widget,
    "/widgets/power-system-canvas" => power_system_canvas_widget,
)

function runtime_job_widget(client)
    return App(; title="Broker job · LineCableModels playground") do session
        steps = NumberInput(
            8.0;
            class="lc-control-input",
            min=1,
            max=40,
            step=1
        )
        interval = NumberInput(
            0.1;
            class="lc-control-input",
            min=0,
            max=1,
            step=0.05
        )
        parameters = () -> Dict{String,Any}(
            "steps" => round(Int, steps.value[]),
            "interval_seconds" => Float64(interval.value[]),
        )
        panel = JobPanel(
            client;
            operation="system.progress",
            parameters,
            input_observables=(steps.value, interval.value),
            session_id="widget-$(uuid4())",
            title="Diagnostic broker lifecycle"
        )
        controls = DOM.div(
            DOM.div(
                control_copy("Progress steps", "Validated diagnostic workload."),
                steps;
                class="lc-control-row"
            ),
            DOM.div(
                control_copy("Step interval", "Delay per step [s]."),
                interval;
                class="lc-control-row"
            ),
            panel;
            class="lc-control-stack"
        )
        return widget_shell("BROKER BOUNDARY", "Asynchronous job panel", controls)
    end
end

function register_widget_routes!(server, client=nothing)
    for (route, factory) in WIDGET_ROUTES
        Bonito.route!(server, route => factory())
    end
    if !isnothing(client)
        Bonito.route!(server, "/widgets/job-panel" => runtime_job_widget(client))
    end
    return server
end
