/* One locally served xterm renderer for gallery, workbench and presentation use. */
(() => {
  "use strict";
  if (globalThis.LineCableModelsTerminal) return;
  const mounted = new Map();
  const node = (tag, text = "", className = "") => {
    const result = document.createElement(tag); result.textContent = text; result.className = className; return result;
  };
  const setText = (element, value) => { if (element.textContent !== value) element.textContent = value; };
  const observer = new MutationObserver(() => {
    for (const [root, control] of mounted) if (!root.isConnected) control.destroy();
  });
  globalThis.addEventListener("pagehide", event => {
    // A frozen page cannot maintain private writer presence. Reconnect is explicit
    // even after BFCache return; a departing iframe also closes its attachment.
    for (const control of [...mounted.values()]) event.persisted ? control.disconnect() : control.destroy();
  });
  function theme(root) {
    const style = getComputedStyle(root), token = name => style.getPropertyValue("--lc-" + name).trim();
    return {background:token("console-bg"),foreground:token("text"),cursor:token("focus"),
      cursorAccent:token("console-bg"),selectionBackground:token("option-selected-bg"),selectionForeground:token("heading"),
      black:token("muted"),brightBlack:token("muted"),white:token("text"),brightWhite:token("heading"),
      red:token("danger"),brightRed:token("danger"),green:token("focus"),brightGreen:token("focus"),
      yellow:token("warning"),brightYellow:token("warning"),blue:token("link"),brightBlue:token("link"),
      magenta:token("code-preprocessor"),brightMagenta:token("code-preprocessor"),cyan:token("code-name"),brightCyan:token("code-name")};
  }
  function mount(root) {
    if (mounted.has(root)) return mounted.get(root);
    const config = JSON.parse(root.dataset.lcmRuntimeTerminal);
    const uuid = /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/;
    if (!config || typeof config.role !== "string" || !/^[a-z0-9][a-z0-9_-]{0,63}$/.test(config.role) ||
        !(config.run_id === null || (typeof config.run_id === "string" && uuid.test(config.run_id))) ||
        !(config.kind === undefined || config.kind === "terminal") ||
        typeof config.title !== "string" || !config.title.length || /[\x00-\x1f\x7f-\x9f]/.test(config.title) || Array.from(config.title).length > 120 ||
        !Number.isInteger(config.rows) || config.rows < 6 || config.rows > 60) throw TypeError("Invalid terminal configuration");
    const {Terminal,FitAddon} = globalThis.LineCableModelsTerminalVendor;
    // Validate before acquiring a client or altering the host element.
    const disposables = [], writes = new Set(); let destroyed = false, fitTimer, themeTimer;
    root.classList.add("lc-runtime-controls", "lc-runtime-terminal");
    const head = node("div", "", "lc-terminal-heading"), heading = node("h3", config.title);
    const badge = node("span", "Disconnected", "lc-terminal-phase lc-status-indicator lc-activity-status lc-activity-stable");
    head.append(heading, badge);
    const actions = node("div", "", "lc-runtime-actions");
    const status = node("p", "", "lc-runtime-note lc-terminal-status lc-status-indicator"); status.setAttribute("role", "status");
    const recovery = node("a", "", "lc-terminal-recovery"); recovery.target = "_top"; recovery.hidden = true;
    const viewport = node("div", "", "lc-terminal-viewport");
    viewport.style.setProperty("--lcm-terminal-rows", String(config.rows));
    viewport.setAttribute("aria-label", config.title); viewport.setAttribute("role", "group");
    // FitAddon reads the parent's computed dimensions. Give it the actual
    // content box, not the surrounding padded/bordered frame.
    const screen = node("div", "", "lc-terminal-screen"); viewport.append(screen);
    const confirmation = node("div", "", "lc-terminal-confirm"); confirmation.hidden = true;
    const warning = node("span", ""), explanation = node("p",
      "Private, disposable Julia session. Stop or restart discards its memory. This page persists no input history.", "lc-runtime-note");
    root.replaceChildren(head,actions,status,recovery,confirmation,viewport,explanation);
    const terminal = new Terminal({cols:80,rows:config.rows,fontSize:13,lineHeight:1.35,
      fontFamily:'"JuliaMono", "SFMono-Regular", Consolas, monospace',scrollback:1000,
      cursorBlink:true,cursorInactiveStyle:"none",disableStdin:true,logLevel:"off",
      minimumContrastRatio:4.5,allowProposedApi:false,windowOptions:{},
      linkHandler:{activate() {},hover() {},leave() {}},theme:theme(root)});
    const fit = new FitAddon(); terminal.loadAddon(fit);
    // No output-controlled page title, hyperlinks, clipboard, theme, or window
    // operation. The parser still handles ordinary Julia ANSI/VT terminal output.
    for (const code of [0,1,2,4,8,10,11,12,52,104,110,111,112])
      disposables.push(terminal.parser.registerOscHandler(code, () => true));
    terminal.open(screen);
    const client = globalThis.LineCableModelsRuntimeClient.acquire(config.run_id);
    const write = bytes => new Promise((resolve,reject) => {
      if (destroyed) return reject(Error("Terminal view closed"));
      const pending = {reject,timer:null}; writes.add(pending);
      const complete = error => {
        if (!writes.delete(pending)) return;
        clearTimeout(pending.timer); error ? reject(error) : resolve();
      };
      pending.timer = setTimeout(() => complete(Error("Terminal renderer stalled")), 3000);
      try { terminal.write(bytes, () => complete()); } catch { complete(Error("Terminal renderer unavailable")); }
    });
    const transport = new globalThis.LineCableModelsTerminalClient.RuntimeTerminal(client, config.role,
      {output:write,reset:() => terminal.reset()});
    const button = (label, action) => {
      const element = node("button",label,"lc-button lc-button-secondary"); element.type = "button";
      element.addEventListener("click", () => {
        // Deliberately exclude keystrokes, submitted Julia and terminal output.
        client.recordActivity("terminal_action", "Terminal · " + label + " requested.");
        action();
      }); return element;
    };
    const connect = button("Connect", () => { void transport.connect().then(ok => {
      if (ok && !destroyed && document.activeElement === connect) terminal.focus();
    }); });
    const disconnect = button("Disconnect", () => transport.disconnect());
    const interrupt = button("Interrupt", () => { if (transport.interrupt()) terminal.focus(); });
    interrupt.title = "Send Ctrl-C; discard input that has not yet been sent";
    let confirmedAction = null;
    const confirm = button("Confirm", () => {
      const action = confirmedAction; confirmation.hidden = true; confirmedAction = null;
      if (action) transport.control(action);
    });
    const cancel = button("Cancel", () => { confirmation.hidden = true; confirmedAction = null; });
    confirmation.append(warning,confirm,cancel);
    const ask = action => {
      confirmedAction = action; warning.textContent = action === "restart" ?
        "Restart creates a clean Julia namespace. Unsaved variables are lost." : "Stop discards this Julia session and its variables.";
      confirm.textContent = action === "restart" ? "Restart Julia" : "Stop Julia";
      confirmation.hidden = false; cancel.focus();
    };
    const stop = button("Stop", () => ask("stop")), restart = button("Restart", () => ask("restart"));
    const clear = button("Clear view", () => terminal.clear());
    const resume = button("Resume input", () => { if (transport.resumeInput()) terminal.focus(); });
    const refresh = button("Refresh status", () => void client.refreshStatus());
    clear.title = "Clear visible scrollback; does not reset Julia or its variables";
    actions.append(connect,disconnect,interrupt,stop,restart,clear,resume,refresh);
    let lastPhase;
    const unsubscribe = transport.subscribe(state => {
      if (root.dataset.terminalPhase !== state.phase) root.dataset.terminalPhase = state.phase;
      setText(badge,state.phase);
      const tone = globalThis.LineCableModelsRuntimeClient.statusTone;
      badge.dataset.tone = tone(state.uncertain ? "uncertain" : state.phase);
      badge.dataset.busy = String(Boolean(state.activity));
      if (lastPhase !== state.phase) {
        lastPhase = state.phase;
        // The shared journal notifies all views. Store the phase before doing
        // so, avoiding recursive recording through the inventory subscription.
        client.recordActivity("terminal_phase", "Terminal · " + state.phase, {tone:tone(state.phase)});
      }
      const availability = globalThis.LineCableModelsRuntimeClient.runAvailability(client.state, client.runId);
      const setup = !availability.accepting ? availability.message : client.state.stale ? "Runtime status is unknown. Refresh status before connecting." :
        client.state.control?.broker !== "online" ? "Broker unavailable. Restore the connection, then refresh status; worker readiness cannot be checked." :
        "Terminal not assigned to this run. Open Workers and preparation, assign the julia-terminal profile, then return here and select Connect.";
      // The initial setup hint must not hide why an attempted connection lost
      // authority. Retain the transport's explicit reason until user recovery.
      setText(status, !transport.writer && !state.available && !state.connected && !state.uncertain ?
        setup :
        state.reviewRequired && state.connected ? "A previous action was unconfirmed. Inspect output before choosing Resume input, or restart Julia for a clean session." :
        state.connected && state.notice ? state.notice : state.message);
      status.dataset.tone = !availability.accepting ? availability.tone : state.uncertain ? "warning" : !state.available ? "warning" : tone(state.phase);
      recovery.hidden = availability.accepting; recovery.href = availability.href; setText(recovery, availability.label);
      connect.disabled = !state.canConnect; setText(connect,transport.session ? "Reconnect" : "Connect");
      connect.dataset.busy = String(state.busy && !state.connected);
      connect.setAttribute("aria-busy", connect.dataset.busy);
      if (state.busy && !state.connected) setText(connect, "Connecting…");
      refresh.disabled = client.state.refreshing || client.state.pending;
      refresh.dataset.busy = String(client.state.refreshing); refresh.setAttribute("aria-busy", refresh.dataset.busy);
      setText(refresh, client.state.refreshing ? "Refreshing…" : "Refresh status");
      disconnect.disabled = !state.canDisconnect; interrupt.disabled = !state.canInput;
      stop.disabled = !state.canControl || ["closing","exited","failed"].includes(state.phase);
      restart.disabled = !state.canControl || state.cleanupPending;
      resume.hidden = !state.reviewRequired || !state.connected;
      resume.disabled = !state.canControl || state.phase !== "ready";
      confirm.disabled = !state.canControl;
      if (terminal.options.disableStdin !== !state.canInput) terminal.options.disableStdin = !state.canInput;
      if (terminal.options.cursorBlink !== state.canInput) terminal.options.cursorBlink = state.canInput;
      if (!state.canControl) { confirmation.hidden = true; confirmedAction = null; }
    });
    disposables.push(terminal.onData(data => transport.enqueue(data)));
    disposables.push(terminal.onBinary(data => {
      if (data.length <= 32768) transport.enqueue(Uint8Array.from(data, char => char.charCodeAt(0) & 255));
    }));
    // Let xterm and native controls handle their keys, then stop deck navigation.
    // Capture-phase deck handlers already recognize editable targets.
    const stopKeys = event => event.stopPropagation();
    root.addEventListener("keydown",stopKeys); root.addEventListener("keyup",stopKeys);
    const resize = () => {
      clearTimeout(fitTimer);
      fitTimer = setTimeout(() => {
        if (destroyed || !root.isConnected || viewport.clientWidth < 40 || viewport.clientHeight < 40) return;
        try { fit.fit(); transport.resize(terminal.cols,terminal.rows); } catch {}
      }, 100);
    };
    const resizeObserver = new ResizeObserver(resize); resizeObserver.observe(viewport);
    const applyTheme = () => {
      clearTimeout(themeTimer);
      themeTimer = setTimeout(() => { if (!destroyed) { terminal.options.theme = theme(root); resize(); } }, 0);
    };
    globalThis.addEventListener("lcm:theme-changed",applyTheme);
    document.fonts?.ready.then(() => { if (!destroyed) resize(); });
    const control = {disconnect:() => transport.disconnect(),destroy() {
      if (destroyed) return; destroyed = true;
      unsubscribe(); transport.destroy(); resizeObserver.disconnect();
      clearTimeout(fitTimer); clearTimeout(themeTimer); globalThis.removeEventListener("lcm:theme-changed",applyTheme);
      root.removeEventListener("keydown",stopKeys); root.removeEventListener("keyup",stopKeys);
      for (const pending of writes) { clearTimeout(pending.timer); pending.reject(Error("Terminal view closed")); }
      writes.clear(); for (const disposable of disposables) disposable.dispose(); terminal.dispose();
      mounted.delete(root); if (!mounted.size) observer.disconnect(); root.replaceChildren();
    }};
    mounted.set(root,control); observer.observe(document.documentElement,{childList:true,subtree:true}); resize();
    return control;
  }
  globalThis.LineCableModelsTerminal = Object.freeze({mount});
})();
