/* One renderer for owned runtime controls in Bonito, decks and the control page. */
(() => {
  "use strict";
  if (globalThis.LineCableModelsRuntimeControls) return;
  const mounted = new Map();
  const occupied = new Set(["reserving", "active", "releasing", "reconciling"]);
  const token = /^[a-z0-9][a-z0-9_-]{0,63}$/;
  const observer = new MutationObserver(() => {
    for (const [root, control] of mounted) if (!root.isConnected) control.destroy();
  });
  globalThis.addEventListener("pagehide", event => {
    // A departing iframe still has an internally connected document. Its
    // observers cannot detect removal from the parent page. BFCache instead
    // freezes the document and retains its controls for pageshow/visibility.
    if (!event.persisted) for (const control of [...mounted.values()]) control.destroy();
  });
  const node = (tag, text = "", className = "") => {
    const element = document.createElement(tag);
    element.textContent = text;
    element.className = className;
    return element;
  };
  const hint = text => node("p", text, "lc-runtime-note");
  const button = (text, action) => {
    const element = node("button", text, "lc-button lc-button-secondary");
    element.type = "button";
    element.addEventListener("click", action);
    return element;
  };
  function field(parent, title) {
    const label = node("label", "", "lc-field");
    const input = node("select", "", "lc-control-select lc-form-control");
    label.append(node("span", title, "lc-field-label"), input);
    parent.append(label);
    return input;
  }
  function options(select, rows, placeholder) {
    // Retain DOM/focus and the user's draft while polling, including an offline
    // pinned worker. Never silently replace a deliberate selection with another.
    const chosen = select.value;
    const values = [{value: "", label: placeholder}, ...rows];
    if (chosen && !values.some(row => row.value === chosen)) {
      values.push({value: chosen, label: chosen + " · no longer available", disabled: true});
    }
    const signature = JSON.stringify(values);
    if (select.dataset.options === signature) return;
    select.dataset.options = signature;
    select.replaceChildren(...values.map(row => {
      const option = node("option", row.label);
      option.value = row.value;
      option.disabled = Boolean(row.disabled);
      return option;
    }));
    select.value = chosen;
  }
  function assignment(state, role) {
    return state.assignments.filter(item => item.role === role)
      .sort((a, b) => b.generation - a.generation)[0] ?? null;
  }
  function table(parent, headings, label) {
    const wrap = node("div", "", "lc-runtime-records");
    const grid = node("table");
    grid.setAttribute("aria-label", label);
    const head = node("thead"), row = node("tr"), body = node("tbody");
    for (const title of headings) { const cell = node("th", title); cell.scope = "col"; row.append(cell); }
    head.append(row); grid.append(head, body); wrap.append(grid); parent.append(wrap);
    return body;
  }
  function record(values) {
    const row = node("tr");
    for (const value of values) row.append(node("td", value == null ? "—" : String(value)));
    return row;
  }
  function selector(root, client, config, action) {
    root.append(node("h3", "Worker · " + config.role));
    const fields = node("div", "", "lc-runtime-fields");
    const profile = field(fields, "Profile"), placement = field(fields, "Placement"), worker = field(fields, "Worker");
    options(placement, ["automatic", "pinned", "dedicated"].map(value => ({value, label:value})), "Choose placement");
    placement.value = "automatic";
    const summary = hint(""); summary.setAttribute("role", "status");
    const actions = node("div", "", "lc-runtime-actions");
    const assign = button("Assign worker", () => {
      const request = {mode:placement.value};
      if (request.mode !== "automatic") request.worker_id = worker.value || null;
      const selected = profile.value;
      void action(id => client.assign(config.role, selected, request, {requestId:id}));
    });
    const release = button("Release assignment", () => {
      const current = assignment(client.state, config.role);
      if (current) void action(id => client.release(current.id, {requestId:id}));
    });
    actions.append(assign, release); root.append(fields, summary, actions);
    let last;
    function update(state, uncertain) {
      last = [state, uncertain];
      const control = state.control;
      options(profile, (control?.profiles ?? []).filter(p => config.profiles.includes(p.id))
        .map(p => ({value:p.id, label:p.id + " · " + p.version})), "Choose profile");
      options(worker, (control?.workers ?? []).filter(w => w.registration.profiles.includes(profile.value))
        .map(w => ({value:w.registration.worker_id,
          label:w.registration.worker_id + " · " + w.registration.state + " / " + w.liveness})),
        placement.value === "pinned" ? "Choose worker (required)" : "Any eligible worker");
      const current = assignment(state, config.role);
      const holdsSlot = current && occupied.has(current.state);
      const locked = state.stale || !control?.enabled || state.pending || uncertain;
      const selectedProfile = control?.profiles.find(p => p.id === profile.value);
      const selectedWorker = control?.workers.find(w => w.registration.worker_id === worker.value);
      const explicitWorker = placement.value !== "automatic" && Boolean(worker.value);
      const workerAvailable = selectedWorker?.registration.state === "approved" && selectedWorker.liveness === "online" &&
        selectedWorker.occupied < selectedWorker.report?.capacity && selectedWorker.report?.profiles.some(p =>
          p.profile_id === selectedProfile?.id && p.version === selectedProfile.version && p.fingerprint === selectedProfile.fingerprint);
      profile.disabled = locked || Boolean(holdsSlot) || !client.runId;
      placement.disabled = profile.disabled;
      worker.disabled = profile.disabled || placement.value === "automatic";
      assign.disabled = locked || !client.runId || Boolean(holdsSlot) || !selectedProfile ||
        control?.broker !== "online" || !placement.value || (placement.value === "pinned" && !worker.value) ||
        (explicitWorker && !workerAvailable);
      release.disabled = locked || !holdsSlot || ["releasing", "reconciling"].includes(current?.state);
      summary.textContent = !client.runId ? "Open an owned application run to assign this role." :
        current ? current.worker_id + " · " + current.state + " · generation " + current.generation +
          (current.usable ? " · lease acknowledged" : " · not available for new work") :
          explicitWorker && !workerAvailable ? "Selected worker is unavailable for this profile. The pinned choice is retained; no fallback will be used." :
          "No assignment for this role.";
    }
    for (const input of [profile, placement, worker]) input.addEventListener("change", () => update(...last));
    return update;
  }
  function preparation(root, client, config, action) {
    root.append(node("h3", "Preparation · " + config.role));
    const status = hint(""); status.setAttribute("role", "status"); root.append(status);
    root.append(hint("An installed profile and an acknowledged lease are not evidence of a prepared executor."));
    const actions = node("div", "", "lc-runtime-actions");
    const prepare = button("Prepare executor", () => {
      const current = assignment(client.state, config.role);
      if (current) void action(id => client.prepare(current.id, config.parameters ?? {}, {requestId:id}));
    });
    const cancel = button("Cancel preparation", () => {
      const current = assignment(client.state, config.role), report = client.state.science[current?.id];
      if (current && report?.current_request_id) {
        const target = report.current_request_id;
        void action(id => client.cancelScientific(current.id, target, {requestId:id}));
      }
    });
    actions.append(prepare, cancel); root.append(actions);
    return (state, uncertain) => {
      const current = assignment(state, config.role);
      const terminalProfile = id => state.control?.profiles.find(profile => profile.id === id)?.kind === "terminal";
      root.hidden = current ? terminalProfile(current.profile) : Boolean(config.profiles?.length && config.profiles.every(terminalProfile));
      if (root.hidden) {
        prepare.disabled = cancel.disabled = true; root.dataset.preparation = "not-applicable"; return;
      }
      const report = state.science[current?.id];
      const usable = current?.usable && !state.stale;
      const value = usable ? report?.preparation ?? "unknown" : "unknown";
      const active = ["starting", "preparing", "executing", "closing"].includes(report?.phase);
      const locked = !usable || !state.control?.preparation_control || state.pending || uncertain;
      // A background status query is not user work and must not permanently
      // disable preparation. The agent still rejects conflicting admission.
      prepare.disabled = locked || active || report?.channel !== "online";
      cancel.disabled = locked || !["starting", "preparing"].includes(report?.phase) || !report?.current_request_id;
      status.textContent = !current ? "Not assigned" : !usable ? "Assignment is not available for preparation." :
        !report ? "Preparation unknown · waiting for executor evidence" :
        report.channel !== "online" ? "Preparation unknown · scientific channel unavailable" :
        value === "ready" ? "Ready · executor generation " + report.executor_generation + " · freshly inspected" :
        active ? report.phase + " · " + Math.round((report.progress ?? 0) * 100) + "% · " +
          (report.elapsed_seconds ?? 0).toFixed(1) + " s · " + (report.output_lines ?? 0) + " output lines" :
        report.failure ? "Preparation " + value + " · " + report.failure :
        report.accepted === false ? "Preparation " + value + " · " + report.reason :
        value === "cold" ? "Cold · explicit preparation required" : "Preparation " + value + " · waiting for fresh evidence";
      root.dataset.preparation = !current ? "unassigned" : value;
    };
  }
  function diagnostics(root) {
    root.append(node("h3", "Worker diagnostics"));
    const workers = table(root, ["Worker", "Registration", "Connection", "Occupied / capacity", "Preparation"], "Worker inventory");
    const logs = node("details", "", "lc-runtime-events");
    logs.append(node("summary", "Control events"));
    const status = hint(""), output = node("pre", "", "lc-runtime-event-log");
    output.tabIndex = 0; output.setAttribute("aria-label", "Structured control event history");
    logs.append(status, output); root.append(logs);
    let inventoryKey, eventKey;
    return state => {
      const rows = state.control?.workers ?? [];
      const key = JSON.stringify(rows);
      if (key !== inventoryKey) {
        inventoryKey = key;
        workers.replaceChildren(...(rows.length ? rows.map(w => record([w.registration.worker_id,
          w.registration.state, w.liveness, w.occupied + " / " + (w.report?.capacity ?? w.registration.capacity),
          "unknown"])) : [record(["No registered workers", "—", "—", "—", "—"])]));
      }
      const events = state.events;
      status.textContent = (state.eventsStale ? "Event connection unavailable. " : "") +
        (events?.gap || events?.localDropped ? "Incomplete history · older events were lost or evicted. " : "") +
        (events ? events.records.length + " retained events (maximum 512)." : "No event history received.");
      const currentKey = JSON.stringify([events?.epoch, events?.cursor, events?.records]);
      if (currentKey !== eventKey) {
        eventKey = currentKey;
        const atEnd = output.scrollTop + output.clientHeight >= output.scrollHeight - 4;
        output.textContent = (events?.records ?? []).map(e => [e.at, e.code, e.worker_id,
          e.run_id && "run=" + e.run_id, e.lease_id && "lease=" + e.lease_id,
          e.generation != null && "generation=" + e.generation,
          e.job_id && "job=" + e.job_id, e.executor_id && "executor=" + e.executor_id,
          e.executor_generation != null && "executor-generation=" + e.executor_generation,
          e.stage && "stage=" + e.stage].filter(Boolean).join(" · ")).join("\n");
        if (atEnd) output.scrollTop = output.scrollHeight;
      }
    };
  }
  function administration(root, client, action) {
    root.append(node("h3", "Worker registration"));
    root.append(hint("Only operator-provisioned identities can be enrolled. Approval does not install or warm an executor. Draining and disabling stop new assignments; they do not kill existing work or revoke broker credentials."));
    const fields = node("div", "", "lc-runtime-fields");
    const provisioned = field(fields, "Provisioned identity"), worker = field(fields, "Registered worker"), stateInput = field(fields, "Registration state");
    options(stateInput, ["approved", "draining", "disabled"].map(value => ({value, label:value})), "Choose state");
    const actions = node("div", "", "lc-runtime-actions");
    const enroll = button("Enroll as pending", () => {
      const selected = provisioned.value; void action(id => client.enroll(selected, {requestId:id}));
    });
    const change = button("Apply registration", () => {
      const selected = client.state.control?.workers.find(w => w.registration.worker_id === worker.value)?.registration;
      const value = stateInput.value;
      if (selected) void action(id => client.registration(selected.worker_id, value, selected.revision, {requestId:id}));
    });
    actions.append(enroll, change); root.append(fields, actions);
    let last;
    function update(state, uncertain) {
      last = [state, uncertain];
      root.hidden = !state.control?.administrator;
      const registered = state.control?.workers ?? [], trusted = state.control?.provisioned ?? [];
      options(provisioned, trusted.filter(w => !registered.some(r => r.registration.worker_id === w.worker_id))
        .map(w => ({value:w.worker_id, label:w.worker_id})), "Choose identity");
      options(worker, registered.filter(w => trusted.some(t => t.worker_id === w.registration.worker_id))
        .map(w => ({value:w.registration.worker_id, label:w.registration.worker_id + " · " + w.registration.state})), "Choose worker");
      const locked = state.stale || state.pending || uncertain || !state.control?.enabled;
      for (const input of [provisioned, worker, stateInput]) input.disabled = locked;
      enroll.disabled = locked || !provisioned.value;
      change.disabled = locked || !worker.value || !stateInput.value;
    }
    for (const input of [provisioned, worker, stateInput]) input.addEventListener("change", () => update(...last));
    return update;
  }
  function execution(root, client, config, job) {
    root.append(node("h3", "Calculation · " + config.operation));
    const status = hint(""); status.setAttribute("role","status");
    const error = hint(""); error.setAttribute("role","status");
    const provenance = hint(""); provenance.setAttribute("aria-label","Result provenance");
    const actions = node("div", "", "lc-runtime-actions");
    const invoke = method => { void job[method]().catch(() => {}); };
    const run = button("Run calculation", () => invoke("run"));
    const cancel = button("Cancel job", () => invoke("cancel"));
    const retry = button("Retry acknowledgement", () => invoke("retry"));
    const refresh = button("Refresh job", () => invoke("refresh"));
    actions.append(run,cancel,retry,refresh);
    const data = node("details"); data.append(node("summary","Result data · preview"));
    const preview = node("pre", "", "lc-runtime-event-log");
    data.append(preview); root.append(status,actions,error,provenance,data);
    let rendered;
    job.subscribe(state => {
      run.disabled = !state.canRun; cancel.disabled = !state.canCancel;
      retry.hidden = !state.canRetry && !job.pending; retry.disabled = !state.canRetry;
      refresh.disabled = !state.receipt;
      status.textContent = !client.runId ? "Open an owned application run to calculate." :
        state.phase + (state.receipt?.cancel_requested ? " · cancellation requested" : "") +
        (state.inputsPending ? " · applying input edits" : "") +
        (state.receipt?.cancel_acknowledged ? " · worker acknowledged cancellation" : "") +
        (state.superseded ? " · inputs or executor changed; this completion cannot replace the view" : "") +
        (state.awaitingEvidence ? " · waiting for matching executor evidence before updating the view" : "");
      error.textContent = state.error || ""; error.hidden = !state.error;
      root.dataset.jobPhase = state.phase; root.dataset.resultCurrent = String(state.current);
      const good = state.lastGood;
      provenance.textContent = !good ? "No successful result yet." :
        (state.current ? "Current result" : "Previous result · not current") + " · " + good.receipt.worker_id +
        " · executor generation " + good.receipt.executor_generation + " · job " + good.receipt.id +
        " · input " + good.receipt.input_hash + " · schema " + good.provenance.result.schema_version;
      data.hidden = !good;
      if (good && good !== rendered) {
        rendered = good;
        const text = JSON.stringify(good.value,null,2);
        preview.textContent = text.slice(0,4096) + (text.length > 4096 ? "\n… preview truncated; the result binding retains the complete value." : "");
      }
    });
    return job;
  }
  function validate(config) {
    if (!config || !["selector", "preparation", "diagnostics", "panel", "execution"].includes(config.kind)) throw new TypeError("Invalid runtime control kind");
    const roles = config.kind === "panel" ? config.roles ?? [] : ["selector", "preparation"].includes(config.kind) ? [config] : [];
    if (!Array.isArray(roles) || roles.length > 64 || new Set(roles.map(r => r.role)).size !== roles.length ||
        roles.some(r => !token.test(r.role) || !Array.isArray(r.profiles) || r.profiles.length > 64 ||
          r.profiles.some(p => !token.test(p)) || new Set(r.profiles).size !== r.profiles.length)) throw new TypeError("Invalid runtime role controls");
    return roles;
  }
  function mount(root, config = JSON.parse(root.dataset.lcmRuntimeControls)) {
    if (mounted.has(root)) return mounted.get(root);
    const roles = validate(config);
    if (!root.isConnected) throw new TypeError("Runtime control must be attached before mounting");
    const client = globalThis.LineCableModelsRuntimeClient.acquire(config.run_id ?? null);
    let job = null;
    try {
      if (config.kind === "execution") job = new globalThis.LineCableModelsRuntimeClient.RuntimeJob(client,config);
    } catch (error) {
      // Invalid construction must not leave an unused client or alter the DOM.
      // A client already used by another component remains owned by that view.
      if (!client.listeners.size) client.close();
      throw error;
    }
    const status = hint("Connecting to runtime control…"); status.setAttribute("role", "status");
    const message = hint(""); message.setAttribute("role", "status");
    const refresh = button("Refresh status", () => void client.refresh());
    let retryAction = null, retryId = null, destroyed = false, unsubscribe;
    const retry = button("Retry same action", () => void action(retryAction, retryId)); retry.hidden = true;
    const dismiss = button("Keep current state", () => {
      retryAction = null; retryId = null; message.textContent = "Unconfirmed action dismissed. Inspect the current assignment before making another change."; render(client.state);
    }); dismiss.hidden = true;
    const commands = node("div", "", "lc-runtime-actions"); commands.append(refresh, retry, dismiss);
    root.replaceChildren(status, commands, message);
    root.classList.add("lc-runtime-controls");
    root.dataset.runtimeKind = config.kind;
    const updates = [];
    const part = () => { const section = node("section", "", "lc-runtime-section"); root.append(section); return section; };
    async function action(callback, requestId = crypto.randomUUID()) {
      if (!callback || destroyed) return;
      retryAction = null; retryId = null; message.textContent = "Waiting for control acknowledgement…";
      try {
        await callback(requestId);
        if (!destroyed) message.textContent = "Control action acknowledged. Status reflects the latest coordinator snapshot.";
      } catch (error) {
        if (destroyed) return;
        message.textContent = error instanceof globalThis.LineCableModelsRuntimeClient.RuntimeRequestError ? error.message : "Control action could not be completed.";
        if (error.uncertain) {
          retryAction = callback; retryId = requestId;
          message.textContent += " Request " + requestId + ". Check current state; retry reuses this exact action identity.";
        }
      }
      if (!destroyed) render(client.state);
    }
    if (config.kind === "selector") updates.push(selector(part(), client, config, action));
    if (config.kind === "preparation") updates.push(preparation(part(), client, config, action));
    if (job) execution(part(),client,config,job);
    if (config.kind === "panel") {
      for (const role of roles) {
        updates.push(selector(part(), client, role, action));
        updates.push(preparation(part(), client, role, action));
      }
      updates.push(administration(part(), client, action));
    }
    if (["panel", "diagnostics"].includes(config.kind)) updates.push(diagnostics(part()));
    function render(state) {
      if (destroyed) return;
      const control = state.control;
      status.textContent = state.stale ? (state.error || "Waiting for runtime status.") + (control ? " Showing last-known values." : "") :
        !control?.enabled ? "Worker control is not configured on this publisher." : "Broker · " + control.broker;
      root.dataset.runtimeStale = String(state.stale);
      refresh.disabled = state.pending;
      retry.hidden = dismiss.hidden = !retryAction;
      message.hidden = !message.textContent;
      retry.disabled = dismiss.disabled = state.pending || state.stale;
      for (const update of updates) update(state, Boolean(retryAction));
    }
    const control = {job, destroy() {
      if (destroyed) return;
      destroyed = true; job?.close(); unsubscribe?.(); mounted.delete(root);
      if (!mounted.size) observer.disconnect();
    }};
    mounted.set(root, control);
    observer.observe(document.documentElement, {childList:true, subtree:true});
    unsubscribe = client.subscribe(render, {withEvents:["panel", "diagnostics"].includes(config.kind),
      withScience:["panel", "preparation"].includes(config.kind)});
    return control;
  }
  globalThis.LineCableModelsRuntimeControls = Object.freeze({mount});
  function start() {
    if (!globalThis.LineCableModelsRuntimeClient) return;
    for (const root of document.querySelectorAll("[data-lcm-runtime-controls]")) {
      try { mount(root); } catch { root.textContent = "Runtime controls could not be initialized."; }
    }
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start, {once:true});
  else start();
})();
