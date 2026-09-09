/* Shared same-origin runtime client. Importing this file starts no request. */
(() => {
  "use strict";
  if (globalThis.LineCableModelsRuntimeClient) return;
  const uuid = /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/;
  const token = /^[a-z0-9][a-z0-9_-]{0,63}$/;
  const operation = /^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$/;
  const digest = /^[a-f0-9]{64}$/;
  const clients = new Map();
  const MAX_RESPONSE = 4 * 1024 * 1024;
  const clock = () => globalThis.performance?.now() ?? Date.now();

  class RuntimeRequestError extends Error {
    constructor(message, {status = 0, uncertain = false, requestId = null} = {}) {
      super(message);
      this.name = "RuntimeRequestError";
      this.status = status;
      this.uncertain = uncertain;
      this.requestId = requestId;
    }
  }
  function identity(value, pattern, label) {
    if (typeof value !== "string" || !pattern.test(value)) throw new TypeError("Invalid " + label);
    return value;
  }
  async function readJSON(response) {
    const advertised = Number(response.headers.get("content-length") || 0);
    if (advertised > MAX_RESPONSE) {
      await response.body?.cancel();
      throw new Error("Response exceeds its bound");
    }
    const reader = response.body?.getReader();
    if (!reader) throw new Error("Missing response body");
    const chunks = [];
    let bytes = 0;
    try {
      for (;;) {
        const {done, value} = await reader.read();
        if (done) break;
        bytes += value.byteLength;
        if (bytes > MAX_RESPONSE) throw new Error("Response exceeds its bound");
        chunks.push(value);
      }
    } finally {
      await reader.cancel().catch(() => {});
      reader.releaseLock();
    }
    const buffer = new Uint8Array(bytes);
    let offset = 0;
    for (const chunk of chunks) { buffer.set(chunk, offset); offset += chunk.byteLength; }
    return JSON.parse(new TextDecoder("utf-8", {fatal: true}).decode(buffer));
  }
  function inventory(value) {
    if (value?.schema_version !== 1 || typeof value.enabled !== "boolean" ||
        !["disabled", "unavailable", "online", "stopped"].includes(value.broker) ||
        !Array.isArray(value.workers) || value.workers.length > 128 ||
        !Array.isArray(value.profiles) || value.profiles.length > 64 ||
        !Array.isArray(value.provisioned) || value.provisioned.length > 128) {
      throw new Error("Incompatible runtime inventory");
    }
    return value;
  }
  function assignments(value) {
    if (!Array.isArray(value) || value.length > 4096 ||
        value.some(item => !uuid.test(item?.id || "") || !uuid.test(item?.run_id || ""))) {
      throw new Error("Incompatible assignment response");
    }
    return value;
  }
  function events(value) {
    if (!uuid.test(value?.epoch || "") || !Number.isSafeInteger(value.cursor) || value.cursor < 0 ||
        typeof value.gap !== "boolean" || !Array.isArray(value.records) || value.records.length > 4096) {
      throw new Error("Incompatible event response");
    }
    return value;
  }
  function scientific(value) {
    if (!value || !["online", "offline", "stopped"].includes(value.channel) ||
        !["idle", "starting", "preparing", "executing", "failed", "closing"].includes(value.phase) ||
        !["unknown", "cold", "preparing", "ready", "failed"].includes(value.preparation) ||
        !Number.isSafeInteger(value.valid_for_ms) || value.valid_for_ms < 0 || value.valid_for_ms > 5000 ||
        typeof value.pending !== "boolean" || ![true, false, null].includes(value.accepted) ||
        !Number.isFinite(value.progress) || value.progress < 0 || value.progress > 1 ||
        !Number.isFinite(value.elapsed_seconds) || value.elapsed_seconds < 0 ||
        !Number.isSafeInteger(value.output_lines) || value.output_lines < 0 || value.output_lines > 1000000 ||
        (value.current_request_id !== null && !uuid.test(value.current_request_id || "")) ||
        (value.preparation === "ready" && (value.channel !== "online" || value.valid_for_ms === 0 ||
          !uuid.test(value.executor_id || "") || !/^[a-f0-9]{64}$/.test(value.preparation_key || "")))) {
      throw new Error("Incompatible scientific status");
    }
    return value;
  }
  function jobReceipt(value, {runId = null, jobId = null, leaseId = null, requestId = null, name = null} = {}) {
    if (!value || ["id", "request_id", "run_id", "lease_id", "worker_boot", "executor_id"]
        .some(key => !uuid.test(value[key] || "")) ||
        !token.test(value.role || "") || !token.test(value.worker_id || "") ||
        !operation.test(value.operation || "") || !digest.test(value.input_hash || "") ||
        !digest.test(value.preparation_key || "") ||
        ![value.generation, value.executor_generation].every(n => Number.isSafeInteger(n) && n > 0) ||
        !["queued", "submitted", "succeeded", "failed", "uncertain", "canceled", "revoked"].includes(value.state) ||
        !["online", "offline", "stopped"].includes(value.channel) ||
        ![value.current_assignment, value.cancel_requested, value.cancel_acknowledged].every(v => typeof v === "boolean") ||
        value.cancel_acknowledged && !value.cancel_requested ||
        !Number.isFinite(Date.parse(value.submitted_at)) || !Number.isFinite(Date.parse(value.deadline)) ||
        Date.parse(value.deadline) <= Date.parse(value.submitted_at) ||
        runId !== null && value.run_id !== runId || jobId !== null && value.id !== jobId ||
        leaseId !== null && value.lease_id !== leaseId || requestId !== null && value.request_id !== requestId ||
        name !== null && value.operation !== name) throw new Error("Incompatible job receipt");
    return value;
  }
  function jobResult(value, expected) {
    if (value?.schema_version !== 1) throw new Error("Incompatible job result");
    const receipt = jobReceipt(value.job, expected), envelope = value.result;
    if (envelope === null) return value;
    const fence = envelope?.fence, target = envelope?.execution, result = envelope?.result;
    if (envelope?.protocol_version !== "2.0" || result?.protocol_version !== "1.0" ||
        !fence || !target || result.job_id !== receipt.id || result.operation !== receipt.operation ||
        result.input_hash !== receipt.input_hash || result.worker_id !== receipt.worker_id ||
        !digest.test(fence.fingerprint || "") || result.environment_fingerprint !== fence.fingerprint ||
        target.executor_id !== receipt.executor_id || target.executor_generation !== receipt.executor_generation ||
        target.preparation_key !== receipt.preparation_key ||
        [["lease_id", "lease_id"], ["run_id", "run_id"], ["worker_id", "worker_id"],
          ["worker_boot", "worker_boot"], ["generation", "generation"], ["role", "role"]]
          .some(([key, source]) => fence[key] !== receipt[source]) ||
        typeof result.schema_version !== "string" || !result.schema_version ||
        !Array.isArray(result.warnings) || result.warnings.some(w => typeof w !== "string") ||
        [result.inline_result, result.artifact, result.failure].filter(v => v != null).length !== 1) {
      throw new Error("Job result provenance differs from its receipt");
    }
    return value;
  }
  function scientificInputs(value) {
    const encoder = new TextEncoder();
    function passive(item, depth = 0) {
      if (depth > 16) throw new TypeError("Scientific inputs exceed their nesting bound");
      if (item === null || typeof item === "boolean") return item;
      if (typeof item === "number" && Number.isFinite(item)) return item;
      if (typeof item === "string" && encoder.encode(item).byteLength <= 131072) return item;
      if (Array.isArray(item) && item.length <= 100000) return item.map(v => passive(v, depth + 1));
      if (item && Object.prototype.toString.call(item) === "[object Object]" &&
          Object.getOwnPropertySymbols(item).length === 0 && Object.keys(item).length <= 512) {
        return Object.fromEntries(Object.keys(item).sort().map(key => {
          const field = Object.getOwnPropertyDescriptor(item, key);
          if (!key || encoder.encode(key).byteLength > 256 || !field || !("value" in field)) {
            throw new TypeError("Scientific input fields must be passive");
          }
          return [key, passive(field.value, depth + 1)];
        }));
      }
      throw new TypeError("Scientific inputs must be finite passive JSON values");
    }
    if (!value || Array.isArray(value) || typeof value !== "object") throw new TypeError("Scientific inputs require an object");
    const normalized = passive(value);
    if (encoder.encode(JSON.stringify(normalized)).byteLength > 65536) throw new TypeError("Scientific inputs exceed 64 KiB");
    return normalized;
  }

  class RuntimeClient {
    constructor(runId = null, {fetcher = globalThis.fetch.bind(globalThis), pollMs = 2000, timeoutMs = 5000} = {}) {
      if (runId !== null) identity(runId, uuid, "run identity");
      if (!Number.isFinite(pollMs) || pollMs < 50 || pollMs > 60000 ||
          !Number.isFinite(timeoutMs) || timeoutMs < 50 || timeoutMs > 15000) {
        throw new TypeError("Invalid runtime polling bounds");
      }
      this.runId = runId;
      this.fetcher = fetcher;
      this.pollMs = pollMs;
      this.timeoutMs = timeoutMs;
      this.state = Object.freeze({control: null, assignments: [], events: null, science: {},
        eventsStale: false, stale: true, error: null, pending: false, updatedAt: null});
      this.listeners = new Map();
      this.controllers = new Set();
      this.refreshing = null;
      this.eventRefreshing = null;
      this.scienceRefreshing = null;
      this.scienceEpoch = 0;
      this.scienceTimer = null;
      this.timer = null;
      this.closed = false;
      this.visibility = () => {
        if (globalThis.document?.visibilityState === "hidden") this.clearScience();
        if (this.listeners.size && globalThis.document?.visibilityState !== "hidden") void this.refresh();
      };
      globalThis.document?.addEventListener("visibilitychange", this.visibility);
    }
    notify(patch) {
      this.state = Object.freeze({...this.state, ...patch});
      for (const listener of this.listeners.keys()) {
        try { listener(this.state); } catch { /* One consumer cannot block the others. */ }
      }
    }
    subscribe(listener, {withEvents = false, withScience = false} = {}) {
      if (this.closed) throw new Error("Runtime client is closed");
      if (typeof listener !== "function") throw new TypeError("A runtime listener is required");
      if (this.listeners.has(listener)) throw new TypeError("Runtime listener is already registered");
      this.listeners.set(listener, {withEvents: Boolean(withEvents), withScience: Boolean(withScience)});
      try { listener(this.state); } catch {}
      void this.refresh();
      let removed = false;
      return () => {
        if (removed) return;
        removed = true;
        this.listeners.delete(listener);
        if (!this.listeners.size) this.close();
      };
    }
    async request(method, path, body = null, timeoutMs = this.timeoutMs) {
      if (this.closed) throw new RuntimeRequestError("Runtime client is closed");
      if (!Number.isFinite(timeoutMs) || timeoutMs < 50 || timeoutMs > 15000) throw new Error("Invalid runtime request timeout");
      const controller = new AbortController();
      this.controllers.add(controller);
      const timer = setTimeout(() => controller.abort(), timeoutMs);
      const mutation = method !== "GET";
      const requestId = body?.request_id ?? null;
      try {
        const response = await this.fetcher("/runtime/api/" + path, {
          method, credentials: "same-origin", cache: "no-store", redirect: "error",
          headers: mutation ? {"Content-Type": "application/json", "X-LCM-Request": "1"} : {"Accept": "application/json"},
          body: mutation ? JSON.stringify(body) : undefined, signal: controller.signal
        });
        const data = await readJSON(response);
        if (!response.ok) {
          const reason = typeof data?.error === "string" ? data.error.slice(0, 240) : "Runtime request was rejected";
          throw new RuntimeRequestError(reason, {status: response.status,
            uncertain: mutation && response.status >= 500, requestId});
        }
        return data;
      } catch (error) {
        if (error instanceof RuntimeRequestError) throw error;
        throw new RuntimeRequestError(mutation ? "Control action was not confirmed; refresh before retrying." :
          "Runtime status is unavailable.", {uncertain: mutation, requestId});
      } finally {
        clearTimeout(timer);
        this.controllers.delete(controller);
      }
    }
    refresh() {
      if (this.closed) return Promise.resolve(this.state);
      if (this.refreshing) return this.refreshing;
      clearTimeout(this.timer);
      const promise = (async () => {
        try {
          const [control, ownedResult] = await Promise.all([
            this.request("GET", "control").then(inventory),
            this.runId ? this.request("GET", "runs/" + this.runId + "/assignments").then(assignments)
              .then(value => ({value}), error => ({error})) : {value:[]}
          ]);
          if (this.closed) return this.state;
          // A publisher without worker control still owns valid UI runs. Its
          // disabled capability must not be obscured by the absent assignment API.
          if (control.enabled && ownedResult.error) throw ownedResult.error;
          const owned = control.enabled ? ownedResult.value : [];
          this.notify({control, assignments: owned, stale: false, error: null, updatedAt: Date.now()});
          this.publishScience(); // Drop revoked/expired evidence before further I/O.
          if (control.enabled && [...this.listeners.values()].some(item => item.withEvents)) void this.refreshEvents();
          if (control.preparation_control && [...this.listeners.values()].some(item => item.withScience)) void this.refreshScience();
        } catch (error) {
          if (!this.closed) { this.notify({stale: true, error: error.message}); this.clearScience(); }
        }
        return this.state;
      })();
      this.refreshing = promise;
      void promise.finally(() => {
        this.refreshing = null;
        if (!this.closed && this.listeners.size) {
          this.timer = setTimeout(() => {
            if (globalThis.document?.visibilityState !== "hidden") void this.refresh();
          }, this.pollMs);
        }
      });
      return promise;
    }
    clearScience() {
      this.scienceEpoch++;
      clearTimeout(this.scienceTimer);
      this.notify({science: {}});
    }
    scientificAssignment(item) {
      return item.usable && this.state.control?.profiles.find(profile => profile.id === item.profile)?.kind !== "terminal";
    }
    publishScience(values = this.state.science) {
      clearTimeout(this.scienceTimer);
      const next = {}, now = clock();
      let expires = Infinity;
      for (const item of this.state.assignments) {
        if (!this.scientificAssignment(item) || this.state.stale || !this.state.control?.preparation_control) continue;
        let value = values[item.id];
        if (!value) continue;
        if (value.preparation === "ready") {
          if (value.readyUntil <= now) value = {...value, preparation: "unknown", preparation_key: null, valid_for_ms: 0};
          else expires = Math.min(expires, value.readyUntil);
        }
        next[item.id] = value;
      }
      this.notify({science: next});
      if (!this.closed && Number.isFinite(expires)) {
        this.scienceTimer = setTimeout(() => this.publishScience(), Math.max(1, Math.ceil(expires - now)));
      }
    }
    refreshScience() {
      if (this.closed || this.state.stale || this.state.pending || !this.state.control?.preparation_control ||
          globalThis.document?.visibilityState === "hidden") return Promise.resolve();
      if (this.scienceRefreshing) return this.scienceRefreshing;
      const epoch = this.scienceEpoch;
      const owned = this.state.assignments.filter(item => this.scientificAssignment(item)).slice(0, 128)[Symbol.iterator]();
      // A slow worker must not stall lease/inventory polling. Bound fan-out and
      // retain at most one current report per assignment, never scientific input.
      const promise = Promise.all(Array.from({length: 4}, async () => {
        for (const item of owned) {
          if (this.closed || epoch !== this.scienceEpoch) break;
          let value;
          const started = clock();
          try {
            value = await this.scientificStatus(item.id);
            value = {...value, readyUntil: started + value.valid_for_ms}; // subtract HTTP latency too
          } catch { value = {channel: "offline", phase: "idle", preparation: "unknown"}; }
          if (!this.closed && epoch === this.scienceEpoch) this.publishScience({...this.state.science, [item.id]: value});
        }
      }));
      this.scienceRefreshing = promise;
      void promise.finally(() => { this.scienceRefreshing = null; });
      return promise;
    }
    scientificStatus(leaseId) {
      identity(leaseId, uuid, "assignment identity");
      return this.request("GET", "assignments/" + leaseId + "/science").then(scientific);
    }
    refreshEvents() {
      if (this.closed || this.eventRefreshing) return this.eventRefreshing;
      const last = this.state.events;
      const query = last ? "?after=" + last.cursor + "&epoch=" + encodeURIComponent(last.epoch) : "";
      const promise = (async () => {
        try {
          const batch = events(await this.request("GET", "control/events" + query));
          if (this.closed) return;
          const previous = last?.epoch === batch.epoch && !batch.gap ? last.records : [];
          const joined = [...previous, ...batch.records];
          this.notify({eventsStale: false, events: {...batch, records: joined.slice(-512),
            gap: batch.gap || Boolean(last?.gap), localDropped: (last?.localDropped || 0) + Math.max(0, joined.length - 512)}});
        } catch {
          // Keep the cursor and last-good records. Reconnect resumes that exact
          // epoch; log transport delays never hold inventory polling.
          if (!this.closed) this.notify({eventsStale: true});
        }
      })();
      this.eventRefreshing = promise;
      void promise.finally(() => { this.eventRefreshing = null; });
      return promise;
    }
    async mutate(method, path, data, requestId = null) {
      if (this.state.pending) throw new RuntimeRequestError("Another control action is pending.");
      if (this.state.stale || !this.state.control?.enabled) throw new RuntimeRequestError("Refresh available runtime status first.");
      const id = identity(requestId ?? globalThis.crypto.randomUUID(), uuid, "request identity");
      this.clearScience(); // An older in-flight query cannot survive a mutation.
      this.notify({pending: true});
      try {
        // Intentionally no automatic retry: callers retain this exact request ID
        // if an explicit retry is appropriate after checking current owned state.
        return await this.request(method, path, {...data, request_id: id});
      } finally {
        // Join any pre-mutation poll, then fetch after the confirmed/uncertain
        // action. Do not mistake an older in-flight snapshot for reconciliation.
        if (this.refreshing) await this.refreshing;
        await this.refresh();
        this.notify({pending: false});
        if ([...this.listeners.values()].some(item => item.withScience)) void this.refreshScience();
      }
    }
    assign(role, profile, placement, {requestId = null} = {}) {
      if (!this.runId) throw new RuntimeRequestError("Open an owned application run before assigning a worker.");
      if (this.state.control?.broker !== "online") throw new RuntimeRequestError("Worker control is unavailable.");
      identity(role, token, "role"); identity(profile, token, "profile");
      if (!placement || !["automatic", "pinned", "dedicated"].includes(placement.mode)) {
        throw new TypeError("Invalid placement");
      }
      const selected = placement.mode === "automatic" ? {mode: "automatic"} :
        {mode: placement.mode, worker_id: placement.worker_id ?? null};
      if (selected.mode === "pinned" || selected.worker_id !== null && selected.mode !== "automatic") {
        identity(selected.worker_id, token, "worker identity");
      }
      return this.mutate("POST", "runs/" + this.runId + "/assignments", {role, profile, placement: selected}, requestId);
    }
    release(leaseId, {requestId = null} = {}) {
      identity(leaseId, uuid, "assignment identity");
      return this.mutate("DELETE", "assignments/" + leaseId, {}, requestId);
    }
    prepare(leaseId, parameters = {}, {requestId = null} = {}) {
      identity(leaseId, uuid, "assignment identity");
      if (!parameters || typeof parameters !== "object" || Array.isArray(parameters) ||
          new TextEncoder().encode(JSON.stringify(parameters)).byteLength > 65536) throw new TypeError("Invalid preparation inputs");
      return this.mutate("POST", "assignments/" + leaseId + "/science", {action: "prepare", parameters}, requestId);
    }
    cancelScientific(leaseId, targetId, {requestId = null} = {}) {
      identity(leaseId, uuid, "assignment identity"); identity(targetId, uuid, "scientific request identity");
      return this.mutate("POST", "assignments/" + leaseId + "/science", {action: "cancel", target_id: targetId}, requestId);
    }
    submitJob(leaseId, name, parameters, {requestId = null} = {}) {
      if (!this.runId) throw new RuntimeRequestError("Open an owned application run before submitting a job.");
      identity(leaseId, uuid, "assignment identity"); identity(name, operation, "registered operation");
      const inputs = scientificInputs(parameters);
      const id = identity(requestId ?? globalThis.crypto.randomUUID(), uuid, "request identity");
      return this.mutate("POST", "assignments/" + leaseId + "/jobs", {operation: name, parameters: inputs}, id)
        .then(value => {
          try { return jobReceipt(value, {runId: this.runId, leaseId, requestId: id, name}); }
          catch { throw new RuntimeRequestError("Job submission receipt could not be verified.", {uncertain: true, requestId: id}); }
        });
    }
    listJobs() {
      if (!this.runId) throw new RuntimeRequestError("An owned application run is required.");
      return this.request("GET", "runs/" + this.runId + "/jobs").then(values => {
        if (!Array.isArray(values) || values.length > 256) throw new Error("Incompatible job history");
        return values.map(value => jobReceipt(value, {runId: this.runId}));
      });
    }
    job(jobId) {
      identity(jobId, uuid, "job identity");
      return this.request("GET", "jobs/" + jobId).then(value => jobReceipt(value, {runId: this.runId, jobId}));
    }
    jobResult(jobId) {
      identity(jobId, uuid, "job identity");
      // A reference remains metadata. Never follow a worker-supplied URL, expose
      // credentials, or fetch a scientific artifact from a public hash route.
      return this.request("GET", "jobs/" + jobId + "/result")
        .then(value => jobResult(value, {runId: this.runId, jobId}));
    }
    cancelJob(jobId, {requestId = null} = {}) {
      identity(jobId, uuid, "job identity");
      const id = identity(requestId ?? globalThis.crypto.randomUUID(), uuid, "request identity");
      return this.mutate("POST", "jobs/" + jobId + "/cancel", {}, id)
        .then(value => {
          try { return jobReceipt(value, {runId: this.runId, jobId}); }
          catch { throw new RuntimeRequestError("Job cancellation receipt could not be verified.", {uncertain: true, requestId: id}); }
        });
    }
    jobArtifact(jobId) {
      identity(jobId, uuid, "job identity");
      // Two bounded object reads plus the durable receipt lookup may exceed
      // an inventory poll's timeout. This independent read never stalls polling.
      return this.request("GET", "jobs/" + jobId + "/artifact", null, 15000).then(value => {
        if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("Incompatible scientific artifact");
        return value;
      });
    }
    enroll(workerId, {requestId = null} = {}) {
      identity(workerId, token, "worker identity");
      return this.mutate("POST", "workers", {worker_id: workerId}, requestId);
    }
    registration(workerId, state, revision, {requestId = null} = {}) {
      identity(workerId, token, "worker identity");
      if (!["approved", "draining", "disabled"].includes(state) || !Number.isSafeInteger(revision) || revision < 1) {
        throw new TypeError("Invalid registration change");
      }
      return this.mutate("PATCH", "workers/" + workerId, {state, expected_revision: revision}, requestId);
    }
    close() {
      if (this.closed) return;
      this.closed = true;
      clearTimeout(this.timer);
      clearTimeout(this.scienceTimer);
      this.scienceEpoch++;
      globalThis.document?.removeEventListener("visibilitychange", this.visibility);
      for (const controller of this.controllers) controller.abort();
      this.listeners.clear();
      if (clients.get(this.runId) === this) clients.delete(this.runId);
    }
  }
  // One view's explicit job intent. Inventory and preparation stay owned by the
  // shared client; this tracker retains at most one job and one successful view.
  class RuntimeJob {
    constructor(client, {role, operation: name, parameters = {}}) {
      if (!(client instanceof RuntimeClient)) throw new TypeError("Expected a shared runtime client");
      this.client = client;
      this.role = identity(role, token, "role");
      this.operation = identity(name, operation, "registered operation");
      if (name.length > 128) throw new TypeError("Operation exceeds its bound");
      this.parameters = scientificInputs(parameters);
      this.inputKey = JSON.stringify(this.parameters);
      this.inputsValid = true;
      this.inputsPending = false;
      this.revision = 0;
      this.active = null;
      this.lastGood = null;
      this.lastGoodInvalidated = false;
      this.pending = null;
      this.phase = "idle";
      this.error = null;
      this.mutating = false;
      this.refreshing = null;
      this.timer = null;
      this.listeners = new Set();
      this.detach = null;
      this.closed = false;
      this.state = null;
      this.publish();
    }
    assignment() {
      return this.client.state.assignments.filter(a => a.role === this.role && a.run_id === this.client.runId)
        .sort((a,b) => b.generation - a.generation)[0] ?? null;
    }
    context() {
      const state = this.client.state, lease = this.assignment(), report = state.science[lease?.id];
      if (state.stale || state.control?.broker !== "online" || !lease?.usable || !uuid.test(lease.worker_boot || "") ||
          !report || report.preparation !== "ready" || report.channel !== "online" ||
          !(report.readyUntil > clock()) || !uuid.test(report.executor_id || "") ||
          !Number.isSafeInteger(report.executor_generation) || report.executor_generation < 1) return null;
      return {run_id:lease.run_id, lease_id:lease.id, role:lease.role, worker_id:lease.worker_id,
        worker_boot:lease.worker_boot, generation:lease.generation, executor_id:report.executor_id,
        executor_generation:report.executor_generation, preparation_key:report.preparation_key};
    }
    sameTarget(left, right) {
      return Boolean(left && right && ["run_id","lease_id","role","worker_id","worker_boot","generation",
        "executor_id","executor_generation","preparation_key"].every(key => left[key] === right[key]));
    }
    markInputsPending() {
      if (this.closed) return;
      this.inputsPending = true;
      this.revision++;
      if (this.active) this.active.superseded = true;
      this.publish();
    }
    setInputs(parameters) {
      if (this.closed) return;
      const wasPending = this.inputsPending;
      this.inputsPending = false;
      let next;
      try { next = scientificInputs(parameters); }
      catch (error) {
        this.inputsValid = false; this.revision++;
        if (this.active) this.active.superseded = true;
        this.error = "Scientific inputs are invalid. Correct the fields before running.";
        this.publish(); throw error;
      }
      const key = JSON.stringify(next);
      if (key === this.inputKey && this.inputsValid) {
        if (wasPending) this.publish();
        return;
      }
      this.inputsValid = true; this.error = null;
      this.parameters = next; this.inputKey = key; this.revision++;
      if (this.active) this.active.superseded = true;
      this.publish();
    }
    publish() {
      if (this.closed) return;
      const context = this.context(), lease = this.assignment(), report = this.client.state.science[lease?.id];
      const active = this.active;
      const replaced = target =>
        Boolean(lease && ["run_id","worker_id","worker_boot","generation"].some(k => lease[k] !== target[k])) ||
        Boolean(lease && lease.id !== target.lease_id) ||
        Boolean(report?.executor_id && (report.executor_id !== target.executor_id ||
          report.executor_generation !== target.executor_generation)) ||
        Boolean(report?.preparation_key && report.preparation_key !== target.preparation_key);
      if (active) {
        // A known replacement is irreversible for this request, even if the
        // user later restores the same inputs. Unknown evidence is not readiness.
        active.superseded ||= active.revision !== this.revision || replaced(active.target) ||
          active.receipt?.current_assignment === false;
        if (active.superseded) active.candidate = null;
        if (active.candidate && this.sameTarget(active.target, context) && active.receipt.current_assignment) {
          this.lastGood = Object.freeze({...active.candidate, revision:active.revision});
          this.lastGoodInvalidated = false;
          active.candidate = null; active.projected = true;
        }
      }
      if (this.lastGood && (replaced(this.lastGood.receipt) ||
          active?.receipt?.id === this.lastGood.receipt.id && active.superseded)) {
        this.lastGoodInvalidated = true;
      }
      const unfinished = Boolean(active && (!active.receipt || ["queued","submitted","uncertain"].includes(active.receipt.state)));
      const available = !this.client.state.stale && !this.client.state.pending && !this.mutating && !this.pending;
      const profile = this.client.state.control?.profiles.find(p => p.id === lease?.profile);
      const current = Boolean(this.inputsValid && !this.inputsPending && this.lastGood && !this.lastGoodInvalidated && this.lastGood.revision === this.revision &&
        this.sameTarget(this.lastGood.receipt, context) && !unfinished && !this.pending);
      this.state = Object.freeze({phase:this.phase, error:this.error, receipt:active?.receipt ?? null,
        lastGood:this.lastGood, current, superseded:Boolean(active?.superseded),
        awaitingEvidence:Boolean(active?.candidate), requestId:this.pending?.requestId ?? active?.requestId ?? null,
        inputsPending:this.inputsPending,
        canRun:Boolean(this.inputsValid && !this.inputsPending && available && !unfinished && context && this.client.runId &&
          this.client.state.control?.assigned_execution && profile?.operations.includes(this.operation)),
        canCancel:Boolean(available && active?.receipt?.current_assignment &&
          ["queued","submitted","uncertain"].includes(active.receipt.state) && !active.receipt.cancel_requested),
        canRetry:Boolean(this.pending && !this.mutating && !this.client.state.pending && !this.client.state.stale)});
      for (const listener of this.listeners) { try { listener(this.state); } catch {} }
    }
    subscribe(listener) {
      if (this.closed || typeof listener !== "function") throw new TypeError("Invalid job listener");
      this.listeners.add(listener);
      if (!this.detach) this.detach = this.client.subscribe(() => this.publish(), {withScience:true});
      try { listener(this.state); } catch {}
      return () => { this.listeners.delete(listener); if (!this.listeners.size) this.close(); };
    }
    async run() {
      if (this.closed) throw new RuntimeRequestError("This calculation control is closed.");
      this.publish();
      const target = this.context();
      if (!this.state.canRun || !target) throw new RuntimeRequestError("Select and prepare an available worker before running.");
      const active = {requestId:globalThis.crypto.randomUUID(), target,
        parameters:scientificInputs(this.parameters), revision:this.revision,
        receipt:null, candidate:null, projected:false, superseded:false};
      this.active = active;
      return this.submit(active);
    }
    async submit(active) {
      this.mutating = true; this.phase = "submitting"; this.error = null;
      this.pending = null; this.publish();
      try {
        const receipt = await this.client.submitJob(active.target.lease_id, this.operation, active.parameters,
          {requestId:active.requestId});
        if (this.closed || this.active !== active) return;
        active.receipt = receipt;
        if (!this.sameTarget(receipt,active.target)) active.superseded = true;
        this.phase = receipt.state;
      } catch (error) {
        if (this.closed) return;
        this.phase = error.uncertain ? "unconfirmed" : "rejected";
        this.error = error instanceof RuntimeRequestError ? error.message : "Job submission could not be confirmed.";
        if (error.uncertain) this.pending = {kind:"submit",requestId:active.requestId,active};
        else this.active = null;
      } finally {
        this.mutating = false; this.publish();
        if (!this.closed && this.active?.receipt) void this.refresh();
      }
    }
    async cancel(requestId = globalThis.crypto.randomUUID()) {
      this.publish();
      if (!this.state.canCancel && this.pending?.kind !== "cancel") throw new RuntimeRequestError("No cancellable job is selected.");
      const active = this.active;
      this.mutating = true; this.pending = null; this.error = null; this.publish();
      try {
        active.receipt = await this.client.cancelJob(active.receipt.id,{requestId});
        if (!this.closed) this.phase = active.receipt.state;
      } catch (error) {
        if (this.closed) return;
        this.error = "Job cancellation was not confirmed. Refresh its status before retrying.";
        if (error.uncertain) this.pending = {kind:"cancel",requestId,active};
      } finally {
        this.mutating = false; this.publish();
        if (!this.closed) void this.refresh();
      }
    }
    retry() {
      this.publish();
      if (!this.state.canRetry) return Promise.reject(new RuntimeRequestError("No unconfirmed action is ready to retry."));
      const pending = this.pending;
      return pending.kind === "submit" ? this.submit(pending.active) : this.cancel(pending.requestId);
    }
    refresh() {
      if (this.closed || !this.active?.receipt) return Promise.resolve();
      if (this.refreshing) return this.refreshing;
      clearTimeout(this.timer);
      const active = this.active;
      const promise = (async () => {
        try {
          // Durable lifecycle is independent of broker/artifact availability.
          // A revoked, unpublished job may never have a result stream at all.
          const receipt = await this.client.job(active.receipt.id);
          if (this.closed || this.active !== active) return;
          if (!this.sameTarget(receipt,active.receipt) || receipt.input_hash !== active.receipt.input_hash ||
              receipt.operation !== this.operation || receipt.request_id !== active.requestId) throw new Error("Changed job receipt");
          active.receipt = receipt; this.phase = receipt.state; this.error = null;
          if (this.pending?.kind === "cancel" && ["succeeded","failed","canceled","revoked"].includes(receipt.state)) this.pending = null;
          this.publish();
          if (!["succeeded","failed","canceled","uncertain"].includes(receipt.state)) return;
          const result = await this.client.jobResult(active.receipt.id);
          if (this.closed || this.active !== active) return;
          if (!this.sameTarget(result.job,active.receipt) || result.job.input_hash !== active.receipt.input_hash ||
              result.job.request_id !== active.requestId) throw new Error("Changed job receipt");
          active.receipt = result.job; this.phase = result.job.state; this.error = null;
          if (this.pending?.kind === "cancel" && ["succeeded","failed","canceled","revoked"].includes(result.job.state)) this.pending = null;
          this.publish();
          const envelope = result.result;
          if (result.job.state === "succeeded" && envelope && !envelope.result.failure &&
              !active.superseded && !active.projected && !active.candidate) {
            const value = envelope.result.artifact ? await this.client.jobArtifact(result.job.id) : envelope.result.inline_result;
            if (this.closed || this.active !== active || active.superseded) return;
            active.candidate = Object.freeze({receipt:result.job, provenance:envelope, value});
          } else if (envelope?.result.failure) {
            this.error = "Job failed · " + String(envelope.result.failure.category).slice(0,128);
          }
          this.publish();
        } catch {
          if (!this.closed && this.active === active) { this.error = "Job result is unavailable. Last successful data is retained."; this.publish(); }
        }
      })();
      this.refreshing = promise;
      void promise.finally(() => {
        this.refreshing = null;
        if (!this.closed && this.active !== active && this.active?.receipt) {
          void this.refresh();
        } else if (!this.closed && this.active === active && (["queued","submitted"].includes(active.receipt?.state) ||
            active.receipt?.state === "succeeded" && !active.projected && !active.superseded)) {
          this.timer = setTimeout(() => void this.refresh(), this.client.pollMs);
        }
      });
      return promise;
    }
    close() {
      if (this.closed) return;
      this.closed = true; clearTimeout(this.timer); this.detach?.(); this.detach = null;
      this.listeners.clear(); this.pending = null; this.active = null; this.lastGood = null;
      this.state = null;
    }
  }
  function acquire(runId = null) {
    if (!clients.has(runId) || clients.get(runId).closed) clients.set(runId, new RuntimeClient(runId));
    return clients.get(runId);
  }
  globalThis.LineCableModelsRuntimeClient = Object.freeze({RuntimeClient, RuntimeJob, RuntimeRequestError, acquire});
})();
