/* Private terminal transport. No broker, Bonito Observable, input log or auto-open. */
(() => {
  "use strict";
  if (globalThis.LineCableModelsTerminalClient) return;
  const UUID = /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/;
  const TOKEN = /^[a-z0-9][a-z0-9_-]{0,63}$/;
  const CHUNK = 8192, FRAME = 65536, INPUT = 32768;
  const encoder = new TextEncoder();
  const clock = () => performance.now();
  const integer = value => Number.isSafeInteger(value) && value >= 0;
  const fields = (value, names) => value && typeof value === "object" && !Array.isArray(value) &&
    Object.keys(value).length === names.length && names.every(name => Object.hasOwn(value, name));
  const reportFields = ["kind", "request_id", "accepted", "reason", "session_id", "phase", "writer_connected",
    "input_sequence", "cursor", "output_sequence", "gap", "bytes", "failure", "cleanup_pending"];
  function decode(message) {
    if (typeof message !== "string" || message.length > FRAME || encoder.encode(message).length > FRAME) throw Error("Invalid terminal frame");
    const value = JSON.parse(message);
    if (value?.kind === "hello" && fields(value, ["kind","schema_version","chunk_bytes","serialized","keepalive_seconds"]) &&
        value.schema_version === 1 && value.chunk_bytes === CHUNK && value.serialized === true && value.keepalive_seconds === 5) return value;
    if (value?.kind === "attached" && fields(value, ["kind","connection_id"]) && UUID.test(value.connection_id)) return value;
    if (value?.kind === "uncertain" && fields(value, ["kind","request_id","action","reason","automatic_retry"]) &&
        UUID.test(value.request_id) && TOKEN.test(value.action) && value.reason === "reply_unavailable" && value.automatic_retry === false) return value;
    if (value?.kind !== "report" || !fields(value, reportFields) || !UUID.test(value.request_id) ||
        typeof value.accepted !== "boolean" || typeof value.reason !== "string" || !TOKEN.test(value.reason) ||
        !(value.session_id === null || UUID.test(value.session_id)) ||
        !["unknown","starting","ready","closing","exited","failed"].includes(value.phase) ||
        typeof value.writer_connected !== "boolean" || typeof value.cleanup_pending !== "boolean" ||
        typeof value.gap !== "boolean" || !(value.failure === null || (typeof value.failure === "string" && TOKEN.test(value.failure))) ||
        ![value.input_sequence,value.cursor,value.output_sequence].every(integer) ||
        !Array.isArray(value.bytes) || value.bytes.length > CHUNK ||
        !value.bytes.every(byte => Number.isInteger(byte) && byte >= 0 && byte <= 255) ||
        value.bytes.length > value.cursor || value.cursor > value.output_sequence ||
        (!value.accepted && (value.bytes.length || value.gap)) ||
        (value.phase !== "unknown" && value.session_id === null)) throw Error("Invalid terminal report");
    return value;
  }

  class RuntimeTerminal {
    constructor(client, role, {socketFactory = url => new WebSocket(url),
      location = globalThis.location, output = async () => {}, reset = () => {}} = {}) {
      if (!client?.subscribe || typeof role !== "string" || !TOKEN.test(role) || typeof output !== "function" || typeof reset !== "function") throw TypeError("Invalid terminal context");
      if (!["https:","http:"].includes(location?.protocol)) throw TypeError("Terminal requires a same-origin HTTP context");
      this.client = client; this.role = role; this.socketFactory = socketFactory;
      this.origin = (location.protocol === "https:" ? "wss://" : "ws://") + location.host;
      this.output = output; this.reset = reset; this.listeners = new Set(); this.detach = null;
      this.target = null; this.writer = null; this.session = null; this.socket = null;
      this.pending = null; this.handshake = null; this.rendering = null; this.timer = null; this.closed = false;
      this.phase = "disconnected"; this.message = "Connect explicitly to an assigned private Julia terminal.";
      this.notice = null;
      this.uncertain = false; this.reviewRequired = false; this.busy = false; this.action = null; this.input = new Uint8Array(0);
      this.sequence = 0; this.cursor = 0; this.columns = 80; this.rows = 18;
      this.dimensions = ""; this.presenceAt = 0; this.lastRead = 0; this.cleanup = false;
      this.publish();
    }
    context() {
      const state = this.client.state;
      const lease = state.assignments.filter(item => item.run_id === this.client.runId && item.role === this.role)
        .sort((a,b) => b.generation - a.generation)[0];
      const profile = state.control?.profiles.find(item => item.id === lease?.profile);
      return globalThis.LineCableModelsRuntimeClient.runAvailability(state, this.client.runId).accepting &&
        !state.stale && state.control?.enabled && state.control.broker === "online" && lease?.usable === true &&
        UUID.test(lease.id) && UUID.test(lease.worker_boot) && integer(lease.generation) && lease.generation > 0 &&
        profile?.kind === "terminal" && profile.isolation === "container" ?
        [lease.id, lease.worker_boot, lease.generation].join(":") : null;
    }
    publish() {
      const available = Boolean(this.context()) && !this.closed;
      const connected = this.socket?.readyState === 1 && !this.handshake;
      this.state = Object.freeze({phase:this.phase, message:this.message, notice:this.notice, available, connected,
        uncertain:this.uncertain, reviewRequired:this.reviewRequired, queuedBytes:this.input.length, cleanupPending:this.cleanup,
        // Transport serialization is busy for every input/read/keepalive. Only
        // lifecycle operations are user-visible loading, not normal REPL I/O.
        activity:this.action || (["connecting","starting","closing","disconnecting"].includes(this.phase) ? this.phase : null),
        canConnect:available && !this.socket && !this.busy,
        canInput:available && connected && this.phase === "ready" && !this.uncertain && !this.reviewRequired && !this.action,
        canControl:available && connected && Boolean(this.session) && !this.uncertain && !this.action,
        canDisconnect:Boolean(this.socket), busy:this.busy || Boolean(this.action)});
      for (const listener of this.listeners) { try { listener(this.state); } catch {} }
    }
    subscribe(listener) {
      if (this.closed || typeof listener !== "function") throw TypeError("Invalid terminal listener");
      this.listeners.add(listener);
      if (!this.detach) this.detach = this.client.subscribe(() => this.inventory());
      listener(this.state);
      return () => this.listeners.delete(listener);
    }
    inventory() {
      const next = this.context();
      if (this.socket && next !== this.target) this.halt("Assignment or live inventory changed; reconnect explicitly.");
      if (next && this.target && next !== this.target) {
        this.writer = null; this.session = null; this.sequence = 0; this.cursor = 0; this.reviewRequired = false; this.reset();
      }
      if (next && !this.socket) this.target = next;
      this.publish();
    }
    send(value) {
      const message = JSON.stringify(value);
      if (!this.socket || this.socket.readyState !== 1 || this.socket.bufferedAmount > FRAME ||
          encoder.encode(message).length > FRAME) throw Error("Terminal send unavailable");
      this.socket.send(message);
    }
    receive(event, socket) {
      if (socket !== this.socket) return;
      try {
        const value = decode(event.data);
        if (this.handshake) {
          if (this.handshake.stage === "hello" && value.kind === "hello") {
            this.handshake.stage = "attached";
            this.send({action:"attach",writer_id:this.writer}); return;
          }
          if (this.handshake.stage === "attached" && value.kind === "attached") {
            const handshake = this.handshake; this.handshake = null;
            clearTimeout(handshake.timer); handshake.resolve(); return;
          }
          throw Error("Unexpected terminal handshake");
        }
        const pending = this.pending;
        if (!pending || value.request_id !== pending.packet.request_id) throw Error("Uncorrelated terminal response");
        if (value.kind === "uncertain") {
          if (value.action !== pending.packet.action) throw Error("Unexpected uncertain response");
          this.halt("Reply lost. Input was not resent; reconnect to inspect the session.", true); return;
        }
        if (value.kind !== "report") throw Error("Unexpected terminal response");
        const packet = pending.packet;
        if (value.accepted) {
          if (!UUID.test(value.session_id) ||
              (!["open","restart"].includes(packet.action) && value.session_id !== packet.session_id) ||
              (packet.action === "restart" && value.session_id === packet.session_id) ||
              (packet.action === "input" && value.input_sequence !== packet.input_sequence) ||
              (packet.action === "read" && (value.cursor < packet.after + value.bytes.length ||
                (!value.gap && value.cursor !== packet.after + value.bytes.length))) ||
              (packet.action !== "read" && (value.cursor || value.bytes.length || value.gap))) throw Error("Terminal response changed its stream");
        }
        clearTimeout(pending.timer); this.pending = null; pending.resolve(value);
      } catch { this.halt("Terminal protocol failed; connection closed without replaying input."); }
    }
    request(action, values = {}) {
      if (this.pending || this.handshake || this.context() !== this.target) return Promise.reject(Error("Terminal is not available"));
      const packet = {action,request_id:crypto.randomUUID(),session_id:action === "open" ? null : this.session,
        input_sequence:0,after:0,columns:0,rows:0,bytes:[],retry:false,...values};
      return new Promise((resolve,reject) => {
        const timer = setTimeout(() => this.halt("Reply deadline expired. Input was not resent; reconnect to inspect.", true),
          action === "restart" ? 32000 : 7000);
        this.pending = {packet,resolve,reject,timer};
        try { this.send(packet); } catch { this.halt("Terminal connection was lost; input was not resent."); }
      });
    }
    async connect() {
      if (!this.state.canConnect) return false;
      this.target = this.context(); this.writer ??= crypto.randomUUID();
      this.busy = true; this.uncertain = false; this.phase = "connecting";
      this.notice = null; this.message = "Connecting to the private terminal…"; this.publish();
      let socket;
      try {
        socket = this.socketFactory(this.origin + "/runtime/api/assignments/" + this.target.split(":")[0] + "/terminal");
        this.socket = socket;
        await new Promise((resolve,reject) => {
          const timer = setTimeout(() => this.halt("Private terminal attachment timed out."), 7000);
          this.handshake = {stage:"hello",resolve,reject,timer};
          socket.addEventListener("message", event => this.receive(event, socket));
          const lost = () => { if (socket === this.socket) this.halt("Terminal disconnected. Reconnect explicitly; queued input was discarded."); };
          socket.addEventListener("close", lost); socket.addEventListener("error", lost);
        });
        const report = await this.request("open", {columns:this.columns,rows:this.rows});
        if (socket !== this.socket) return false;
        if (!this.accept(report)) return false;
        if (this.session !== report.session_id) { this.cursor = 0; this.reset(); }
        this.session = report.session_id; this.sequence = report.input_sequence;
        this.presenceAt = clock(); this.dimensions = "";
        return true;
      } catch { if (!socket || socket === this.socket) this.halt("Private terminal connection could not be established."); return false; }
      finally { this.busy = false; this.publish(); this.schedule(); }
    }
    accept(report) {
      if (!report.accepted) {
        this.halt("Terminal action rejected: " + report.reason + ". No input was replayed."); return false;
      }
      this.phase = this.action === "disconnect" ? "disconnecting" : report.phase; this.cleanup = report.cleanup_pending;
      this.message = this.action === "disconnect" ? "Disconnecting the private writer…" :
        report.failure ? "Terminal " + report.phase + " · " + report.failure :
        report.phase === "ready" ? "Julia REPL running · input acknowledgements confirm queuing, not evaluation." :
        report.phase === "starting" ? "Starting Julia… input is disabled until the REPL is ready." :
        "Terminal " + report.phase + (report.cleanup_pending ? " · cleanup pending" : "");
      return true;
    }
    enqueue(value) {
      if (!this.state.canInput) return false;
      // Bound before encoding/copying; reject a whole oversized paste, never half a command.
      if (!(typeof value === "string" || value instanceof Uint8Array) || value.length > INPUT) {
        this.notice = "Input exceeds the 32 KiB queue; nothing from that paste was queued."; this.publish(); return false;
      }
      const bytes = typeof value === "string" ? encoder.encode(value) : value;
      if (bytes.length + this.input.length > INPUT) {
        this.notice = "Input queue is full; additional input was not queued."; this.publish(); return false;
      }
      const next = new Uint8Array(this.input.length + bytes.length);
      next.set(this.input); next.set(bytes,this.input.length); this.input = next;
      this.notice = null;
      this.publish(); this.schedule(0); return true;
    }
    interrupt() {
      if (!this.state.canInput) return false;
      this.input = new Uint8Array(0); return this.enqueue(new Uint8Array([3]));
    }
    resumeInput() {
      if (!this.state.canControl || this.phase !== "ready" || !this.reviewRequired) return false;
      this.reviewRequired = false; this.publish(); return true;
    }
    resize(columns, rows) {
      if (![columns,rows].every(value => Number.isInteger(value) && value >= 1 && value <= 1000)) return false;
      this.columns = columns; this.rows = rows; this.schedule(0); return true;
    }
    control(action) {
      if (!["stop","restart"].includes(action) || !this.state.canControl) return false;
      this.input = new Uint8Array(0); this.action = action; this.publish(); this.schedule(0); return true;
    }
    disconnect() {
      this.input = new Uint8Array(0);
      if (!this.socket || !this.session || this.handshake || this.uncertain) {
        this.halt("Disconnected; unsent input was discarded."); return;
      }
      // Serialize intentional detach after the current request. Do not offer a
      // new writer attachment while the old one still owns the remote session.
      this.action = "disconnect"; this.phase = "disconnecting";
      this.message = "Disconnecting the private writer…"; this.publish(); this.schedule(0);
    }
    schedule(delay = 100) {
      clearTimeout(this.timer); this.timer = null;
      if (!this.closed && this.socket && !this.busy && this.session && !this.uncertain)
        this.timer = setTimeout(() => { this.timer = null; void this.pump(); }, delay);
    }
    async pump() {
      if (this.busy || this.closed || !this.socket || this.context() !== this.target) return;
      this.busy = true; const socket = this.socket;
      try {
        let action, values = {};
        if (this.action) {
          action = this.action;
          if (action === "restart") values = {columns:this.columns,rows:this.rows};
        } else if (["starting","ready"].includes(this.phase) && clock() - this.presenceAt >= 5000) {
          action = "keepalive";
        } else if (this.phase === "ready" && clock() - this.lastRead >= 250) {
          action = "read"; values = {after:this.cursor};
        } else if (this.phase === "ready" && this.input.length) {
          action = "input"; const bytes = this.input.slice(0,CHUNK); this.input = this.input.slice(bytes.length);
          values = {input_sequence:this.sequence + 1,bytes:Array.from(bytes)};
        } else if (this.phase === "ready" && this.dimensions !== this.columns + ":" + this.rows) {
          action = "resize"; values = {columns:this.columns,rows:this.rows};
        } else if (this.phase === "ready") { action = "read"; values = {after:this.cursor}; }
        else if (clock() - this.lastRead >= 250) action = "status";
        else return;
        const report = await this.request(action, values);
        if (socket !== this.socket || !this.accept(report)) return;
        if (action === "disconnect") {
          if (report.writer_connected) throw Error("Terminal writer was not released");
          this.halt("Disconnected. Reconnect within the worker's grace period; unsent input was discarded.");
          return;
        }
        if (["read","status"].includes(action)) this.lastRead = clock();
        if (["keepalive","input","resize","restart"].includes(action)) this.presenceAt = clock();
        if (["input","restart"].includes(action)) this.sequence = report.input_sequence;
        if (action === "resize") this.dimensions = values.columns + ":" + values.rows;
        if (action === "restart") {
          this.session = report.session_id; this.cursor = 0; this.dimensions = ""; this.reviewRequired = false; this.notice = null; this.reset();
        }
        if (action === "read") {
          if (report.gap) { this.reset(); this.notice = "Earlier output was evicted; showing the retained terminal tail."; }
          // The renderer's completion callback releases backpressure. At most one
          // chunk is decoded/rendered; no output poll while a write is outstanding.
          await this.renderOutput(Uint8Array.from(report.bytes), report.gap);
          if (socket !== this.socket) return;
          this.cursor = report.cursor;
        }
        if (this.action === action) this.action = null;
      } catch { if (socket === this.socket) this.halt("Terminal processing failed; input was not resent."); }
      finally { this.busy = false; this.publish(); this.schedule(); }
    }
    renderOutput(bytes, gap) {
      return new Promise((resolve,reject) => {
        const pending = {reject,timer:null}; this.rendering = pending;
        const complete = error => {
          if (this.rendering !== pending) return;
          this.rendering = null; clearTimeout(pending.timer); error ? reject(error) : resolve();
        };
        pending.timer = setTimeout(() => complete(Error("Terminal renderer stalled")), 3000);
        Promise.resolve().then(() => this.output(bytes,gap)).then(() => complete(), () => complete(Error("Terminal renderer failed")));
      });
    }
    halt(message, uncertain = false) {
      const rendering = this.rendering; this.rendering = null;
      if (rendering) { clearTimeout(rendering.timer); rendering.reject(Error("Terminal view disconnected")); }
      const pending = this.pending;
      this.uncertain ||= uncertain || Boolean(pending && !["open","status","read"].includes(pending.packet.action));
      this.reviewRequired ||= this.uncertain;
      this.pending = null;
      if (pending) { clearTimeout(pending.timer); pending.reject(Error("Terminal reply unavailable")); }
      const handshake = this.handshake; this.handshake = null;
      if (handshake) { clearTimeout(handshake.timer); handshake.reject(Error("Terminal attachment unavailable")); }
      const socket = this.socket; this.socket = null;
      try { socket?.close(1000); } catch {}
      clearTimeout(this.timer); this.timer = null; this.input = new Uint8Array(0); this.action = null;
      this.phase = this.uncertain ? "uncertain" : "disconnected";
      this.message = this.uncertain ? message + " Execution may already have occurred; do not retype blindly." : message;
      this.publish();
    }
    destroy() {
      if (this.closed) return;
      this.closed = true; this.halt("Terminal view closed."); this.detach?.(); this.detach = null;
      this.listeners.clear(); this.session = null; this.writer = null; this.output = async () => {}; this.reset = () => {};
    }
  }
  globalThis.LineCableModelsTerminalClient = Object.freeze({RuntimeTerminal});
})();
