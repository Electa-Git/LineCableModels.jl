/* Explicit browser-only transport fixture. No real runtime or REPL is claimed. */
(() => {
  const fixture = {sent:[],sockets:[],renderers:[],uncertain:false,phase:"starting",output:[],inputSequence:0,session:null,keys:0};
  const encoder = new TextEncoder();
  fixture.append = text => fixture.output.push(...encoder.encode(text));
  class Socket extends EventTarget {
    readyState = 1; bufferedAmount = 0;
    constructor(url) { super(); fixture.sockets.push(this); this.url = url;
      queueMicrotask(() => this.emit({kind:"hello",schema_version:1,chunk_bytes:8192,serialized:true,keepalive_seconds:5})); }
    emit(value) { this.dispatchEvent(new MessageEvent("message",{data:JSON.stringify(value)})); }
    close() { this.readyState = 3; this.dispatchEvent(new Event("close")); }
    send(raw) {
      const packet = JSON.parse(raw); fixture.sent.push(packet);
      setTimeout(() => {
        if (this.readyState !== 1) return;
        if (packet.action === "attach") return this.emit({kind:"attached",connection_id:crypto.randomUUID()});
        if (packet.action === "input" && fixture.uncertain) {
          fixture.uncertain = false;
          return this.emit({kind:"uncertain",request_id:packet.request_id,action:"input",reason:"reply_unavailable",automatic_retry:false});
        }
        if (packet.action === "open") fixture.session ??= crypto.randomUUID();
        if (packet.action === "restart") {
          fixture.session = crypto.randomUUID(); fixture.output = []; fixture.inputSequence = 0;
          fixture.phase = "starting";
        }
        if (packet.action === "stop") fixture.phase = "exited";
        if (packet.action === "input") {
          fixture.inputSequence = packet.input_sequence; fixture.output.push(...packet.bytes);
        }
        const bytes = packet.action === "read" ? fixture.output.slice(packet.after,packet.after+8192) : [];
        this.emit({kind:"report",request_id:packet.request_id,accepted:true,reason:"accepted",session_id:fixture.session,
          phase:fixture.phase,writer_connected:packet.action!=="disconnect",input_sequence:fixture.inputSequence,
          cursor:packet.action === "read" ? packet.after+bytes.length : 0,output_sequence:fixture.output.length,
          gap:false,bytes,failure:null,cleanup_pending:false});
      }, fixture.replyDelay || 0);
    }
  }
  const Vendor = LineCableModelsTerminalVendor.Terminal;
  globalThis.LineCableModelsTerminalVendor = {...LineCableModelsTerminalVendor, Terminal:class extends Vendor {
    constructor(...args) {super(...args);fixture.renderers.push(this);}
  }};
  const Transport = LineCableModelsTerminalClient.RuntimeTerminal;
  globalThis.LineCableModelsTerminalClient = {RuntimeTerminal:class extends Transport {
    constructor(client,role,options) {super(client,role,{...options,socketFactory:url=>new Socket(url)});fixture.transport=this;}
  }};
  fixture.mount = () => LineCableModelsTerminal.mount(document.querySelector("#terminal"));
  fixture.mount(); globalThis.__terminalFixture = fixture;
  fixture.button = label => [...document.querySelectorAll("#terminal button")].find(button => button.textContent === label);
  window.addEventListener("keydown", () => fixture.keys++);
})();
