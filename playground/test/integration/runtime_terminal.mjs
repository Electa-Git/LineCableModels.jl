#!/usr/bin/env node
// Serialized browser transport contracts; explicit in-memory socket fixture, not a REPL.
import assert from "node:assert/strict";
import {readFile} from "node:fs/promises";
import {randomUUID,webcrypto} from "node:crypto";
import {setTimeout as delay} from "node:timers/promises";
globalThis.crypto ??= webcrypto;
await import(new URL("../../assets/runtime-client.js",import.meta.url));
await import(new URL("../../assets/runtime-terminal-client.js",import.meta.url));
const {RuntimeTerminal} = globalThis.LineCableModelsTerminalClient;
const until = async predicate => {
  for (let i=0;i<300;i++) { if (predicate()) return; await delay(10); }
  throw Error("Terminal fixture deadline exceeded");
};
let checks = 0;
function fixture(options = {}) {
  const run = randomUUID(), lease = randomUUID();
  const listeners = new Set(), sent = [], sockets = [], rendered = [];
  let session = null, inputSequence = 0, tail = [], writer = null, resets = 0;
  let hook = null, connected = false;
  const client = {runId:run,state:{stale:false,run:{id:run,state:"running"},runStale:false,control:{enabled:true,broker:"online",profiles:[{id:"repl",kind:"terminal",isolation:"container"}]},
    assignments:[{id:lease,run_id:run,role:"terminal",profile:"repl",worker_boot:randomUUID(),generation:1,usable:true}]},
    subscribe(fn) { listeners.add(fn); fn(this.state); return () => listeners.delete(fn); }};
  class Socket extends EventTarget {
    readyState = 1; bufferedAmount = 0;
    constructor(url) { super(); this.url = url; sockets.push(this); queueMicrotask(() => this.emit({kind:"hello",schema_version:1,chunk_bytes:8192,serialized:true,keepalive_seconds:5})); }
    emit(value) { this.dispatchEvent(new MessageEvent("message",{data:JSON.stringify(value)})); }
    close() { this.readyState = 3; connected = false; this.dispatchEvent(new Event("close")); }
    send(raw) {
      const packet = JSON.parse(raw); sent.push(packet);
      queueMicrotask(() => {
        if (this.readyState !== 1) return;
        if (packet.action === "attach") {
          writer = packet.writer_id; this.emit({kind:"attached",connection_id:randomUUID()}); return;
        }
        if (hook?.(packet,this)) return;
        let cursor = 0, bytes = [], gap = false;
        if (packet.action === "open") { session ??= randomUUID(); connected = true; }
        if (packet.action === "disconnect") connected = false;
        if (packet.action === "restart") { session = randomUUID(); inputSequence = 0; tail = []; }
        if (packet.action === "input") { inputSequence = packet.input_sequence; tail.push(...packet.bytes); }
        if (packet.action === "read") { bytes = tail.slice(packet.after,packet.after+8192); cursor = packet.after+bytes.length; }
        this.emit({kind:"report",request_id:packet.request_id,accepted:true,reason:"accepted",session_id:session,
          phase:"ready",writer_connected:connected,input_sequence:inputSequence,cursor,output_sequence:tail.length,
          gap,bytes,failure:null,cleanup_pending:false});
      });
    }
  }
  const terminal = new RuntimeTerminal(client,"terminal",{location:{protocol:"https:",host:"owned.example"},
    socketFactory:url => new Socket(url),output:async bytes => {rendered.push(...bytes); if (options.output) await options.output(bytes);},
    reset:() => {resets++;rendered.length=0;}});
  const unsubscribe = terminal.subscribe(() => {});
  return {client,terminal,sent,sockets,rendered,listeners,get writer(){return writer;},get resets(){return resets;},
    hook:fn => {hook = fn;},notify:() => {for (const fn of listeners) fn(client.state);},
    close:() => {unsubscribe(); terminal.destroy();}};
}

assert.throws(() => new RuntimeTerminal({subscribe(){}},undefined,{location:{protocol:"https:",host:"x"}})); checks++;
const f = fixture();
try {
  assert.equal(f.sockets.length,0); assert.equal(f.terminal.state.canConnect,true); checks+=2;
  assert.equal(await f.terminal.connect(),true); checks++;
  assert.match(f.sockets[0].url,/^wss:\/\/owned.example\/runtime\/api\/assignments\/[a-f0-9-]+\/terminal$/); checks++;
  assert.equal(f.terminal.state.canInput,true); checks++;
  assert.equal(f.terminal.enqueue("α = 3\r"),true);
  await until(() => new TextDecoder().decode(Uint8Array.from(f.rendered)) === "α = 3\r"); checks++;
  const first = f.sent.find(p => p.action === "input");
  assert.equal(first.input_sequence,1); assert.equal(first.retry,false);
  assert.deepEqual(Object.keys(first).sort(),["action","request_id","session_id","input_sequence","after","columns","rows","bytes","retry"].sort()); checks+=3;
  assert.equal(f.terminal.enqueue("x".repeat(32769)),false);
  assert.equal(f.terminal.enqueue("α".repeat(32768)),false); checks+=2;
  assert.equal(f.terminal.resize(124,27),true); assert.equal(f.terminal.resize(1001,10),false);
  await until(() => f.sent.some(p => p.action === "resize" && p.columns === 124 && p.rows === 27)); checks+=3;
  const session = f.terminal.session, writer = f.writer;
  f.terminal.disconnect();
  assert.equal(f.terminal.state.canConnect,false); assert.equal(f.terminal.state.canInput,false); checks+=2;
  await until(() => f.terminal.state.canConnect);
  assert.equal(f.sent.filter(p=>p.action==="disconnect").length,1); checks++;
  assert.equal(await f.terminal.connect(),true); assert.equal(f.terminal.session,session); assert.equal(f.writer,writer); checks+=3;
  assert.equal(f.terminal.interrupt(),true);
  await until(() => f.sent.some(p => p.action === "input" && p.bytes.length === 1 && p.bytes[0] === 3)); checks+=2;
  assert.equal(f.terminal.control("restart"),true);
  await until(() => f.terminal.session !== session);
  assert.equal(f.terminal.sequence,0); assert.equal(f.resets,2); checks+=3;
  f.hook((packet,socket) => {
    if (packet.action !== "input") return false;
    socket.emit({kind:"uncertain",request_id:packet.request_id,action:"input",reason:"reply_unavailable",automatic_retry:false}); return true;
  });
  f.terminal.enqueue("possibly_executed()\r");
  await until(() => f.terminal.state.uncertain && !f.terminal.busy);
  const count = f.sent.filter(p => p.action === "input").length;
  await delay(200);
  assert.equal(f.sent.filter(p => p.action === "input").length,count); assert.equal(f.terminal.state.canInput,false);
  assert.equal(f.terminal.input.length,0); assert.match(f.terminal.state.message,/not retype blindly/); checks+=4;
  f.hook(null); assert.equal(await f.terminal.connect(),true);
  assert.equal(f.sent.filter(p => p.action === "input").length,count); checks+=2;
  assert.equal(f.terminal.state.canInput,false); assert.equal(f.terminal.state.reviewRequired,true);
  assert.equal(f.terminal.enqueue("do_not_replay\r"),false);
  assert.equal(f.terminal.resumeInput(),true);assert.equal(f.terminal.state.canInput,true);checks+=5;
  f.client.state.stale = true; f.notify();
  assert.equal(f.terminal.state.canInput,false); assert.equal(f.terminal.socket,null); checks+=2;
  f.client.state.stale = false; f.client.state.assignments[0].id = randomUUID(); f.notify();
  assert.equal(f.terminal.session,null); assert.equal(f.terminal.writer,null); checks+=2;
  f.client.state.control.profiles[0].isolation = "native"; f.notify();
  assert.equal(f.terminal.state.canConnect,false); checks++;
} finally {f.close();}
assert.equal(f.listeners.size,0); checks++;

const detaching=fixture();let releaseDetach;
try {
  await detaching.terminal.connect();
  detaching.hook((packet,socket)=>{
    if(packet.action!=="disconnect")return false;
    releaseDetach=()=>socket.emit({kind:"report",request_id:packet.request_id,accepted:true,reason:"accepted",
      session_id:detaching.terminal.session,phase:"ready",writer_connected:false,input_sequence:0,
      cursor:0,output_sequence:0,gap:false,bytes:[],failure:null,cleanup_pending:false});
    return true;
  });
  detaching.terminal.disconnect();await until(()=>Boolean(releaseDetach));
  await delay(250);
  assert.equal(detaching.terminal.state.canConnect,false);
  assert.equal(detaching.terminal.state.canInput,false);
  assert.equal(detaching.terminal.state.phase,"disconnecting");checks++;
  assert.equal(await detaching.terminal.connect(),false);
  assert.equal(detaching.sent.filter(p=>p.action==="disconnect").length,1);checks+=4;
  releaseDetach();await until(()=>detaching.terminal.state.canConnect);
  assert.equal(detaching.terminal.phase,"disconnected");checks++;
}finally{detaching.close();}

let release;
const blocked = new Promise(resolve => {release = resolve;});
const slow = fixture({output:() => blocked});
try {
  await slow.terminal.connect();
  await until(() => slow.sent.some(p => p.action === "read"));
  const count = slow.sent.length;
  await delay(250); assert.equal(slow.sent.length,count,"render callback must bound further output reads"); checks++;
  slow.terminal.enqueue("queued\r"); slow.terminal.disconnect(); release();
  await until(() => !slow.terminal.busy);
  assert.equal(slow.terminal.input.length,0); assert.equal(slow.sent.some(p => p.action === "input"),false); checks+=2;
} finally {release();slow.close();}

for (const malformed of [{kind:"unknown"}, {kind:"attached",connection_id:randomUUID()}, "x".repeat(65537)]) {
  const bad = fixture();
  try {
    await bad.terminal.connect();
    if (typeof malformed === "string") bad.sockets[0].dispatchEvent(new MessageEvent("message",{data:malformed}));
    else bad.sockets[0].emit(malformed);
    assert.equal(bad.terminal.socket,null); checks++;
  } finally {bad.close();}
}
// Static security bounds supplement the behavioral fixture; no input persistence.
const source = await readFile(new URL("../../assets/runtime-terminal-client.js",import.meta.url),"utf8");
assert.doesNotMatch(source,/localStorage|sessionStorage|console\.|innerHTML|\.notify\(/); checks++;
console.log(`Runtime terminal transport: ${checks} assertions passed (explicit socket fixture).`);
