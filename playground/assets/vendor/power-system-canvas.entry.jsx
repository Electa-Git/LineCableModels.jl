import React, { useEffect, useState } from "react";
import { createRoot } from "react-dom/client";
import {
  OneLineEditor,
  compile,
  fitToContentSoon,
  useEditorStore,
} from "sldeditor";
import "sldeditor/style.css";

const SNAPSHOT_SCHEMA = "lcm.power-system-topology/1";

// The Julia wrapper primes this key before this module loads, so the gallery
// opens on the specimen rather than sldeditor's first-run card. Remove the key
// immediately: the embedded editor does not own persistent browser state.
try {
  window.localStorage.removeItem("ole-onboarding-dismissed");
} catch (_) {
  // Storage may be disabled; the editor itself supports that environment.
}

function cloneDiagram(diagram) {
  return globalThis.structuredClone
    ? globalThis.structuredClone(diagram)
    : JSON.parse(JSON.stringify(diagram));
}

function resolvedTheme() {
  return document.documentElement.dataset.lcmResolvedTheme === "light"
    ? "light"
    : "dark";
}

function Editor({ diagram, readOnly }) {
  const [theme, setTheme] = useState(resolvedTheme);

  useEffect(() => {
    const observer = new MutationObserver(() => setTheme(resolvedTheme()));
    observer.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["data-lcm-resolved-theme"],
    });
    return () => observer.disconnect();
  }, []);

  return (
    <OneLineEditor
      className="h-full w-full"
      diagram={diagram}
      readOnly={readOnly}
      locale="en"
      theme={theme}
    />
  );
}

function topologySnapshot(diagram) {
  const model = compile(diagram);
  return {
    schema: SNAPSHOT_SCHEMA,
    diagram,
    topology: {
      nodes: Array.from(model.nodes.values(), (node) => ({
        id: node.id,
        terminals: Array.from(node.terminals),
      })),
      terminal_to_node: Array.from(model.terminalToNode.entries(), ([terminal, node]) => ({
        terminal,
        node,
      })),
      diagnostics: model.diagnostics.map((diagnostic) => ({ ...diagnostic })),
    },
  };
}

function summaryText(snapshot) {
  const diagram = snapshot.diagram;
  const diagnostics = snapshot.topology.diagnostics;
  const diagnosticText = diagnostics.length
    ? ` · ${diagnostics.length} diagnostic${diagnostics.length === 1 ? "" : "s"}`
    : "";
  return `${diagram.buses?.length || 0} buses · ${diagram.elements.length} devices · ` +
    `${snapshot.topology.nodes.length} electrical nodes${diagnosticText}`;
}

function mount(editorElement, configuration, topologyPayload, summaryElement) {
  const initialDiagram = cloneDiagram(configuration.diagram);
  let destroyed = false;
  let scheduledFrame = 0;
  let lastDiagram = null;

  // sldeditor deliberately uses one Zustand store per browser document. Bonito
  // gives every live widget its own document, and this adapter replaces the
  // package's localStorage persistence with a volatile store for that document.
  useEditorStore.persist.setOptions({
    storage: {
      getItem: () => null,
      setItem: () => undefined,
      removeItem: () => undefined,
    },
  });
  useEditorStore.getState().setDiagram(cloneDiagram(initialDiagram));

  const publish = (diagram) => {
    if (destroyed) return;
    const snapshot = topologySnapshot(diagram);
    summaryElement.dataset.kind = snapshot.topology.diagnostics.some(
      (diagnostic) => diagnostic.severity === "error",
    )
      ? "error"
      : "ready";
    summaryElement.textContent = summaryText(snapshot);
    topologyPayload.notify(JSON.stringify(snapshot));
  };

  const schedulePublish = (state) => {
    if (state.diagram === lastDiagram) return;
    lastDiagram = state.diagram;
    if (scheduledFrame) cancelAnimationFrame(scheduledFrame);
    scheduledFrame = requestAnimationFrame(() => {
      scheduledFrame = 0;
      publish(state.diagram);
    });
  };

  const unsubscribe = useEditorStore.subscribe(schedulePublish);
  const reactRoot = createRoot(editorElement);
  reactRoot.render(
    <Editor diagram={initialDiagram} readOnly={Boolean(configuration.readOnly)} />,
  );
  schedulePublish(useEditorStore.getState());
  requestAnimationFrame(() => fitToContentSoon());

  return {
    reset() {
      useEditorStore.getState().setDiagram(cloneDiagram(initialDiagram));
      requestAnimationFrame(() => fitToContentSoon());
    },
    snapshot() {
      return topologySnapshot(useEditorStore.getState().diagram);
    },
    destroy() {
      if (destroyed) return;
      destroyed = true;
      if (scheduledFrame) cancelAnimationFrame(scheduledFrame);
      unsubscribe();
      reactRoot.unmount();
    },
  };
}

globalThis.LCMPowerSystemCanvas = Object.freeze({ mount });
