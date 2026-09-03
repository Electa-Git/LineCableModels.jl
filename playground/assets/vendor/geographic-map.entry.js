import Feature from "ol/Feature.js";
import KML from "ol/format/KML.js";
import Point from "ol/geom/Point.js";
import TileLayer from "ol/layer/Tile.js";
import VectorLayer from "ol/layer/Vector.js";
import Map from "ol/Map.js";
import View from "ol/View.js";
import { toLonLat, fromLonLat } from "ol/proj.js";
import OSM from "ol/source/OSM.js";
import VectorSource from "ol/source/Vector.js";
import CircleStyle from "ol/style/Circle.js";
import Fill from "ol/style/Fill.js";
import Stroke from "ol/style/Stroke.js";
import Style from "ol/style/Style.js";
import Text from "ol/style/Text.js";
import { unzipSync } from "fflate";

const KML_MIME = "application/vnd.google-earth.kml+xml";
const MAP_PROJECTION = "EPSG:3857";
const DATA_PROJECTION = "EPSG:4326";

function randomId() {
  if (globalThis.crypto && typeof globalThis.crypto.randomUUID === "function") {
    return globalThis.crypto.randomUUID();
  }
  return `reference-${Date.now()}-${Math.random().toString(16).slice(2)}`;
}

function archiveKml(bytes) {
  const archive = unzipSync(new Uint8Array(bytes));
  const names = Object.keys(archive).filter((name) => /\.kml$/i.test(name));
  const root = names.find((name) => name.toLowerCase() === "doc.kml");
  const selected = root || (names.length === 1 ? names[0] : null);
  if (!selected) {
    throw new Error(
      names.length === 0
        ? "KMZ archive does not contain a KML document"
        : "KMZ archive has multiple KML documents and no root doc.kml",
    );
  }
  return new TextDecoder("utf-8", { fatal: true }).decode(archive[selected]);
}

function sourceText(bytes, filename) {
  const data = new Uint8Array(bytes);
  const zipped = /\.kmz$/i.test(filename) ||
    (data.length >= 4 && data[0] === 0x50 && data[1] === 0x4b);
  return zipped
    ? archiveKml(bytes)
    : new TextDecoder("utf-8", { fatal: true }).decode(data);
}

function parseKml(bytes, filename) {
  const text = sourceText(bytes, filename);
  const document = new DOMParser().parseFromString(text, "application/xml");
  if (document.querySelector("parsererror")) {
    throw new Error("The KML document is not well-formed XML");
  }
  if (document.doctype) {
    throw new Error("KML documents containing a DTD are not accepted");
  }

  const format = new KML({ extractStyles: false, showPointNames: false });
  const parsed = format.readFeatures(text, {
    dataProjection: DATA_PROJECTION,
    featureProjection: MAP_PROJECTION,
  });
  const supported = [];
  let lineCount = 0;
  let pointCount = 0;
  let ignoredCount = 0;

  for (const feature of parsed) {
    const geometry = feature.getGeometry();
    const type = geometry && geometry.getType();
    if (type !== "LineString" && type !== "Point") {
      ignoredCount += 1;
      continue;
    }
    if (type === "LineString") lineCount += 1;
    if (type === "Point") pointCount += 1;
    const prefix = type === "LineString" ? "line" : "point";
    const ordinal = type === "LineString" ? lineCount : pointCount;
    const sourceId = String(feature.getId() ?? `${prefix}-${ordinal}`);
    feature.setId(sourceId);
    feature.set("lcmSourceId", sourceId, true);
    feature.set("lcmSourceType", type, true);
    supported.push(feature);
  }

  if (supported.length === 0) {
    throw new Error("KML contains no supported Point or LineString placemarks");
  }
  return { features: supported, lineCount, pointCount, ignoredCount };
}

const importedLineStyle = new Style({
  stroke: new Stroke({ color: "#159fc3", width: 4 }),
});

function importedPointStyle(feature) {
  const label = String(feature.get("name") || "").trim();
  return new Style({
    image: new CircleStyle({
      radius: 5,
      fill: new Fill({ color: "#dbdee0" }),
      stroke: new Stroke({ color: "#243031", width: 1.5 }),
    }),
    text: label
      ? new Text({
          text: label,
          offsetY: -14,
          font: "600 12px Lato, sans-serif",
          fill: new Fill({ color: "#f2f2f2" }),
          stroke: new Stroke({ color: "rgba(20, 26, 26, 0.9)", width: 3 }),
        })
      : undefined,
  });
}

function importedStyle(feature) {
  return feature.getGeometry()?.getType() === "LineString"
    ? importedLineStyle
    : importedPointStyle(feature);
}

function referenceStyle(feature) {
  const identifier = String(feature.get("identifier") || "").trim();
  return new Style({
    image: new CircleStyle({
      radius: feature.get("draft") ? 7 : 6,
      fill: new Fill({ color: feature.get("draft") ? "#f3bf61" : "#1abc9c" }),
      stroke: new Stroke({ color: "#102523", width: 2 }),
    }),
    text: identifier
      ? new Text({
          text: identifier,
          offsetY: -17,
          font: "700 12px JuliaMono, monospace",
          fill: new Fill({ color: "#f2f2f2" }),
          backgroundFill: new Fill({ color: "rgba(23, 27, 27, 0.88)" }),
          padding: [3, 5, 3, 5],
        })
      : undefined,
  });
}

function referenceRecord(feature) {
  const coordinate = toLonLat(feature.getGeometry().getCoordinates());
  return {
    uid: String(feature.get("uid")),
    source_feature: String(feature.get("sourceFeature")),
    longitude: Number(coordinate[0]),
    latitude: Number(coordinate[1]),
    identifier: String(feature.get("identifier")),
    tag: String(feature.get("tag")),
  };
}

function mountGeographicMap(root, configuration, referencePayload) {
  const mapElement = root.querySelector('[data-map-role="canvas"]');
  const addButton = root.querySelector('[data-map-role="add-reference"]');
  const status = root.querySelector('[data-map-role="status"]');
  const emptyState = root.querySelector('[data-map-role="empty"]');
  const editor = root.querySelector('[data-map-role="editor"]');
  const editorTitle = root.querySelector('[data-map-role="editor-title"]');
  const identifierInput = root.querySelector('[data-map-role="identifier"]');
  const tagInput = root.querySelector('[data-map-role="tag"]');
  const editorError = root.querySelector('[data-map-role="editor-error"]');
  const saveButton = root.querySelector('[data-map-role="save"]');
  const cancelButton = root.querySelector('[data-map-role="cancel"]');
  const deleteButton = root.querySelector('[data-map-role="delete"]');
  const referenceList = root.querySelector('[data-map-role="reference-list"]');
  const referenceCount = root.querySelector('[data-map-role="reference-count"]');
  const uploadElement = root.querySelector(".lc-upload-field");
  const defaultTag = configuration.tags[0]?.value || "reference";

  const importedSource = new VectorSource();
  const referenceSource = new VectorSource();
  const importedLayer = new VectorLayer({ source: importedSource, style: importedStyle });
  const referenceLayer = new VectorLayer({ source: referenceSource, style: referenceStyle });
  const layers = [];
  if (configuration.basemap === "osm") {
    layers.push(new TileLayer({ source: new OSM({ crossOrigin: "anonymous" }) }));
  }
  layers.push(importedLayer, referenceLayer);

  const map = new Map({
    target: mapElement,
    layers,
    view: new View({ center: fromLonLat([4.35, 50.85]), zoom: 7 }),
  });

  let placing = false;
  let activeReference = null;
  let destroyed = false;
  let loadGeneration = 0;
  let lastPublishedPayload = "";

  const setStatus = (message, kind = "ready") => {
    status.textContent = message;
    status.dataset.kind = kind;
  };

  const committedReferences = () =>
    referenceSource.getFeatures().filter((feature) => !feature.get("draft"));

  const publishReferences = () => {
    const payload = JSON.stringify(committedReferences().map(referenceRecord));
    lastPublishedPayload = payload;
    referencePayload.notify(payload);
  };

  const closeEditor = ({ discardDraft = true } = {}) => {
    if (discardDraft && activeReference?.get("draft")) {
      referenceSource.removeFeature(activeReference);
    }
    activeReference = null;
    editor.hidden = true;
    editorError.textContent = "";
    identifierInput.value = "";
    renderReferenceList();
  };

  const openEditor = (feature) => {
    if (activeReference && activeReference !== feature) closeEditor();
    activeReference = feature;
    const draft = Boolean(feature.get("draft"));
    editor.hidden = false;
    editorTitle.textContent = draft ? "New reference" : "Edit reference";
    identifierInput.value = String(feature.get("identifier") || "");
    tagInput.value = String(feature.get("tag") || defaultTag);
    deleteButton.hidden = draft;
    editorError.textContent = "";
    requestAnimationFrame(() => identifierInput.focus());
  };

  const editReference = (uid) => {
    const feature = referenceSource
      .getFeatures()
      .find((candidate) => candidate.get("uid") === uid);
    if (feature) openEditor(feature);
  };

  const renderReferenceList = () => {
    const references = committedReferences();
    referenceCount.textContent = String(references.length);
    referenceList.replaceChildren();
    if (references.length === 0) {
      const empty = document.createElement("p");
      empty.className = "lc-map-reference-empty";
      empty.textContent = "No references placed.";
      referenceList.append(empty);
      return;
    }
    for (const feature of references) {
      const button = document.createElement("button");
      button.type = "button";
      button.className = "lc-map-reference-row";
      const identity = document.createElement("strong");
      identity.textContent = feature.get("identifier");
      const tag = document.createElement("span");
      const option = configuration.tags.find(
        (candidate) => candidate.value === feature.get("tag"),
      );
      tag.textContent = option?.label || feature.get("tag");
      button.append(identity, tag);
      button.addEventListener("click", () => editReference(feature.get("uid")));
      referenceList.append(button);
    }
  };

  const replaceReferences = (payload) => {
    if (!payload || payload === lastPublishedPayload) return;
    let records;
    try {
      records = JSON.parse(payload);
      if (!Array.isArray(records)) return;
    } catch (_) {
      return;
    }
    closeEditor();
    referenceSource.clear();
    for (const record of records) {
      if (!Number.isFinite(record.longitude) || !Number.isFinite(record.latitude)) continue;
      const feature = new Feature({
        geometry: new Point(fromLonLat([record.longitude, record.latitude])),
      });
      feature.setProperties({
        uid: String(record.uid),
        sourceFeature: String(record.source_feature),
        identifier: String(record.identifier),
        tag: String(record.tag),
        draft: false,
      });
      referenceSource.addFeature(feature);
    }
    renderReferenceList();
  };

  const lineAtPixel = (pixel) =>
    map.forEachFeatureAtPixel(
      pixel,
      (feature, layer) => {
        if (layer !== importedLayer) return undefined;
        return feature.getGeometry()?.getType() === "LineString" ? feature : undefined;
      },
      { layerFilter: (layer) => layer === importedLayer, hitTolerance: 8 },
    );

  const referenceAtPixel = (pixel) =>
    map.forEachFeatureAtPixel(
      pixel,
      (feature, layer) => (layer === referenceLayer ? feature : undefined),
      { layerFilter: (layer) => layer === referenceLayer, hitTolerance: 7 },
    );

  const setPlacing = (enabled) => {
    placing = enabled;
    addButton.setAttribute("aria-pressed", String(enabled));
    addButton.classList.toggle("is-active", enabled);
    mapElement.classList.toggle("is-placing", enabled);
    setStatus(
      enabled ? "Click a route line to place a reference" : "Browse the imported geometry",
    );
  };

  const clearReferences = (notify = true) => {
    closeEditor();
    referenceSource.clear();
    renderReferenceList();
    if (notify) publishReferences();
  };

  const loadSource = async (url, filename) => {
    const generation = ++loadGeneration;
    setStatus("Reading geographic source…", "loading");
    try {
      const response = await fetch(url, { cache: "no-store", credentials: "same-origin" });
      if (!response.ok) throw new Error(`Source request failed (${response.status})`);
      const parsed = parseKml(await response.arrayBuffer(), filename);
      if (destroyed || generation !== loadGeneration) return;

      importedSource.clear();
      importedSource.addFeatures(parsed.features);
      clearReferences();
      emptyState.hidden = true;
      const extent = importedSource.getExtent();
      if (extent.every(Number.isFinite)) {
        map.getView().fit(extent, {
          padding: [36, 36, 36, 36],
          duration: 180,
          maxZoom: 18,
        });
      }
      const ignored = parsed.ignoredCount ? ` · ${parsed.ignoredCount} unsupported ignored` : "";
      setStatus(
        `${parsed.lineCount} line${parsed.lineCount === 1 ? "" : "s"} · ` +
          `${parsed.pointCount} point${parsed.pointCount === 1 ? "" : "s"}${ignored}`,
      );
    } catch (error) {
      if (destroyed || generation !== loadGeneration) return;
      setStatus(error instanceof Error ? error.message : "Unable to load geographic source", "error");
    }
  };

  addButton.addEventListener("click", () => setPlacing(!placing));
  cancelButton.addEventListener("click", () => closeEditor());
  deleteButton.addEventListener("click", () => {
    if (!activeReference) return;
    referenceSource.removeFeature(activeReference);
    activeReference = null;
    editor.hidden = true;
    publishReferences();
    renderReferenceList();
  });
  saveButton.addEventListener("click", () => {
    if (!activeReference) return;
    const identifier = identifierInput.value.trim();
    if (!identifier) {
      editorError.textContent = "ID is required.";
      identifierInput.focus();
      return;
    }
    const duplicate = committedReferences().some(
      (feature) =>
        feature !== activeReference &&
        String(feature.get("identifier")).toLocaleLowerCase() === identifier.toLocaleLowerCase(),
    );
    if (duplicate) {
      editorError.textContent = "ID must be unique.";
      identifierInput.focus();
      return;
    }
    activeReference.set("identifier", identifier);
    activeReference.set("tag", tagInput.value);
    activeReference.set("draft", false);
    activeReference.changed();
    editor.hidden = true;
    activeReference = null;
    editorError.textContent = "";
    publishReferences();
    renderReferenceList();
  });
  identifierInput.addEventListener("keydown", (event) => {
    if (event.key === "Enter") {
      event.preventDefault();
      saveButton.click();
    } else if (event.key === "Escape") {
      event.preventDefault();
      cancelButton.click();
    }
  });

  map.on("singleclick", (event) => {
    const reference = referenceAtPixel(event.pixel);
    if (reference) {
      openEditor(reference);
      return;
    }
    if (!placing) {
      closeEditor();
      return;
    }
    const line = lineAtPixel(event.pixel);
    if (!line) {
      setStatus("No route line within the placement tolerance", "warning");
      return;
    }
    closeEditor();
    const coordinate = line.getGeometry().getClosestPoint(event.coordinate);
    const feature = new Feature({ geometry: new Point(coordinate) });
    feature.setProperties({
      uid: randomId(),
      sourceFeature: String(line.get("lcmSourceId")),
      identifier: "",
      tag: defaultTag,
      draft: true,
    });
    referenceSource.addFeature(feature);
    openEditor(feature);
  });

  map.on("pointermove", (event) => {
    if (event.dragging) return;
    const overReference = Boolean(referenceAtPixel(event.pixel));
    const overLine = placing && Boolean(lineAtPixel(event.pixel));
    mapElement.classList.toggle("has-map-target", overReference || overLine);
  });

  root.addEventListener("lcm:upload-committed", (event) => {
    loadSource(event.detail.url, event.detail.file.original_name);
  });
  root.addEventListener("lcm:upload-removed", () => {
    loadGeneration += 1;
    importedSource.clear();
    clearReferences();
    emptyState.hidden = false;
    setStatus("Choose a KML or KMZ file", "empty");
  });

  referencePayload.on(replaceReferences);
  replaceReferences(referencePayload.value);
  renderReferenceList();

  const resize = new ResizeObserver(() => {
    requestAnimationFrame(() => requestAnimationFrame(() => map.updateSize()));
  });
  resize.observe(mapElement);

  if (uploadElement?.dataset.uploadState === "ready" && uploadElement.dataset.uploadName) {
    loadSource(uploadElement.dataset.uploadUrl, uploadElement.dataset.uploadName);
  }

  return {
    map,
    destroy() {
      destroyed = true;
      loadGeneration += 1;
      resize.disconnect();
      map.setTarget(undefined);
    },
  };
}

globalThis.LCMGeographicMap = Object.freeze({ mountGeographicMap, KML_MIME });

