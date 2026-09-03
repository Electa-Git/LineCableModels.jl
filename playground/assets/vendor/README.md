# Vendored browser assets

## Geographic map

`geographic-map.bundle.js` is a browser-ready bundle generated from
`geographic-map.entry.js` with these pinned dependencies:

- OpenLayers 10.10.0
- fflate 0.8.2

The generated bundle is committed so running the Julia playground does not
require Node.js, a package registry, or a JavaScript CDN. `openlayers.css` is
copied from the same OpenLayers release. Upstream license texts live in
`licenses/`.

To rebuild in a temporary npm prefix:

```sh
npm install --prefix /tmp/lcm-geographic-assets --no-save --ignore-scripts \
  ol@10.10.0 fflate@0.8.2 esbuild@0.25.9

/tmp/lcm-geographic-assets/node_modules/.bin/esbuild \
  assets/vendor/geographic-map.entry.js \
  --bundle --format=iife --target=es2022 --minify \
  --alias:ol=/tmp/lcm-geographic-assets/node_modules/ol \
  --alias:fflate=/tmp/lcm-geographic-assets/node_modules/fflate \
  --outfile=assets/vendor/geographic-map.bundle.js
```

## Power-system canvas

`power-system-canvas.bundle.js` and `power-system-canvas.bundle.css` are built
from `power-system-canvas.entry.jsx` with these pinned dependencies:

- sldeditor 0.23.0
- React 19.1.1
- React DOM 19.1.1
- esbuild 0.25.9

The generated assets are committed so the publisher neither downloads code nor
requires Node.js at runtime. The adapter deliberately mounts one editor per
browser document and replaces sldeditor's persistent Zustand storage with a
volatile store; each Bonito widget route is already isolated in its own iframe.

To rebuild in a temporary npm prefix:

```sh
npm install --prefix /tmp/lcm-sld-assets --no-save --ignore-scripts \
  sldeditor@0.23.0 react@19.1.1 react-dom@19.1.1 esbuild@0.25.9

/tmp/lcm-sld-assets/node_modules/.bin/esbuild \
  assets/vendor/power-system-canvas.entry.jsx \
  --bundle --format=iife --target=es2022 --minify \
  --alias:sldeditor/style.css=/tmp/lcm-sld-assets/node_modules/sldeditor/dist/style.css \
  --alias:sldeditor=/tmp/lcm-sld-assets/node_modules/sldeditor \
  --alias:react=/tmp/lcm-sld-assets/node_modules/react \
  --alias:react-dom=/tmp/lcm-sld-assets/node_modules/react-dom \
  --outfile=assets/vendor/power-system-canvas.bundle.js
```
