# Vendored browser assets

## Private Julia terminal

`runtime-terminal.bundle.js` and `runtime-terminal.bundle.css` are generated from
`runtime-terminal.entry.js` using xterm.js **6.0.0**, its fit addon **0.11.0** and
esbuild **0.25.9**. They contain the renderer and fitting utility only, not the
upstream demo/attach server or addon. The owned runtime component supplies the
authorized socket, bounded I/O, theme projection and cleanup. The shared browser
component is `JuliaTerminal`; these assets alone do not enable terminal execution.

The generated assets are served locally, without runtime npm or CDN access.
Upstream MIT notices are retained in `licenses/xterm-LICENSE.txt` and
`licenses/xterm-addon-fit-LICENSE.txt`. On dependency upgrades review the owned
clipboard, link, escape-sequence, keyboard and backpressure policy against the
[upstream security guidance](https://xtermjs.org/docs/guides/security/) and rerun
the complete terminal browser/relay gates.

To rebuild from `playground/` in a fresh temporary npm prefix:

```sh
LCM_TERMINAL_BUILD=$(mktemp -d /tmp/lcm-terminal-build.XXXXXXXX)
npm install --prefix "$LCM_TERMINAL_BUILD" --no-save --ignore-scripts --no-audit --no-fund \
  @xterm/xterm@6.0.0 @xterm/addon-fit@0.11.0 esbuild@0.25.9
"$LCM_TERMINAL_BUILD/node_modules/.bin/esbuild" assets/vendor/runtime-terminal.entry.js \
  --bundle --format=iife --target=es2022 --minify \
  --alias:@xterm/xterm="$LCM_TERMINAL_BUILD/node_modules/@xterm/xterm" \
  --alias:@xterm/addon-fit="$LCM_TERMINAL_BUILD/node_modules/@xterm/addon-fit" \
  --outfile=assets/vendor/runtime-terminal.bundle.js
```

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
