# Presentation starter

Copy this entire directory to a URL-safe slug under `presentations/`.
Lowercase hyphenated names are conventional; uppercase is supported for
acronyms:

```sh
cp -a showcase/bonito/presentations/_template \
  showcase/bonito/presentations/my-tutorial
```

The `_template` directory is deliberately ignored by normal presentation
discovery. Once copied to `my-tutorial`, restart the server and the playground
will discover it automatically.

`content/010_starter.jl` is executable documentation. It demonstrates:

- Julia Markdown, code, tables, inline mathematics, and display mathematics;
- native sliders, dropdowns, buttons, mapped values, and action callbacks;
- an explicit process-wide preflight resource with `cold → preparing → hot`
  state;
- one, two, and three-column layouts;
- full-width top and bottom regions;
- control/plot, two-row, sidebar/main, and four-panel layouts;
- a persistent reactive WGLMakie figure; and
- presentation-local images and captions.

Rename the module, deck ID, group, title, and page IDs. Delete the examples you
do not need and replace the remaining page bodies with real content. Keep
resource construction in `setup(session)` when multiple pages share it.

`PREFLIGHT_RESOURCE` is a separate pattern for expensive server-process state
that must be ready before **Open** becomes available. Replace its activation
body with real preparation, or remove `resources = (PREFLIGHT_RESOURCE,)` when
the presentation has nothing to warm. Do not put session widgets or figures in
preflight.

Do not copy layout implementation or CSS into the presentation. The stable
layout functions live in `PageAuthoring`; `webgrid` remains available for a
composition that genuinely falls outside the supplied masters.
