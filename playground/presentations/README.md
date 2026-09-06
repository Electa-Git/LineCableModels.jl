# LCM presentation authoring

The presentation layer is a Quarto format. It does not import Julia packages,
open broker connections, or execute document code. Reveal supplies controller
features; the `lcm-deck` extension owns layout and live-view lifecycle.

## Start a deck

Use `starter.qmd` as the smallest complete source. Its required format is:

```yaml
format: lcm-deck-revealjs
```

The gallery separates **complete decks** from **reusable slide layouts**.
The six layout cards are building blocks for individual slides, not separate
applications. Each card links to its rendered specimen slide and to a copyable
Markdown block in `layouts.qmd` (`/presentations/layouts.html` in the site).
The starter uses full canvas and feature + sidebar; it can use any of the six.

Create slides with level-two headings. Keep the sequence flat. Put one named
layout below each heading and use the exact number of direct `.lcm-slot`
children declared below.

| Layout class | Slots | Geometry |
|---|---:|---|
| `.lcm-layout-full-canvas` | 1 | one complete working region |
| `.lcm-layout-balanced` | 2 | equal columns |
| `.lcm-layout-feature-sidebar` | 2 | dominant region plus sidebar |
| `.lcm-layout-top-split` | 3 | top band over two columns |
| `.lcm-layout-dashboard-grid` | 4 | two-by-two grid |
| `.lcm-layout-media-story` | 2 | media plus narrative |

Add `.lcm-slot-panel` when a slot needs the standard panel surface. Use
ordinary Markdown inside slots. A compact in-panel title is written as
`[Title]{.lcm-panel-title}`; nested Markdown headings are deliberately avoided
because Pandoc may materialize them as nested HTML sections.

## Embed a live view

Only the validated same-origin shortcode may create a live boundary:

```markdown
{{< bonito route="/widgets/control-panel"
    title="Cable controls"
    height="100%"
    public-url="https://example.org/playground/widgets/control-panel" >}}
```

`route` is the same-origin runtime route. `public-url` is the durable link used
by speaker previews and PDF output. The audience iframe mounts on first slide
entry and then remains mounted. Static modes never activate it.

## Presenter controls

- Home icon in the control bar, or **Menu → Tools → Playground home**: return to the playground landing page in the same tab.
- Arrow keys or space: navigate when focus is outside a form control.
- `S`: open Reveal speaker notes.
- `O`: toggle overview.
- `F`: request fullscreen through Reveal.
- `L`: toggle the red LCM pointer (the same red in both themes).
- `E`: open PDF preview; use **Print / Save PDF** there, or **Return to slides**.

The compact status bar reports presentation readiness, live-view loading, and
resize settling separately from the laser's on/off state. The laser works
immediately, including over same-origin live views, and does not advance slides.
The theme selector uses the same saved preference and palette as the playground
and workbenches; loaded widgets follow changes without replacing their sessions.

Focused inputs, selects, buttons, and sliders keep their navigation keys.
Resize, orientation, overview, and fullscreen changes enter a short settling
state before live children receive final viewport dimensions.
Overview retains each live session but shows its inert linked placeholder in
the transformed thumbnail. Click a thumbnail to select it, or press O/Escape to
return. The current slide is centered; arrow keys move along Reveal's overview.
Browser printing temporarily substitutes placeholders and restores the live
frame on return. PDF preview is a separate static page: its return button reloads
the audience page at the same slide. It does not promise session retention.

The PDF menu, E, `?view=print`, and `?print-pdf` all enter the same LCM preview.
It uses normal vertical scrolling with fitted 16:9 pages; printing produces one
16:9 page per slide. Live views are links, not screenshots. Use absolute
`public-url` values for published PDFs so links remain valid away from the
local server. The CLI export uses a temporary server, so relative links are
only appropriate for browser previews, not distributed CLI-generated PDFs.

Restart the publisher after rebuilding: Quarto can generate new hashed CSS
filenames, and the publisher registers static routes when it starts.

## Commands

```sh
lcm presentation check presentations/starter.qmd
lcm presentation build presentations/starter.qmd
lcm presentation start presentations/starter.qmd
lcm presentation export presentations/starter.qmd --pdf
```

Run `test/integration/run-presentation.sh` before migrating a scientific deck.
The gate exercises five display sizes and all layouts, real overview/PDF menu
clicks, exact canvas clicks, focused input keys, delayed-child readiness, laser
tracking inside frames, cached/live theme changes, session persistence, static
speaker receivers, and actual Chrome PDF output.
