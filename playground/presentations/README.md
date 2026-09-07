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

- Home icon in the footer control bar: return to the playground landing page in the same tab. Tools contains presentation actions only.
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

## Incremental lists

Use Quarto's native `.incremental` fenced Div around an ordinary bullet or
numbered list. No LCM ordering parameter is needed: items, including nested
items, reveal in document order; a second list follows the first one.

Try **One clear argument** in `starter.qmd` for a nested bullet list, or
**Top band plus two supporting regions** in `specimen.qmd` for two lists on
one slide (nested bullets on the left, numbered steps on the right). Use the
next/previous presentation controls to reveal/hide one item at a time. Hidden
items keep their allocated space; print output shows every item.

```markdown
::: {.incremental}
- State the question.
- Show the evidence.
  - Reveal supporting detail next.
- Explain the conclusion.
:::
```

## Equation explanations

Mark an exact term using MathJax's `\cssId`, then associate a Markdown block
with the same ID. Both must be on the same slide, inside a layout slot:

```markdown
$$
L = \frac{\cssId{impedance-imag}{\Im(Z)}}{\omega}
$$

::: {.lcm-math-note target="impedance-imag"}
**Imaginary part of impedance**

Dividing by angular frequency $\omega$ gives the inductance.
:::
```

IDs must be unique across the deck: start with a letter, then use letters,
digits, hyphens or underscores. Use literal `\cssId{...}{...}` declarations,
not IDs generated by another TeX macro. Notes accept ordinary paragraphs,
emphasis, links, lists, code and inline math; no raw HTML, images, nested
notes or Markdown headings. Use bold text for an in-note title.

The outlined term becomes interactive after MathJax finishes typesetting.
The footer reports loading, ready, or unavailable equation notes. Click a term
to toggle its explanation; another term switches explanations. Outside click,
Escape, the close button, slide changes and overview close the callout. Tab
focuses terms; Enter/Space opens a non-modal explanation, and Escape returns
focus. Selecting text does not activate the term. No hover trigger is needed
in this first version; moving the laser does not open notes.

Callouts reserve no slide space and consume shared light/dark palette tokens.
Quarto's bundled Popper positions them within the slide and updates on resize;
LCM adds only the term marker, bracket/leader and interaction lifecycle. The
component lives in `_extensions/lcm-deck/math-notes.{js,css}`, with authoring
validation in `deck.lua`. It does not import application state or inspect
live-widget documents. The current adapter targets MathJax 2's typesetting
queue, matching this format's existing renderer, without changing equations.
Do not remove Quarto's bundled Popper script when customizing the format.

Overview, speaker receivers, PDF preview and browser/CLI printing omit
transient explanations and term outlines. The equation remains unchanged;
all native incremental items are visible in print. Explanation prose is not
included in PDFs: put essential information on the slide or in a separate
handout. A missing renderer/dependency reports unavailability and exposes the
unbound explanation text on the audience slide instead of silently losing it.

See **Balanced evidence** in `specimen.qmd` for two working term annotations.

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
