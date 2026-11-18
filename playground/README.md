# LineCableModels playground

The playground is a Quarto website of presentations. Each presentation is a
folder `presentations/<name>/` with an `index.qmd`, and the landing page,
`index.qmd`, lists them. Reveal supplies the slide controller. The `lcm-deck`
extension in `_extensions/lcm-deck/` sets the layouts and the presentation
controls, and `assets/` holds the palette.

## Render and serve

From the folder that holds `playground/`:

```sh
julia --project=playground -m Playground
```

The command renders the project with the Quarto of `quarto_jll` and serves
`_site`. The site is then at `http://127.0.0.1:8080/`. The options `--host` and
`--port` set the address, and `--no-render` serves the site of the last render.
Every render passes `-M engine:markdown`, which tells Quarto to leave every code
block as text.

## Start a presentation

Copy `presentations/starter/` to a new folder of `presentations/`. Its front
matter selects the format:

```yaml
format: lcm-deck-revealjs
```

The landing page shows the `title`, `subtitle`, `description` and `image` of
each presentation's front matter. Quarto takes the description and the image
of a presentation from its front matter only, and each presentation states both.

Start each slide with a level-two heading, and nest no slide inside another.
Below each heading, put one named layout with the exact number of direct
`.lcm-slot` children that the table gives.

| Layout class | Slots | Geometry |
| --- | ---: | --- |
| `.lcm-layout-full-canvas` | 1 | one complete working region |
| `.lcm-layout-balanced` | 2 | equal columns |
| `.lcm-layout-feature-sidebar` | 2 | dominant region plus sidebar |
| `.lcm-layout-top-split` | 3 | top band over two columns |
| `.lcm-layout-dashboard-grid` | 4 | two-by-two grid |
| `.lcm-layout-media-story` | 2 | media plus narrative |

Add `.lcm-slot-panel` to give a slot the standard panel surface. Write
ordinary Markdown inside slots. Write a compact title inside a panel as
`[Title]{.lcm-panel-title}`. Avoid Markdown headings inside a slot, because
Pandoc can turn them into nested HTML sections.

`presentations/specimen/` shows every layout.

## Presenter controls

- The home icon at the bottom right returns to the landing page in the same tab.
- The arrow keys and the space bar move between slides while the focus is outside a form control.
- `S` opens the speaker notes of Reveal.
- `O` turns the overview on and off.
- `F` asks Reveal for full screen.
- `L` turns the red pointer on and off. It is the same red in both themes.
- `E` opens the PDF preview. A click on **Print / Save PDF** prints the slides, and a click on **Return to slides** goes back.

The status bar at the bottom reports whether the presentation is ready and
whether a resize is settling. The pointer leaves the slide where it is. The
theme selector chooses between dark, light and the system's choice. The
landing page and every presentation share the saved choice.

A focused input, select, button or slider keeps its arrow keys. A resize, a
change of orientation, the overview and full screen each start a short
settling state. In the overview, a click on a thumbnail selects that slide, and
`O` or Escape returns to it. The arrow keys move along the overview of Reveal.
The PDF preview is a separate static page. Its return button reloads the
presentation at the same slide.

The PDF menu, `E`, `?view=print` and `?print-pdf` all open the same preview.
It scrolls vertically through fitted 16:9 pages, and printing produces one
16:9 page per slide.

## Incremental lists

Put Quarto's `.incremental` fenced Div around an ordinary bullet or numbered
list. The items reveal in the order of the source, nested items included, and
a second list follows the first.

The slide **One clear argument** of the starter has a nested bullet list. The
slide **Top band plus two supporting regions** of the specimen has two lists,
nested bullets on the left and numbered steps on the right. Each step forward
reveals one item, and each step back hides it again. A hidden item keeps its
space, and print shows every item.

```markdown
::: {.incremental}
- State the question.
- Show the evidence.
  - Reveal supporting detail next.
- Explain the conclusion.
:::
```

## Equation explanations

Mark an exact term with MathJax's `\cssId`, then give a Markdown block the same
ID. Both belong to the same slide, inside a layout slot:

```markdown
$$
L = \frac{\cssId{impedance-imag}{\Im(Z)}}{\omega}
$$

::: {.lcm-math-note target="impedance-imag"}
**Imaginary part of impedance**

Dividing by angular frequency $\omega$ gives the inductance.
:::
```

An ID is unique across the presentation. It starts with a letter, followed by
letters, digits, hyphens or the character `_`. Write each `\cssId{...}{...}`
literally, never through another TeX macro. A note holds paragraphs, emphasis,
links, lists, code and inline math. `deck.lua` refuses raw HTML, images, nested
notes and Markdown headings in a note. Bold text gives a note its title.

The outlined term responds once MathJax has typeset the equation. The footer
reports whether the equation notes are loading, ready or unavailable. A click
on a term opens its explanation, and a second click closes it. A click on
another term switches to its explanation. A click outside, Escape, the close
button, a change of slide and the overview each close the explanation.

Tab moves the focus between terms. Enter or Space opens the explanation of the
focused term, and Escape returns the focus to the term. A text selection
leaves the term as it is, and so does the pointer.

An explanation takes up none of the slide's space and follows the palette of
the light and the dark theme. Popper, which Quarto bundles, places it within
the slide and moves it on resize. The deck adds the term's outline, the line
to its explanation and the response to clicks and keys. The code of the notes
is `_extensions/lcm-deck/math-notes.js` and `math-notes.css`, and `deck.lua`
checks the source of each note. The notes read the typesetting queue of
MathJax 2, the version that Quarto's Reveal format loads.

Quarto adds Popper while `citations-hover` or `footnotes-hover` is on, which is
their default. A presentation that turns both off shows its notes as plain
text.

The overview, the speaker view, the PDF preview and printing show the equation
without the outlines and the explanations. Print shows every incremental item.
A PDF leaves out the text of the explanations. Put on the slide itself
what the audience needs. When MathJax or Popper is missing, the footer reports
the notes unavailable and the slide shows the text of each explanation.

The slide **Balanced evidence** of the specimen has two notes.

## Code cells

A cell is a fenced code block of class `lcm-cell` inside a layout slot:

````markdown
::: {.lcm-slot}
```{.julia .lcm-cell autorun="true" echo="false"}
Z = 0.05 + 0.31im
```
:::
````

The braces start with a dot, so Quarto selects no engine and leaves the block
to `deck.lua`. A cell takes two attributes, each `true` or `false`: `autorun`,
`false` when absent, and `echo`, `true` when absent.

`deck.lua` turns each cell into a live frame that loads `live/cells/<n>`
relative to the page, where `n` counts the cells of the presentation from 1.
The cell's highlighted source takes the frame's place in the print view, the
speaker view and the overview. It also shows when the live view is unavailable.
Frames in one panel share the height that the panel's other content leaves.

Beside the rendered page, `deck.lua` writes `cells.json`, one entry per cell in
the order of the source: `index`, `slide`, `source`, `autorun` and `echo`.
`slide` counts the level-two headings from 1. A presentation without cells
gets an empty list.

`deck.lua` refuses the following, with a message that gives the slide:

- a block written `{julia}`, which belongs to the Julia engine of Quarto
- an attribute other than `autorun` and `echo`
- a value other than `true` and `false`
- a cell outside a `.lcm-slot`

The slide **Balanced evidence** of the specimen has two cells. The first
defines `Z`, and the second computes `R` and `L` from it.
