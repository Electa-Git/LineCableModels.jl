# Manual plotting checks

The automated Cairo suite covers deterministic rendering, responsive docks,
callbacks, and live-state SVG export. Changes to interactive window creation,
display, layout, resizing, or callbacks also require a real-display GL check.

These checks are deliberately absent from CI. Use the installed development
dependencies without replacing the active project:

```sh
export JULIA_LOAD_PATH="@:$PWD/test/visual:$PWD/test/manual/plotting:@v1.12:@stdlib"
```

Run the GL and SVG inspection under a real display:

```sh
julia --project=. test/manual/plotting/manual_gl.jl
julia --project=. test/manual/plotting/manual_gl.jl --backends-first
```

Both runs begin with CairoMakie installed but unloaded: GL controls work and
SVG controls are absent. They then explicitly import CairoMakie, export the old
window through the API, and recreate a GL window with a working SVG button.
Export preserves the live view and does not switch display backends.

Run the complete interactive gallery under a real display:

```sh
julia --project=. test/manual/plotting/manual_gl_gallery.jl
```

Run the focused cable-collection preview gallery:

```sh
julia --project=. test/manual/plotting/manual_gl_cable_collection.jl
```

It opens an automatically arranged 2×3 canvas and an explicit 1×4 canvas.
Confirm that each cable id appears as its subplot title, neither window has a
legend, and each window has one shared set of three material colorbars.

Run the focused Monte Carlo distribution gallery:

```sh
julia --project=. test/manual/plotting/manual_gl_monte_carlo.jl
```

It opens five native Makie views: histogram, density, empirical CDF, model CDF,
and Q–Q. Its synthetic completed fixture needs no new Monte Carlo calculation.

Resize the gallery windows from their initial size to a tall window and back.
Confirm that:

- axes and canvas remain inside the window;
- the legend stays in its selected native layout dock;
- automatic `:show_all` legends wrap without omitting entries;
- all three material colorbars stay inside the side dock;
- the toolbar and status row retain their positions;
- reset, log toggles, and visibility toggles remain functional;
- SVG export is available when CairoMakie was loaded before plot construction.

For the comparison grid, run:

```sh
julia --project=. test/manual/plotting/manual_gl_comparison.jl
```

Resize each window and confirm that the matrix grid and responsive legend fit
without clipping. Toggle one legend entry and confirm that the same scenario
is hidden in every panel.

Before the refactor, `blocks` paginated matrices and `layout` could pair
quantities. Now `layout` is the only panel capacity and each quantity has its
own figures. Result plots use `fig_size`; previews and `plotwindow` retain
`size`. The complete script and saved-run disposition is in [../README.md](../README.md).
