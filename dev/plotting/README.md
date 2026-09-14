# Plotting smoke checks

The automated Cairo suite covers deterministic rendering, responsive docks,
callbacks, and live-state SVG export. Changes to interactive window creation,
display, layout, resizing, or callbacks also require a real-display GL check.

Instantiate the development environment once:

```sh
julia --project=dev/plotting -e 'using Pkg; Pkg.instantiate()'
```

Run the automated GL and SVG gate under a real display or Xvfb:

```sh
julia --project=dev/plotting dev/plotting/manual_gl.jl
```

Run the complete interactive gallery under a real display:

```sh
julia --project=dev/plotting dev/plotting/manual_gl_gallery.jl
```

Run the focused cable-collection preview gallery:

```sh
julia --project=dev/plotting dev/plotting/manual_gl_cable_collection.jl
```

It opens an automatically arranged 2×3 canvas and an explicit 1×4 canvas.
Confirm that each cable id appears as its subplot title, neither window has a
legend, and each window has one shared set of three material colorbars.

Run the focused Monte Carlo distribution gallery:

```sh
julia --project=dev/plotting dev/plotting/manual_gl_monte_carlo.jl
```

It opens the native Makie histogram, density, empirical-CDF, model-CDF, and Q-Q
verbs, a composed samples/model page, and one indexed line-parameter marginal.
Read the comments beside each call to follow observable publication, direct
Makie primitives, and composition in the standard shell.

Resize the compact cable-preview window from its initial size to a tall window
and back to a short window. Confirm that:

- axes and canvas remain inside the window;
- the legend stays in its selected native layout dock;
- hidden legend entries return when vertical space increases;
- `(...)` appears when entries overflow and disappears when they fit;
- all three material colorbars stay inside the side dock;
- the toolbar and status row retain their positions;
- reset, log toggles, visibility toggles, and SVG export remain functional.

For the comparison grid, run:

```sh
julia --project=dev/plotting dev/plotting/manual_gl_comparison.jl
```

Resize each window and confirm that the matrix grid and responsive legend fit
without clipping. Toggle one legend entry and confirm that the same scenario
is hidden in every panel.
