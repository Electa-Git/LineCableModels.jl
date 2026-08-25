# Plotting smoke checks

The automated Cairo suite covers deterministic rendering, responsive docks,
callbacks, and SVG reconstruction. Changes to interactive window creation,
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

Resize the compact cable-preview window from its initial size to a tall window
and back to a short window. Confirm that:

- axes and canvas remain inside the window;
- the legend stays in the 220-pixel side dock;
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
