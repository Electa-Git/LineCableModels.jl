# PlotBuilder architecture decision

PlotBuilder follows five non-negotiable rules:

1. `make_render` is one generic pipeline. Domain recipes specialize accessors;
   they do not replace the pipeline.
2. Recipe variation uses Julia dispatch, including `Val` dispatch for modes,
   grouping, and primitive rendering.
3. `PlotRecipe`, `RenderSpec`, layouts, axes, series, views, and page components
   are backend-neutral. Makie objects exist only in the Makie extension.
4. Scientific selection, units, labels, grouping, and page identity remain
   declarative. Interactive state is reconstructed into the same typed model
   for export.
5. Layouts are named grid trees. Recipes provide defaults and callers may
   supply a preset or `LayoutSpec` without changing a recipe implementation.

These rules apply to every maintained plot family. The user-facing `plot`,
`preview`, `show_material_scale`, and `export_svg` methods remain narrow
adapters over the same recipe and renderer path.

The developer API may evolve before LineCableModels 1.0, but changes must keep
the separation between domain recipes, backend-neutral specifications, and
backend rendering. Architectural tests enforce the single pipeline and the
absence of Makie from core specification construction.
