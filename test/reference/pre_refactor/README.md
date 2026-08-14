# Pre-refactor plotting references

These images were rendered from commit
`abe355170add20cc0c841e218fed82fe8ffbab51`, immediately before the result and
PlotBuilder consolidation. They are retained as the human visual-parity
baseline; the PNG files in the parent directory are the executable tolerant
goldens for the consolidated renderer.

The preserved line-parameter renderer contained two unqualified calls to
`value` that prevented it from rendering in its recorded environment. For
baseline capture only, the temporary checkout changed those calls to
`Measurements.value`. No plot data, specification, theme, layout, or style was
changed.

- `line_rlcg.png`: series-resistance page.
- `line_zy_cartesian.png`: Cartesian series-impedance page.
- `line_zy_polar.png`: polar series-impedance-magnitude page.
- `cable_preview.png`: cable-design preview.
