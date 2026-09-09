```@meta
EditURL = "../literate/gauntlet.jl"
```

# Gauntlet

Gauntlet runs cases manually and retains numerical results. This page is its
only report: a compact comparison of the **applicable baseline formulations**
for each case, read from explicitly configured, persisted benchmarks.
Documentation generation never starts PSCAD, FEM, analytical calculations or
uncertainty propagation. Full-catalogue and UQ results remain stored for
analysis through the ordinary result, observation and plotting APIs.

## Recorded baseline comparisons

Both errors are calculated element-wise by
[`compare`](@ref LineCableModels.Engine.compare), with **A = the benchmark's
reference** and **B = its candidate**. Neither role is inferred from a backend.

```math
\mathrm{NRMSE}=100\sqrt{\frac{\sum_k|B_k-A_k|^2}{\sum_k|A_k|^2}},
\qquad
\mathrm{RMS}_{\mathrm{pointwise}}=100\sqrt{\frac1N\sum_k\left|\frac{B_k-A_k}{A_k}\right|^2}.
```

**One row is one benchmark.** Reference and candidate columns identify the
backend and requested formulation selection. Z and Y stay side by side;
each cell contains **maximum error in percent (response, excitation)**.
The maximum is across matrix entries, not a whole-matrix error. The two
normalizations can attain their maxima at different entries. Expand the
benchmark identities below each full-band table to see the terminal order.

**The full-band summary comes first.** Frequency slices follow in separate
sections, with the same column layout and their actual stored ranges.
No interpolation or additional simulations are performed.

46 explicitly configured benchmarks across 19 cases. Only completed, checksum-verified operands are included; no backend is selected as a reference by this page.

Bracketed numbers identify the requested formulation selections:

```@raw html
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">selection</th><th style = "text-align: left;">formulation</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: right;">1</td><td style = "text-align: left;">all slots :default</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: right;">2</td><td style = "text-align: left;">insulation_admittance=Ametani2004, semicon_admittance=Ametani2004; remaining slots :default</td></tr></tbody></table></div>
```

### Full-band Z/Y summary

Stored range: **0.1–1.0e6 Hz**, **101 samples**.

```@raw html
<h4>132 kV 630 mm² cables in flat horizontal formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">4.482 (4, 4)</td><td style = "text-align: left;">4.797 (4, 4)</td><td style = "text-align: left;">1107.0 (9, 8)</td><td style = "text-align: left;">1107.0 (8, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.501 (7, 1)</td><td style = "text-align: left;">0.3885 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">4.482 (4, 4)</td><td style = "text-align: left;">4.797 (4, 4)</td><td style = "text-align: left;">1213.0 (9, 8)</td><td style = "text-align: left;">3093.0 (5, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.501 (7, 1)</td><td style = "text-align: left;">0.3885 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_132kv_630mm2_flathor</code>.</p><ol>
<li><code>dielectric_baseline/cable_132kv_630mm2_flathor__default__fem_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_630mm2_flathor__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_630mm2_flathor__ametani2004__fem_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_630mm2_flathor__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2, 3, 4): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:jacket</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:jacket</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:jacket</code>.</p>
</details>

```

```@raw html
<h4>CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">6.757 (2, 1)</td><td style = "text-align: left;">3.753 (2, 1)</td><td style = "text-align: left;">100.0 (3, 6)</td><td style = "text-align: left;">100.0 (4, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7179 (3, 1)</td><td style = "text-align: left;">0.1831 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">12.9 (2, 1)</td><td style = "text-align: left;">4.827 (2, 1)</td><td style = "text-align: left;">100.0 (5, 3)</td><td style = "text-align: left;">100.0 (1, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7179 (3, 1)</td><td style = "text-align: left;">0.1831 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_132kv_cigre_tb880_case0_630cu_trefoil</code>.</p><ol>
<li><code>dielectric_baseline/cable_132kv_cigre_tb880_case0_630cu_trefoil__default__fem_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_cigre_tb880_case0_630cu_trefoil__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_cigre_tb880_case0_630cu_trefoil__ametani2004__fem_reference</code></li>
<li><code>dielectric_baseline/cable_132kv_cigre_tb880_case0_630cu_trefoil__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2, 3, 4): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>, <code>5=cable:3:core</code>, <code>6=cable:3:sheath</code>.</p>
</details>

```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">7.837 (1, 1)</td><td style = "text-align: left;">5.419 (1, 1)</td><td style = "text-align: left;">711.0 (6, 5)</td><td style = "text-align: left;">711.0 (5, 6)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.1789 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">7.838 (1, 1)</td><td style = "text-align: left;">5.419 (1, 1)</td><td style = "text-align: left;">797.2 (6, 5)</td><td style = "text-align: left;">1908.0 (5, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.1789 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_18kv_1000mm2_trefoil</code>.</p><ol>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil__default__fem_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil__ametani2004__fem_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2, 3, 4): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:jacket</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:jacket</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:jacket</code>.</p>
</details>

```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil — homogenized assembly</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">6.54 (2, 1)</td><td style = "text-align: left;">3.482 (3, 1)</td><td style = "text-align: left;">100.0 (1, 3)</td><td style = "text-align: left;">109.6 (9, 7)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.1789 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">6.54 (2, 1)</td><td style = "text-align: left;">3.482 (3, 1)</td><td style = "text-align: left;">100.0 (9, 7)</td><td style = "text-align: left;">435.3 (1, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.1789 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_18kv_1000mm2_trefoil_homogenized</code>.</p><ol>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil_homogenized__default__fem_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil_homogenized__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil_homogenized__ametani2004__fem_reference</code></li>
<li><code>dielectric_baseline/cable_18kv_1000mm2_trefoil_homogenized__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2, 3, 4): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:jacket</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:jacket</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:jacket</code>.</p>
</details>

```

```@raw html
<h4>220 kV Milliken 2500 mm² Al / 252 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.8387 (4, 1)</td><td style = "text-align: left;">0.2139 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.8387 (4, 1)</td><td style = "text-align: left;">0.2139 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_220kv_eaxecew_1x2500_252_trefoil</code>.</p><ol>
<li><code>dielectric_baseline/cable_220kv_eaxecew_1x2500_252_trefoil__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_220kv_eaxecew_1x2500_252_trefoil__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:foil</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:foil</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:foil</code>.</p>
</details>

```

```@raw html
<h4>18/30 kV NA2XS2Y 630 mm² Al cables in touching trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.6601 (3, 1)</td><td style = "text-align: left;">0.1684 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.6601 (3, 1)</td><td style = "text-align: left;">0.1684 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_30kv_na2xs2y_630mm2_trefoil</code>.</p><ol>
<li><code>dielectric_baseline/cable_30kv_na2xs2y_630mm2_trefoil__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_30kv_na2xs2y_630mm2_trefoil__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>, <code>5=cable:3:core</code>, <code>6=cable:3:sheath</code>.</p>
</details>

```

```@raw html
<h4>320 kV armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.022 (4, 1)</td><td style = "text-align: left;">2.065 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.022 (4, 1)</td><td style = "text-align: left;">2.065 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_320kv_armoured_dc_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_320kv_armoured_dc_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_320kv_armoured_dc_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:armor</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:armor</code>.</p>
</details>

```

```@raw html
<h4>320 kV unarmoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.019 (3, 1)</td><td style = "text-align: left;">2.064 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.019 (3, 1)</td><td style = "text-align: left;">2.064 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_320kv_no_armour_dc_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_320kv_no_armour_dc_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_320kv_no_armour_dc_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>.</p>
</details>

```

```@raw html
<h4>380 kV 2000 mm² cables in flat vertical formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.529 (7, 1)</td><td style = "text-align: left;">0.3983 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.529 (7, 1)</td><td style = "text-align: left;">0.3983 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_380kv_2000mm2_flatver</code>.</p><ol>
<li><code>dielectric_baseline/cable_380kv_2000mm2_flatver__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_380kv_2000mm2_flatver__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:jacket</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:jacket</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:jacket</code>.</p>
</details>

```

```@raw html
<h4>380 kV armoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">2.594 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">2.594 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_380kv_armoured_ac_flat</code>.</p><ol>
<li><code>dielectric_baseline/cable_380kv_armoured_ac_flat__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_380kv_armoured_ac_flat__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:armor</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:armor</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:armor</code>.</p>
</details>

```

```@raw html
<h4>380 kV unarmoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">2.593 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">2.593 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_380kv_no_armour_ac_flat</code>.</p><ol>
<li><code>dielectric_baseline/cable_380kv_no_armour_ac_flat__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_380kv_no_armour_ac_flat__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>, <code>5=cable:3:core</code>, <code>6=cable:3:sheath</code>.</p>
</details>

```

```@raw html
<h4>525 kV 1600 mm² armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.577 (4, 1)</td><td style = "text-align: left;">0.4093 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.577 (4, 1)</td><td style = "text-align: left;">0.4093 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_525kv_1600mm2_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_525kv_1600mm2_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_525kv_1600mm2_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:armor</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:armor</code>.</p>
</details>

```

```@raw html
<h4>525 kV unarmoured land cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">2.594 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">2.594 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_525kv_land_no_armour_ac_flat</code>.</p><ol>
<li><code>dielectric_baseline/cable_525kv_land_no_armour_ac_flat__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_525kv_land_no_armour_ac_flat__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>, <code>5=cable:3:core</code>, <code>6=cable:3:sheath</code>.</p>
</details>

```

```@raw html
<h4>525 kV unarmoured land cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.022 (3, 1)</td><td style = "text-align: left;">2.065 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.022 (3, 1)</td><td style = "text-align: left;">2.065 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_525kv_land_no_armour_dc_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_525kv_land_no_armour_dc_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_525kv_land_no_armour_dc_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:2:core</code>, <code>4=cable:2:sheath</code>.</p>
</details>

```

```@raw html
<h4>525 kV armoured subsea cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">2.594 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">2.594 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_525kv_subsea_armoured_ac_flat</code>.</p><ol>
<li><code>dielectric_baseline/cable_525kv_subsea_armoured_ac_flat__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_525kv_subsea_armoured_ac_flat__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:armor</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:armor</code>, <code>7=cable:3:core</code>, <code>8=cable:3:sheath</code>, <code>9=cable:3:armor</code>.</p>
</details>

```

```@raw html
<h4>525 kV armoured subsea cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.023 (4, 1)</td><td style = "text-align: left;">2.065 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.023 (4, 1)</td><td style = "text-align: left;">2.065 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_525kv_subsea_armoured_dc_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_525kv_subsea_armoured_dc_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_525kv_subsea_armoured_dc_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:armor</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:armor</code>.</p>
</details>

```

```@raw html
<h4>640 kV 2000 mm² cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.538 (4, 1)</td><td style = "text-align: left;">0.3985 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.538 (4, 1)</td><td style = "text-align: left;">0.3985 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>cable_640kv_2000mm2_bipole</code>.</p><ol>
<li><code>dielectric_baseline/cable_640kv_2000mm2_bipole__default__pscad_reference</code></li>
<li><code>dielectric_baseline/cable_640kv_2000mm2_bipole__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:1:sheath</code>, <code>3=cable:1:jacket</code>, <code>4=cable:2:core</code>, <code>5=cable:2:sheath</code>, <code>6=cable:2:jacket</code>.</p>
</details>

```

```@raw html
<h4>Single 1000 mm² solid conductor</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.5348 (1, 1)</td><td style = "text-align: left;">0.1398 (1, 1)</td><td style = "text-align: left;">15.84 (1, 1)</td><td style = "text-align: left;">4.166 (1, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.5348 (1, 1)</td><td style = "text-align: left;">0.1398 (1, 1)</td><td style = "text-align: left;">15.84 (1, 1)</td><td style = "text-align: left;">4.166 (1, 1)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>solid_1000mm2_single</code>.</p><ol>
<li><code>dielectric_baseline/solid_1000mm2_single__default__pscad_reference</code></li>
<li><code>dielectric_baseline/solid_1000mm2_single__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>.</p>
</details>

```

```@raw html
<h4>Two buried wires with 1 mm insulation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">9.618 (1, 1)</td><td style = "text-align: left;">2.77 (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">9.618 (1, 1)</td><td style = "text-align: left;">2.77 (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr></tbody></table></div>
<details><summary>Benchmark identities and terminal order</summary>
<p>Case: <code>two_insulated_wires</code>.</p><ol>
<li><code>dielectric_baseline/two_insulated_wires__default__pscad_reference</code></li>
<li><code>dielectric_baseline/two_insulated_wires__ametani2004__pscad_reference</code></li>
</ol>
<p>Terminal order (benchmark rows 1, 2): <code>1=cable:1:core</code>, <code>2=cable:2:core</code>.</p>
</details>

```

### Frequency slices — Z/Y

#### Band `dc`

Stored range: **0.1–102.33 Hz**, **44 samples**.

```@raw html
<h4>132 kV 630 mm² cables in flat horizontal formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">4.493 (7, 7)</td><td style = "text-align: left;">5.046 (7, 7)</td><td style = "text-align: left;">1107.0 (9, 8)</td><td style = "text-align: left;">1107.0 (9, 8)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06422 (1, 1)</td><td style = "text-align: left;">0.03417 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">4.493 (7, 7)</td><td style = "text-align: left;">5.046 (7, 7)</td><td style = "text-align: left;">2003.0 (6, 5)</td><td style = "text-align: left;">2003.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06422 (1, 1)</td><td style = "text-align: left;">0.03417 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.233 (2, 1)</td><td style = "text-align: left;">0.6276 (4, 3)</td><td style = "text-align: left;">8.864 (6, 2)</td><td style = "text-align: left;">7.421 (2, 6)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06621 (3, 3)</td><td style = "text-align: left;">0.03533 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.234 (2, 1)</td><td style = "text-align: left;">0.628 (4, 3)</td><td style = "text-align: left;">100.0 (6, 1)</td><td style = "text-align: left;">100.0 (4, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06621 (3, 3)</td><td style = "text-align: left;">0.03533 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">2.753 (7, 7)</td><td style = "text-align: left;">4.181 (7, 7)</td><td style = "text-align: left;">711.0 (6, 5)</td><td style = "text-align: left;">711.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06135 (4, 4)</td><td style = "text-align: left;">0.03337 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">2.753 (7, 7)</td><td style = "text-align: left;">4.181 (7, 7)</td><td style = "text-align: left;">1608.0 (2, 3)</td><td style = "text-align: left;">1608.0 (3, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06135 (4, 4)</td><td style = "text-align: left;">0.03337 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil — homogenized assembly</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.812 (6, 5)</td><td style = "text-align: left;">0.5 (1, 1)</td><td style = "text-align: left;">8.868 (3, 9)</td><td style = "text-align: left;">7.435 (3, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06135 (4, 4)</td><td style = "text-align: left;">0.03337 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.812 (6, 5)</td><td style = "text-align: left;">0.5 (1, 1)</td><td style = "text-align: left;">8.868 (3, 9)</td><td style = "text-align: left;">23.3 (6, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06135 (4, 4)</td><td style = "text-align: left;">0.03337 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>220 kV Milliken 2500 mm² Al / 252 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.07385 (4, 4)</td><td style = "text-align: left;">0.04461 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.07385 (4, 4)</td><td style = "text-align: left;">0.04461 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18/30 kV NA2XS2Y 630 mm² Al cables in touching trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05424 (3, 3)</td><td style = "text-align: left;">0.03026 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05424 (3, 3)</td><td style = "text-align: left;">0.03026 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05889 (1, 1)</td><td style = "text-align: left;">0.03046 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05889 (1, 1)</td><td style = "text-align: left;">0.03046 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV unarmoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05951 (1, 1)</td><td style = "text-align: left;">0.03208 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05951 (1, 1)</td><td style = "text-align: left;">0.03208 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV 2000 mm² cables in flat vertical formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06675 (7, 7)</td><td style = "text-align: left;">0.04463 (7, 7)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06675 (7, 7)</td><td style = "text-align: left;">0.04463 (7, 7)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV armoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05889 (1, 1)</td><td style = "text-align: left;">0.03046 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05889 (1, 1)</td><td style = "text-align: left;">0.03046 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV unarmoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05951 (1, 1)</td><td style = "text-align: left;">0.03208 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05951 (1, 1)</td><td style = "text-align: left;">0.03208 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV 1600 mm² armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06045 (1, 1)</td><td style = "text-align: left;">0.03078 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06045 (1, 1)</td><td style = "text-align: left;">0.03078 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05536 (1, 1)</td><td style = "text-align: left;">0.04136 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05536 (1, 1)</td><td style = "text-align: left;">0.04136 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05536 (1, 1)</td><td style = "text-align: left;">0.04136 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05536 (1, 1)</td><td style = "text-align: left;">0.04136 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.0604 (1, 1)</td><td style = "text-align: left;">0.04037 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.0604 (1, 1)</td><td style = "text-align: left;">0.04037 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.0604 (1, 1)</td><td style = "text-align: left;">0.04037 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.0604 (1, 1)</td><td style = "text-align: left;">0.04037 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>640 kV 2000 mm² cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06355 (1, 1)</td><td style = "text-align: left;">0.04288 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06355 (1, 1)</td><td style = "text-align: left;">0.04288 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>Single 1000 mm² solid conductor</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06355 (1, 1)</td><td style = "text-align: left;">0.03402 (1, 1)</td><td style = "text-align: left;">0.00429 (1, 1)</td><td style = "text-align: left;">0.001651 (1, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06355 (1, 1)</td><td style = "text-align: left;">0.03402 (1, 1)</td><td style = "text-align: left;">0.00429 (1, 1)</td><td style = "text-align: left;">0.001651 (1, 1)</td></tr></tbody></table></div>
```

```@raw html
<h4>Two buried wires with 1 mm insulation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.05052 (1, 1)</td><td style = "text-align: left;">0.06414 (1, 1)</td><td style = "text-align: left;">3.341e-5 (2, 2)</td><td style = "text-align: left;">1.305e-5 (2, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.05052 (1, 1)</td><td style = "text-align: left;">0.06414 (1, 1)</td><td style = "text-align: left;">3.351e-5 (2, 2)</td><td style = "text-align: left;">1.305e-5 (2, 2)</td></tr></tbody></table></div>
```

#### Band `harmonic`

Stored range: **53.703–2570.4 Hz**, **25 samples**.

```@raw html
<h4>132 kV 630 mm² cables in flat horizontal formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">5.142 (4, 4)</td><td style = "text-align: left;">4.862 (4, 4)</td><td style = "text-align: left;">1107.0 (9, 8)</td><td style = "text-align: left;">1107.0 (8, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02161 (1, 1)</td><td style = "text-align: left;">0.0516 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">5.142 (4, 4)</td><td style = "text-align: left;">4.862 (4, 4)</td><td style = "text-align: left;">2003.0 (6, 5)</td><td style = "text-align: left;">2003.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02161 (1, 1)</td><td style = "text-align: left;">0.0516 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">4.203 (2, 2)</td><td style = "text-align: left;">3.088 (2, 1)</td><td style = "text-align: left;">11.16 (2, 6)</td><td style = "text-align: left;">9.884 (6, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02055 (3, 3)</td><td style = "text-align: left;">0.0503 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">4.213 (2, 2)</td><td style = "text-align: left;">3.094 (2, 1)</td><td style = "text-align: left;">100.0 (4, 1)</td><td style = "text-align: left;">100.0 (5, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02055 (3, 3)</td><td style = "text-align: left;">0.0503 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">5.426 (7, 7)</td><td style = "text-align: left;">4.164 (7, 7)</td><td style = "text-align: left;">711.0 (6, 5)</td><td style = "text-align: left;">711.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02414 (4, 4)</td><td style = "text-align: left;">0.05514 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">5.426 (7, 7)</td><td style = "text-align: left;">4.164 (7, 7)</td><td style = "text-align: left;">1608.0 (3, 2)</td><td style = "text-align: left;">1608.0 (3, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02414 (4, 4)</td><td style = "text-align: left;">0.05514 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil — homogenized assembly</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">3.619 (3, 1)</td><td style = "text-align: left;">2.426 (3, 1)</td><td style = "text-align: left;">11.12 (3, 9)</td><td style = "text-align: left;">9.868 (9, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02414 (4, 4)</td><td style = "text-align: left;">0.05514 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">3.619 (3, 1)</td><td style = "text-align: left;">2.426 (3, 1)</td><td style = "text-align: left;">11.12 (3, 9)</td><td style = "text-align: left;">9.868 (9, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02414 (4, 4)</td><td style = "text-align: left;">0.05514 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>220 kV Milliken 2500 mm² Al / 252 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.01881 (4, 1)</td><td style = "text-align: left;">0.04299 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.01881 (4, 1)</td><td style = "text-align: left;">0.04299 (4, 4)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18/30 kV NA2XS2Y 630 mm² Al cables in touching trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02744 (3, 3)</td><td style = "text-align: left;">0.05632 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02744 (3, 3)</td><td style = "text-align: left;">0.05632 (3, 3)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02389 (4, 1)</td><td style = "text-align: left;">0.04327 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02389 (4, 1)</td><td style = "text-align: left;">0.04327 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV unarmoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02386 (3, 1)</td><td style = "text-align: left;">0.04291 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02386 (3, 1)</td><td style = "text-align: left;">0.04291 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV 2000 mm² cables in flat vertical formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.01866 (5, 5)</td><td style = "text-align: left;">0.03765 (7, 7)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.01866 (5, 5)</td><td style = "text-align: left;">0.03765 (7, 7)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV armoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02426 (7, 1)</td><td style = "text-align: left;">0.04327 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02426 (7, 1)</td><td style = "text-align: left;">0.04327 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV unarmoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02423 (5, 1)</td><td style = "text-align: left;">0.04291 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02423 (5, 1)</td><td style = "text-align: left;">0.04291 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV 1600 mm² armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.01785 (3, 3)</td><td style = "text-align: left;">0.03809 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.01785 (3, 3)</td><td style = "text-align: left;">0.03809 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02426 (5, 1)</td><td style = "text-align: left;">0.0305 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02426 (5, 1)</td><td style = "text-align: left;">0.0305 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02389 (3, 1)</td><td style = "text-align: left;">0.0305 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02389 (3, 1)</td><td style = "text-align: left;">0.0305 (1, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02427 (7, 1)</td><td style = "text-align: left;">0.03416 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02427 (7, 1)</td><td style = "text-align: left;">0.03416 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.0239 (4, 1)</td><td style = "text-align: left;">0.03416 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.0239 (4, 1)</td><td style = "text-align: left;">0.03416 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>640 kV 2000 mm² cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.01834 (3, 1)</td><td style = "text-align: left;">0.03546 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.01834 (3, 1)</td><td style = "text-align: left;">0.03546 (1, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>Single 1000 mm² solid conductor</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.02245 (1, 1)</td><td style = "text-align: left;">0.05267 (1, 1)</td><td style = "text-align: left;">0.08528 (1, 1)</td><td style = "text-align: left;">0.04354 (1, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.02245 (1, 1)</td><td style = "text-align: left;">0.05267 (1, 1)</td><td style = "text-align: left;">0.08528 (1, 1)</td><td style = "text-align: left;">0.04354 (1, 1)</td></tr></tbody></table></div>
```

```@raw html
<h4>Two buried wires with 1 mm insulation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.06242 (2, 1)</td><td style = "text-align: left;">0.03808 (1, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.06242 (2, 1)</td><td style = "text-align: left;">0.03808 (1, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr></tbody></table></div>
```

#### Band `narrow`

Stored range: **977.24–1.0e6 Hz**, **44 samples**.

```@raw html
<h4>132 kV 630 mm² cables in flat horizontal formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">4.482 (4, 4)</td><td style = "text-align: left;">4.532 (4, 4)</td><td style = "text-align: left;">1107.0 (9, 8)</td><td style = "text-align: left;">1107.0 (9, 8)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.501 (7, 1)</td><td style = "text-align: left;">0.5883 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">4.482 (4, 4)</td><td style = "text-align: left;">4.532 (4, 4)</td><td style = "text-align: left;">1155.0 (9, 8)</td><td style = "text-align: left;">1682.0 (9, 8)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.501 (7, 1)</td><td style = "text-align: left;">0.5883 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>CIGRE TB 880 Case 0 — 132 kV 630 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">6.757 (2, 1)</td><td style = "text-align: left;">5.436 (2, 2)</td><td style = "text-align: left;">100.0 (3, 6)</td><td style = "text-align: left;">100.0 (5, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7179 (3, 1)</td><td style = "text-align: left;">0.2769 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">12.9 (2, 1)</td><td style = "text-align: left;">7.119 (2, 1)</td><td style = "text-align: left;">100.0 (5, 3)</td><td style = "text-align: left;">100.0 (1, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7179 (3, 1)</td><td style = "text-align: left;">0.2769 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">7.837 (1, 1)</td><td style = "text-align: left;">6.756 (1, 1)</td><td style = "text-align: left;">711.0 (6, 5)</td><td style = "text-align: left;">711.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.2706 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">7.838 (1, 1)</td><td style = "text-align: left;">6.757 (1, 1)</td><td style = "text-align: left;">749.7 (6, 5)</td><td style = "text-align: left;">1285.0 (6, 5)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.2706 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18 kV 1000 mm² cables in trefoil — homogenized assembly</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">fem [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">6.54 (2, 1)</td><td style = "text-align: left;">5.15 (2, 2)</td><td style = "text-align: left;">100.0 (1, 3)</td><td style = "text-align: left;">100.0 (9, 7)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.2706 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">3</td><td style = "text-align: left;">fem [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">6.54 (2, 1)</td><td style = "text-align: left;">5.15 (2, 2)</td><td style = "text-align: left;">100.0 (9, 7)</td><td style = "text-align: left;">100.0 (5, 9)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">4</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.7017 (4, 1)</td><td style = "text-align: left;">0.2706 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>220 kV Milliken 2500 mm² Al / 252 mm² Cu cables in trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.8387 (4, 1)</td><td style = "text-align: left;">0.3237 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.8387 (4, 1)</td><td style = "text-align: left;">0.3237 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>18/30 kV NA2XS2Y 630 mm² Al cables in touching trefoil</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.6601 (3, 1)</td><td style = "text-align: left;">0.2547 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.6601 (3, 1)</td><td style = "text-align: left;">0.2547 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.022 (4, 1)</td><td style = "text-align: left;">3.129 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.022 (4, 1)</td><td style = "text-align: left;">3.129 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>320 kV unarmoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.019 (3, 1)</td><td style = "text-align: left;">3.127 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.019 (3, 1)</td><td style = "text-align: left;">3.127 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV 2000 mm² cables in flat vertical formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.529 (7, 1)</td><td style = "text-align: left;">0.6032 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.529 (7, 1)</td><td style = "text-align: left;">0.6032 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV armoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">3.93 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">3.93 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>380 kV unarmoured cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">3.929 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">3.929 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV 1600 mm² armoured cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.577 (4, 1)</td><td style = "text-align: left;">0.6199 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.577 (4, 1)</td><td style = "text-align: left;">0.6199 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">3.93 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (5, 1)</td><td style = "text-align: left;">3.93 (5, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV unarmoured land cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.022 (3, 1)</td><td style = "text-align: left;">3.129 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.022 (3, 1)</td><td style = "text-align: left;">3.129 (3, 1)</td><td style = "text-align: left;">Inf (4, 2)</td><td style = "text-align: left;">Inf (4, 2)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cables in AC flat formation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">3.93 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">10.04 (7, 1)</td><td style = "text-align: left;">3.93 (7, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>525 kV armoured subsea cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">8.023 (4, 1)</td><td style = "text-align: left;">3.129 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">8.023 (4, 1)</td><td style = "text-align: left;">3.129 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>640 kV 2000 mm² cable DC bipole</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">1.538 (4, 1)</td><td style = "text-align: left;">0.6036 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">1.538 (4, 1)</td><td style = "text-align: left;">0.6036 (4, 1)</td><td style = "text-align: left;">Inf (6, 3)</td><td style = "text-align: left;">Inf (6, 3)</td></tr></tbody></table></div>
```

```@raw html
<h4>Single 1000 mm² solid conductor</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">0.5348 (1, 1)</td><td style = "text-align: left;">0.2065 (1, 1)</td><td style = "text-align: left;">15.84 (1, 1)</td><td style = "text-align: left;">6.312 (1, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">0.5348 (1, 1)</td><td style = "text-align: left;">0.2065 (1, 1)</td><td style = "text-align: left;">15.84 (1, 1)</td><td style = "text-align: left;">6.312 (1, 1)</td></tr></tbody></table></div>
```

```@raw html
<h4>Two buried wires with 1 mm insulation</h4>
<div class = "data-frame" style = "overflow-x: scroll;"><table class = "data-frame" style = "margin-bottom: 6px;"><thead><tr class = "columnLabelRow"><th class = "stubheadLabel" style = "font-weight: bold; text-align: right;">Row</th><th style = "text-align: left;">reference</th><th style = "text-align: left;">candidate</th><th style = "text-align: left;">Z NRMSE</th><th style = "text-align: left;">Z pointwise</th><th style = "text-align: left;">Y NRMSE</th><th style = "text-align: left;">Y pointwise</th></tr></thead><tbody><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">1</td><td style = "text-align: left;">pscad [1]</td><td style = "text-align: left;">coaxial [1]</td><td style = "text-align: left;">9.618 (1, 1)</td><td style = "text-align: left;">4.197 (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr><tr class = "dataRow"><td class = "rowLabel" style = "font-weight: bold; text-align: right;">2</td><td style = "text-align: left;">pscad [2]</td><td style = "text-align: left;">coaxial [2]</td><td style = "text-align: left;">9.618 (1, 1)</td><td style = "text-align: left;">4.197 (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td><td style = "text-align: left;">Inf (2, 1)</td></tr></tbody></table></div>
```

#### Band `wide`

No stored samples or comparisons for this band; errors are missing.

Independent cold/warmed timings: not recorded. Stored elapsed-at-completion values describe calculation batches, not per-selection execution.



## Comparison conventions

No denominator floor or solver-data modification is applied. Absolute RMS
retains the measured difference. If the reference trace lies within the recorded
observable tolerance, relative RMS is unavailable with a per-cell reason.
Pointwise normalization is also unavailable when any selected reference sample
is numerically zero; no samples are omitted. Defaults are 1e-10 Ω/m for R,
1e-15 H/m for L, 1e-12 S/m for G and 1e-16 F/m for C. Z and Y thresholds follow
R + 2πfL and G + 2πfC. Override them with `compare(...; atol=(G=..., C=...))`.
Unsupported observables and empty bands retain their separate reasons.

Saved comparisons retain the policy used when they were calculated. Historical
tables can therefore contain zero or infinite ratios under the previous policy;
recompute comparisons from their retained raw operands to apply the current
policy. The report renderer does not change stored metrics.

The deterministic summary displays **Z and Y only**. Shunt conductance
`G = real(Y)` remains available for an explicitly requested loss study through
`compare(reference, candidate, G; ...)`; it is not an additional default KPI.
Saved G comparisons are retained, but do not expand the Z/Y summary.
UQ moment benchmarks, when selected, have separate mean and standard-deviation
sections; their R/L/C/G observables are not mixed into deterministic Z/Y tables.

Two profiles are included: all-`:default`, and `:Ametani2004` for both
insulation and semicon with all other slots at `:default`. Defaults
are contextual and may implement different approximations in each backend.
PSCAD's scalar export matches the selected radial dielectric admittance at
its base frequency (50 Hz in this campaign), before the native loss-tangent
cap of 10. That fit is not an arbitrary broadband constitutive law. The owned
engine and FEM evaluate the retained dielectric constituents at each frequency;
homogeneous radial equivalence does not assert full-geometry field equivalence.
In this campaign, each benchmark explicitly selects PSCAD or FEM as reference
and coaxial as candidate. This is a choice in its definition, not a Gauntlet
rule: same-backend and other cross-backend pairs are equally valid.
No backend is ground truth or
an approved CI reference. Inputs, terminal order, basis and frequencies must
match; there is no implicit interpolation or conversion. No detailed HTML
pages, plots, input-object dumps or numerical-file copies are published here.

## Run and select stored results

```bash
lcm gauntlet run --directory /path/to/campaign \
  --backends coaxial,fem,pscad --formulas default --frequency-range 0.1,1e6
lcm gauntlet status --directory /path/to/campaign
lcm gauntlet resume --directory /path/to/campaign
lcm gauntlet run --directory /path/to/lossy-campaign \
  --backends coaxial,fem,pscad --formulas default --dielectric Ametani2004 \
  --frequency-range 0.1,1e6

lcm gauntlet compare --definition /path/to/benchmarks.toml \
  --output /path/to/comparisons
LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/comparisons \
  julia --project=docs docs/make.jl
```

Omit `--cases` to run the indexed catalogue, or pass comma-separated case IDs.
`--frequency-range` sets the same 101 logarithmic samples for every case;
without it, each case retains its own extent, with a 0.1 Hz minimum.
`--formulas catalogue` retains the broader formula sweep without expanding
this summary. Formula grids, explicit lossy selections and uncertainty runs
are documented in the
[Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/gauntlet/README.md).
FEM Monte Carlo remains disabled.

A benchmark definition names exactly one reference and one candidate artifact
(paths and SHA-256 checksums), plus quantities, bands and normalizations. See the
CLI guide for the TOML format. Comparing saved files never reruns their solvers.
Missing operands are errors; they never trigger a replacement reference.
Multiple comparison directories can be selected with the platform path-list
separator (`:` on Unix, `;` on Windows); identical benchmark records appear once.
This completed-only campaign retains 86 calculations and 46 benchmark pairs:
38 PSCAD-reference and 8 FEM-reference. The 32 unfinished FEM selections remain
excluded. No simulations are resumed by this page. Batch elapsed-at-completion
values are not advertised as per-selection cold or warmed execution timings.

## Numerical references for CI

Stored Gauntlet results are not automatically approved CI references. The
separate [numerical-reference gate](https://github.com/Electa-Git/LineCableModels.jl/tree/main/test/numerical)
requires reviewed results, explicit tolerances and artifact bindings. It never
starts a Gauntlet campaign or refreshes a reference. Approval remains separate
from collecting results or displaying this summary.

## Benchmark data API

```@docs
LineCableModels.Engine.RMSError
LineCableModels.Engine.LineParametersBenchmark
LineCableModels.Engine.compare
LineCableModels.Engine.absolute_error
LineCableModels.Engine.relative_error
```
