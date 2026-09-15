# Registered frequency-dependent soil relations

The frequency-dependent earth-property registry contains the following
explicit relations. Their names identify the principal source contributor;
the source references remain part of each implementation's documentation.

| Formula identifier | Compact description | Long description |
| --- | --- | --- |
| `:constant` | Constant | Frequency-independent earth-material pass-through |
| `:alipio2014` | Alipio | Alipio–Visacro causal soil dispersion (2014) |
| `:cigre2019` | CIGRE | CIGRE WG C4.33 recommended soil dispersion (2019) |
| `:datsios2019` | Datsios | Datsios–Mikropoulos two-limit soil fit (2019) |
| `:longmire1975` | Longmire | Longmire–Smith 13-term dielectric relaxation (1975) |
| `:messier1985` | Messier | Messier square-root soil dispersion (1985) |
| `:portela1999` | Portela | Portela power-law soil dispersion (1999) |
| `:scott1967` | Scott | Scott–Carroll–Cunningham empirical moist-soil fit (1967) |
| `:visacro1987` | Visacro | Visacro–Portela empirical soil dispersion (1987) |
| `:visacro2012` | Visacro | Visacro–Alipio empirical soil dispersion (2012) |

`:default` is a routing alias for `:constant` and is not a second soil law.
The concrete formulas accept their source-specific physical parameters through
the `parameters` field of the formula selection.

[Back to Earth properties](../../earth_properties.md)
