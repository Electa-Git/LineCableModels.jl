# Benchmark cable cases

These scripts declare the benchmark-system cable cases for the wideband model
assessment using the construction API of LineCableModels. Values follow the
benchmark specification, Appendix A, with the selected core representation below.

| Script | Specification case |
| --- | --- |
| `cable_320kv_armoured_dc_bipole.jl` | DC320-A |
| `cable_320kv_no_armour_dc_bipole.jl` | DC320-U |
| `cable_380kv_armoured_ac_flat.jl` | AC380-A |
| `cable_380kv_no_armour_ac_flat.jl` | AC380-U |
| `cable_525kv_land_no_armour_dc_bipole.jl` | DC525-L |
| `cable_525kv_land_no_armour_ac_flat.jl` | AC525-L |
| `cable_525kv_subsea_armoured_dc_bipole.jl` | DC525-S |
| `cable_525kv_subsea_armoured_ac_flat.jl` | AC525-S |

Run a script from the repository root:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl
```

Use the same command with any script in the table. Each baseline computes 227
frequencies from 1 mHz to 1 MHz, including 50 Hz. Outputs go to
`.linecablemodels/benchmarks/<script name>/<yyyymmdd-HHMMSS>/`: `problem.json`
contains the physical problem. The `baseline.json` file retains the full Z/Y
matrices as R/X and G/B components, without clipping, in units per metre.

The 1600 mm² baseline uses an explicit stranded core with uncorrected materials.
The package derives the corrections. Land and subsea cases use solid equivalent
cores. Temperatures use °C. Other builder inputs use SI units. The builders
convert catalog voltage and area to kV and mm² for the design metadata.

These scripts are manual benchmark cases, separate from the CI fixtures.

After the baseline, each script runs three sensitivity sections. The
`one_at_a_time` section computes each ±10 % level independently for every diameter,
thickness, material resistivity, relative permeability, relative permittivity,
temperature coefficient and lay ratio. Counts, positions, line length, voltage,
conductor cross section, reference temperature and frequencies are unchanged.
The `temperature` section uses one Gridspace with 20, 40 and 90 °C against the
70 °C baseline. The `earth` section uses one with 100 and 10000 Ω·m against 1000 Ω·m.

Each output folder also contains `one_at_a_time.json`, `temperature.json` and
`earth.json`, which keep observations and native comparisons in the all, DC,
harmonic and narrow bands. Armour bedding uses the larger of its declared thickness
and the thickness required to fit the wires, and moves all outer layers with it.

The scaled deviation e at each frequency is the largest absolute matrix-entry difference
divided by the largest absolute baseline diagonal entry.
The scripts print rankings and save `ranking.txt`, ordered by the all-band RMS of e
for Y, then Z, largest first.
Each quantity also has the RMS of e in each native band and its maximum e with frequency.
