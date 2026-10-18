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

After `earth`, `representation` compares `conductor`, `semicons`, `armour` (armoured
cases), `dielectric` and `shunt` with the baseline, saves `representation.json`, and
ranks scaled deviations. The arguments `stranded_core` (default `true`) and
`core_diameter` (46.3 mm) select the 1600 mm² core, with WP1 corrected materials for
its solid variant. The argument `hollow_core` (default `false`) selects annular land
or subsea cores of metal area `conductor_cross_section`, with resistivity set to keep
the 20 °C DC resistance. The argument `armour_as_tube` (default `false`) selects a
steel annulus of equal metal area after the declared bedding, without the wire fit
rule. Layers of zero thickness are absent, and `semicons` folds inner semicon, outer
semicon and water blocking into insulation at its baseline permittivity.
The dielectric and shunt variants apply `insulation_admittance = :lossless` with
`semicon_admittance = :lossless`, and `shunt_model = :boundary`, respectively, to the
baseline problem.

The final `uncertainty` section gives the `swept` inputs independent uniform ±10 % distributions
(`10 / sqrt(3)` % standard uncertainty), with other inputs unchanged.
Monte Carlo uses 200 successful trials and seed 20261016, with retained diagnostics
and at most 20 failures, alongside linear error propagation (LEP).
The KPIs scale the largest entry standard deviation (s) or absolute Monte Carlo
mean shift (m) by the largest absolute baseline diagonal entry.
For Z and Y, outputs give band RMS (all, DC, harmonic, narrow), maxima with frequency,
and LEP/MC ratios of all-band RMS spread.
The `ranking.txt` file records KPIs, trials, failures and wall times, and
`uncertainty.json` retains both methods' observations.
The section dominates a script's run time.

`sensitivity_tables.jl` reads the newest complete `ranking.txt` of each script and
writes LaTeX tables (longtable, booktabs, siunitx) without computing any case:

```sh
julia benchmarks/sensitivity_tables.jl [output folder]
```

The default output folder is `.linecablemodels/benchmarks/tables/`.

`julia --project=test benchmarks/cable_320kv_no_armour_dc_bipole.jl wp22` selects the DC320-U WP2.2 study after the baseline.
DC320-A uses the same `wp22` and optional `pscad` arguments, with separate sheath and armour bonds and a pulse window from its modal delays.
No argument selects WP2.1.
The 48 ULM variants vary Yc and H orders (8, 16), time rows, passivity enforcement and DC correction.
The study also fits default FDCM and Marti models.

Three NLT step references use single-point sheath bonds of 1, 10 and 20 Ω.
Every baseline fit runs at 10 Ω for the single-point and both-ends step scenarios.
The third scenario uses a double exponential.
The 1.2/50 µs pulse has a 320 kV peak, single-point bonds and an open receiving end.
Its 3 ms window has 12000 steps of 0.25 µs and includes the first modal arrival.
The pulse reference uses N = 8192 and T = 3.75 ms, with a highest frequency of 2.1844 MHz.

A second baseline uses 101 logarithmic frequencies from 1 Hz to 1 MHz plus 50 Hz.
The ULM 16/16 fits use time rows, passivity enforcement and each DC choice, alongside default Marti.
The study compares their fit errors and single-point step deviations with the 227-frequency fits.

The run folder holds `wp22_segment.json`, `wp22_fits.json`, `wp22_references.json`, `wp22_response_<index>.json`, `wp22_ii_response_<index>.json`, `wp22_iii_response_<index>.json`, `wp22_sampling_baseline.json`, `wp22_sampling.json` and `wp22_fit_<index>_<quantity>_1_fitting_error.xlsx` workbooks.
Archives keep fit diagnostics, DC targets and achieved values, core 1 KPIs, wall times and transient RMS deviations scaled by the reference peak.
Growth ratios compare peaks in the final two 10 % windows.
A growth flag identifies a ratio that exceeds one and 110 % of its reference ratio.

Append `pscad` to run the optional PSCAD section after the other studies: `julia --project=test benchmarks/cable_320kv_no_armour_dc_bipole.jl wp22 pscad`.
Set `LINECABLEMODELS_PSCAD_CONFIG` to the station configuration path.
The section runs FD-Phase orders 10, 20 and 30 with each DC and passivity setting, and default FD-Mode.
It covers the single-point step at 10 Ω and the pulse.
It saves station identity, responses, comparisons, statuses with failure messages, and native projects, logs and outputs under `wp22_pscad_<scenario>_<index>/` in the run folder.
A station failure produces a status file, and execution continues.
