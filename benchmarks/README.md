# Benchmark cable cases

These scripts declare the benchmark-system cable cases for the wideband model
assessment using the construction API of LineCableModels, with the selected core
representation below.

| Script | Physical description |
| --- | --- |
| `cable_320kv_armoured_dc_bipole.jl` | 320 kV DC, armoured, land, bipole, stranded core |
| `cable_320kv_no_armour_dc_bipole.jl` | 320 kV DC, unarmoured, land, bipole, stranded core |
| `cable_380kv_armoured_ac_flat.jl` | 380 kV AC, armoured, land, flat arrangement, stranded core |
| `cable_380kv_no_armour_ac_flat.jl` | 380 kV AC, unarmoured, land, flat arrangement, stranded core |
| `cable_525kv_land_no_armour_dc_bipole.jl` | 525 kV DC, land, unarmoured, bipole, solid equivalent core |
| `cable_525kv_land_no_armour_ac_flat.jl` | 525 kV AC, land, unarmoured, flat arrangement, solid equivalent core |
| `cable_525kv_subsea_armoured_dc_bipole.jl` | 525 kV DC, subsea, armoured, bipole, solid equivalent core |
| `cable_525kv_subsea_armoured_ac_flat.jl` | 525 kV AC, subsea, armoured, flat arrangement, solid equivalent core |

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
`core_diameter` (46.3 mm) select the 1600 mm² core, with corrected materials for
its solid variant. The argument `hollow_core` (default `false`) selects annular land
or subsea cores of metal area `conductor_cross_section`, with resistivity set to keep
the 20 °C DC resistance. The argument `armour_as_tube` (default `false`) selects a
steel annulus of equal metal area after the declared bedding, without the wire fit
rule. Layers of zero thickness are absent, and `semicons` folds inner semicon, outer
semicon and water blocking into insulation at its baseline permittivity.
The dielectric and shunt variants apply `insulation_admittance = :lossless` with
`semicon_admittance = :lossless`, and `shunt_model = :boundary`, respectively, to the
baseline problem.

Without arguments, the original joint uncertainty section follows the sensitivity
comparisons. It uses 200 uniform trials at ±10 % and Measurements-based LEP.
Its `uncertainty.json` and `ranking.txt` keep the original matrix comparisons.

Dedicated uncertainty studies cover the 320 kV DC bipoles with and without armour.
They also cover the 525 kV subsea armoured DC bipole and 380 kV armoured AC line.

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl uncertainty
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl uncertainty reduced
```

The full study uses 2000 Monte Carlo trials at ±10 % and 500 at ±5 % and ±20 %.
Each swept input has independent uniform uncertainty, declared with
`Grid(nominal, (percent / sqrt(3),))`. The seed is 20261016.
Monte Carlo composes each trial with `modal = :default`.

Both LEP methods compute phase Z/Y directly.
Central LEP also runs the modal composition for the modal observables.
Its timing includes both computations.

Central perturbations use the package default of 0.01 times each input standard uncertainty.
Near a radius transition, a central perturbation can select different branches
of the builder. Measurements follows the branch selected at the nominal point.
These methods can give different spreads near the wire-fit limit.

Each method runs once for compilation and once for warm timing.
The first timing includes computation. Estimated compilation overhead equals the
first timing minus warm timing, with zero as the lower bound.

The solves use the declared 227 frequencies. Trial archives keep complete ordered
coefficients at 1 mHz, 1 Hz, 50 Hz, 1 kHz, 10 kHz, 100 kHz and 1 MHz.
Modal tables use the last six frequencies.

For Monte Carlo, the generator reconstructs phase Z/Y from raw modal R/L/C/G coefficients and paired bases.
The LEP archives supply phase Z/Y from the direct computation.
The archives also keep trial attenuation, velocity and surge impedance.
The generator derives impedance magnitudes and phase travel time for each trial.

Mode assignment uses the nominal reference at each frequency.
The package matches numerically degenerate reference modes through their cluster spans
and rotates each assigned basis onto the reference basis.
Assignment overlaps and permutations have axes `(mode, frequency, trial)` or `(mode, frequency, evaluation)`.
Reference cluster sizes have axes `(mode, frequency)`.
The archives keep these records and their full frequency axes.

The diagnostic table gives accepted counts and the minimum overlap across all modes and frequencies.
It also gives the count and values of clustered frequencies.
Archives that predate cluster records have empty cluster columns.
The package conditions modal results on successful mode assignment.

Armoured cases compare the original nominal with a second bedding nominal under
the ordinary maximum of stack and wire-fit radii. Each comparison uses 1000 trials
at ±10 % and both LEP methods for phase Z/Y. Central LEP also computes modal
observables. Each nominal has its own baseline and modal reference.

The CSV records `nominal_id` and `bedding_thickness_m` for each result.
Each mean shift uses its own nominal value as the reference.
The LEP/MC ratio compares spreads from the same nominal.

The second design uses 1.05 times the minimum nominal thickness needed to clear
the wire-fit rule throughout the full ±10 % input support. The minimum uses
the smallest stack dimensions and `bedding_thickness`, with the largest wire diameter.
An increase of 5 % over the fit thickness at nominal dimensions does not clear this support.
The geometry table includes that literal calculation as a diagnostic.

The geometry archive records the nominal clearance above the wire-fit radius.
It also records minimum clearances over the inner-input and full ±10 % supports.
Those bounds include the tolerance of `bedding_thickness`. The generator writes these margins to
`uncertainty_kink_geometry.csv` and `.tex` and the summary.
A positive minimum means the builder uses the stack radius throughout that support.

The inner-input fixed-radius comparison remains a partial control.
The builder keyword `armour_wire_centre_radius` defaults to `nothing` for the
ordinary rule. A supplied radius determines the bedding thickness and fixes
the armour radius. Infeasible geometry produces a status archive, and the study
continues with the declared input distribution.

The inner group includes bedding and the material properties of layers under the armour.
Shared polypropylene properties belong to this group and also affect the jacket.
The outer group contains armour properties, wire diameter, lay ratio and jacket thickness.
Both groups also keep ordinary-rule comparisons at the original nominal.

Run only the kink item, with 100 trials per Monte Carlo population, with:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl uncertainty kink reduced
```

Without `reduced`, this item uses 1000 trials per population. It skips the other
tolerance levels and transient computations. The full `uncertainty` command
includes these revised kink controls.

The three DC cases compute voltage transient KPIs at ±10 % with 500 trials and central LEP.
A core step uses 10 Ω bonds for sheaths and armour at the source end.
The remote end is open. NLT uses N = 4096 and Tfactor = 1.25.
The nominal window is twenty times the fastest modal delay, with 4000 time steps.
These settings match the prior step study.

The transient archive keeps the remote core peak and its time.
It also keeps the 10 %-of-peak arrival time and sheath and armour peaks at both ends.
The sampled time grid can produce zero or discontinuous timing derivatives.
These results remain visible in the tables.

The `reduced` argument uses 200 trials at ±10 % and 100 for each other comparison.
Physical and solver settings remain the same. The package `save` API writes JSON.
`uncertainty_run.json` marks a complete run and records input sets and wall times.

The `frequency_domain` selector runs every configuration except transient Monte Carlo.
The DC cases keep the transient baseline and central LEP.
Each invocation writes a new run folder.

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl uncertainty frequency_domain
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl uncertainty frequency_domain reduced
```

Generate CSV and LaTeX tables from complete run folders:

```sh
julia --project=test benchmarks/uncertainty_tables.jl <run folder>... <output folder>
```

To combine a frequency-domain run with archived transient Monte Carlo for the same cable, use:

```sh
julia --project=test benchmarks/uncertainty_tables.jl <frequency-domain folder> <earlier full folder> <output folder>
```

Supply at most one frequency-domain folder per cable.
The generator takes transient trials from the first supplied folder for that cable that contains the Monte Carlo archive.
The frequency-domain folder supplies the transient baseline and central LEP when present.
The generator requires equal nominal inputs in the two source declarations.
The combined output uses the frequency-domain run folder name.
Other complete runs keep separate output directories.

Rows identify their Monte Carlo, LEP and baseline source folders.
The text summary and `uncertainty_sources.csv` list the source folders for each table.

Nested trial prefixes use N = 100, 200,
500, 1000 and 2000, up to the accepted count, plus the complete accepted sample.
Each prefix uses 2000 percentile bootstrap resamples with seed 20261018 plus N.
The same trial resample serves every observable.

Tables report signed mean shifts, spread and 95 % intervals in physical units.
Phase coefficients have real, imaginary and magnitude rows.
Modal rows give attenuation in dB/km, velocity in m/µs, surge impedance in Ω and travel time in µs.
For each ratio, bounds equal LEP spread divided by the opposite bounds of the Monte Carlo interval.

Each comparison has a separate class for spread and mean.
Spread is `equivalent` when the whole ratio interval lies within [0.9, 1.1].
It is `different` when the whole interval lies below 0.9 or above 1.1.
An interval that overlaps a bound is `inconclusive`.
Spread is `insensitive` when Monte Carlo sigma is below 0.1 % of the absolute nominal value or exactly zero.
Unsupported and nonfinite LEP quantities have separate classes.

The mean has a `mean shift` flag when its interval excludes zero and the absolute
point shift exceeds 1 % of the absolute nominal value.
A `within 1%` mean has its whole shift interval within that tolerance.
Other means are `inconclusive`.

Zero nominal values keep absolute shifts. Intervals are marginal bootstrap estimates.
The CSV keeps the original combined `agrees` flag.
The summaries use the separate classes.

The structural-zero threshold is one billionth of the largest complex entry in each Z or Y matrix.
The test compares magnitudes at the same frequency.
It applies to the complex entry for all its component rows.
It preserves small loss components of a substantial complex entry.
The generator lists structural zeros and excludes them from comparison counts.

Output prefixes `uncertainty`, `uncertainty_kink` and `uncertainty_transient` each
have CSV and LaTeX tables. Separate CSV files hold timing and diagnostic data.
The class columns identify self and mutual entries within each cable by conductor pair.
They also identify entries between cables and modal quantities by mode.
Transient classes identify the terminal and line end.

The files `uncertainty_classes` and `uncertainty_observable_classes` give CSV and
LaTeX summaries by physical class, study, nominal and method at the complete sample.
The latter also separates observables.
Each summary gives spread counts, mean counts, excluded counts and the warm speed-up.
Counts include directed entries, cable copies and observable components.
These counts describe table rows, which can share the same input draws.

The central-to-Measurements ratio in `uncertainty_lep_methods.csv` marks a method difference outside [0.9, 1.1].
Its support columns identify which method has spread equivalent to Monte Carlo.
Class totals appear in `uncertainty_lep_method_classes.csv` and `.tex`.
The file `uncertainty_structural_zeros.csv` lists excluded rows at the complete sample.

The files `uncertainty_kink_effects.csv` and `.tex` give paired changes in Monte
Carlo means, defined as fixed radius minus ordinary rule. Paired comparisons use
identical input sets and complete populations with zero rejected trials.

A text summary records the source, wall times and estimated full computation cost.
Speed-up compares the measured Monte Carlo population with the measured LEP pass.
A separate estimate scales the time for Monte Carlo to each nested prefix.

Full-run estimates scale the sample cost with trial count.
LEP evaluation counts remain constant. CSV files keep numeric precision.
LaTeX uses three significant digits with longtable and booktabs.

`sensitivity_tables.jl` reads the newest complete `ranking.txt` of each script and
writes LaTeX tables (longtable, booktabs, siunitx) without computing any case:

```sh
julia benchmarks/sensitivity_tables.jl [output folder]
```

The default output folder is `.linecablemodels/benchmarks/tables/`.

All eight cable scripts also accept `modal` alone:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl modal
```

This argument measures how cable construction, temperature, soil and layout
change the propagation modes. The scans use a line-parameter Gridspace and
`compute(space, formulation; modal = :default)` over the same 227 frequencies.
The case script saves native observations directly with the package API.
Its run folder contains `modal_run.json` with the declared inputs and wall times,
plus `modal_<study>_<index>.json.gz` for each computed point.
The package writes JSON, and gzip compresses each point archive.
These archives keep full propagation constant, velocity and characteristic impedance
curves at all 227 frequencies. Both modal bases use the six report frequencies.
The phase Z/Y coefficients use only the lowest frequency for the DC estimates.
The archives also keep the modal diagnostics.

Run the separate generator with one or more run folders, then an output folder:

```sh
julia --project=test benchmarks/modal_tables.jl <run folder>... <output folder>
```

The generator reads the archives and derives the tables without a physical computation.
It writes each run under `<output folder>/<script name>/<timestamp>/`.
The output prefixes describe these physical comparisons:

| Output prefix | Physical comparison |
| --- | --- |
| `modal_baseline` | Current patterns, attenuation, phase velocity, characteristic impedance and phase travel time for the declared cable. |
| `modal_one_at_a_time` | Each swept input at nine levels from 90 % to 110 % of nominal, with changes in physical units and percent. |
| `modal_temperature_earth` | The product of temperature and soil resistivity. Temperature spans 20 to 90 °C at 10 °C intervals. Resistivity spans 10 to 10000 Ω·m at quarter-decade intervals. |
| `modal_layout` | Cable separation in metres: 0.25, 0.5, 0.75, 1, 1.5 and 2. Depth in metres: 0.8, 1, 1.25, 1.5 and 2. |
| `modal_corrections` | Solid core, core lay, folded dielectric layers, armour tube, armour lay, PE permittivity and the complete corrected set, each against the declared baseline. |
| `modal_earth_formulas` | Supported pairs of buried-conductor earth impedance and admittance formulas against the Unified pair. The pairs include ideal earth admittance. |

Each prefix has a JSON archive compressed with gzip. It holds full propagation
constant, velocity and characteristic impedance curves from native modal observations.
The CSV table gives attenuation in dB/km, velocity in m/µs, impedance magnitude
in Ω and phase travel time in µs. It samples 1 Hz, 50 Hz, 1 kHz, 10 kHz, 100 kHz
and 1 MHz. Divide the line length by the phase velocity to calculate the phase travel time.

The CSV also gives the sheath and armour skin transition frequencies.
These use the material resistivity at the selected temperature and permeability.
For armour wires, the skin length is wire diameter. For a tube, it is wall thickness.
Each `_conductors.csv` file records the conductor resistance diagonal in Ω/km
and the complete shunt conductance matrix in S/km at the lowest frequency.
These phase coefficients give an estimate of DC resistance.

For each input, the nonlinearity fraction measures the largest departure of its
nine-level observable curve from the line through its endpoints.
The fraction uses the absolute change between endpoints as its denominator.
A zero change with nonzero departure gives NaN. A percentage with a zero
baseline has an empty field.

Mode classes follow the current amplitudes and relative phases at each frequency.
Opposite core and sheath currents identify the coaxial family.
Opposite sheath and armour currents identify the sheath-armour family.
Net current identifies the dominant return conductor. Cable current correlation
distinguishes common currents from currents that flow against another cable.

For three cables, the same rule aggregates terminal amplitudes across cables.
It tests the cable pairs whose current norms exceed 20 % of the total norm.
The generated JSON archives keep these classes at the six report frequencies
and keep the study inputs.
Mode assignment uses maximum normalized current-basis overlap at 10 kHz with
Hungarian assignment. Tables keep the baseline mode, assigned overlap and
largest two candidate overlaps.

The `modal_summary.txt` file lists ambiguous matches, infeasible layouts,
correction availability and wall times for each study.
A match is ambiguous when its best or assigned overlap is below 0.9.
It is also ambiguous when its two largest overlaps differ by at most 0.05.
Mode numbers keep the package order.

The `_status.csv` file records corrections that are not applicable or whose
corrected material values are absent from the case representation.
The `_registry.csv` file lists each registered earth formula and its availability.
Each case script declares six pairs. Unified, Saad 1996 and Wedepohl-Wilcox 1973
impedance each use Unified or ideal admittance.
The registry metadata retains the availability recorded for these comparisons.

The `modal_fallback_continuity.csv` file compares velocity and attenuation jumps
on each side of a fallback with the same numbered baseline mode.
A relative jump equals the absolute change divided by the absolute previous value.
The reference is the largest adjacent jump within five baseline samples centred
on the fallback index. The table flags a jump above three times that reference.

The `modal_mode_order.csv` file records mode assignments from the current bases
at 1 Hz, 50 Hz and 1 kHz. It flags a changed mode number or overlap below 0.9.
The summary records these checks for each point.

The three DC study scripts select frequency-domain sensitivity when run without arguments.
They also accept `transient_sensitivity` alone:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl transient_sensitivity
```

This study measures the transient consequence of physical input variations with
the exact NLT reference. It uses the reference core step with single-point bonds
at 10 Ω and the impulse. Both keep their nominal windows and NLT settings.
Each line-problem Gridspace lifts into a transient problem through the package API.

Before each transient calculation, a line-parameter Gridspace checks Z and Y on
the real frequency axis with the same formulation. The band limits come from
the scenario's NLT settings: `1 / (2T)` to `(2N - 1) / (2T)`, with
`T = Tfactor * tmax`. The log grid has at least 400 intervals per decade and
includes both limits.

The check uses the smallest eigenvalue of `(A + Aᴴ) / 2`
for each matrix. Any strictly negative eigenvalue marks the point
`not_computable`. Its record gives the first negative sample, the previous
sample and the offending quantity. A point without a negative sample records
`passive` for that band.

The `*_passivity.json` archives keep both eigenvalue curves and the frequency grid.
Flagged points keep the passivity evidence and status. They have no transient
calculation, waveform observables or waveform archive.

Every swept input uses 90, 100 and 110 % of nominal. Separate grids use temperatures
of 20, 40, 70 and 90 °C and earth resistivities of 100, 1000 and 10000 Ω·m.
The correction steps and their availability match the modal study.

The `transient_sensitivity.json` archive keeps core 1 peak voltage, peak time and
arrival time at the `receiving` end. Arrival uses 10 % of the peak.
It also keeps voltage peaks for every sheath and armour terminal at both ends.
The core 1 current peak uses the `sending` end.
The RMS deviation of core 1 voltage at the `receiving` end uses its nominal peak
as the denominator.

Voltages use kV, times use µs and currents use kA. Peaks are absolute magnitudes.
The RMS ratio is dimensionless. The package writes these observations with `save`.
Each point also records its window, NLT controls and passivity summary.
The `transient_sensitivity_run.json` archive records inputs, correction availability,
retained waveform filenames and wall times.

For each scenario, full voltage and current waveforms remain for the nominal point.
They also remain for both extreme levels of the five inputs with the largest peak
change at the `receiving` end.
This selection considers computable levels of the swept inputs, temperature and earth resistivity.
An extreme level without a computable response has no waveform archive.
Other points keep scalar observations.

Generate CSV and LaTeX tables from one or more completed run folders:

```sh
julia --project=test benchmarks/transient_sensitivity_tables.jl <run folder>... <output folder>
```

The generator writes tables under `<output folder>/<script name>/<timestamp>/`.
Each scenario has changes in physical units and percent relative to nominal.
The ±10 % inputs follow their largest absolute peak change at the `receiving` end,
largest first.
Temperature, earth and correction rows follow. A separate table states which
corrections are available. Percent changes with a zero nominal value have an empty CSV field
and a dash in LaTeX.

Flagged rows give the status, onset frequency and quantity with a negative eigenvalue.
Numerical observables and changes remain empty for these points.
The `transient_passivity.csv` table lists the status and passivity summary of every point.
The generator needs archives with passivity records.

`time_domain` runs line-model fits and their time-domain comparisons without PSCAD.
`pscad` runs only the baseline, the NLT references for the single-point step and impulse scenarios at 10 Ω, and PSCAD.
`time_domain pscad` runs both time-domain comparisons and PSCAD.
`model_settings` alone runs the two baselines, the single-point step reference at 10 Ω, Marti transformation-frequency comparisons and coarse ULM residue-limit and interpolation comparisons in the three cases.
Run the 320 kV unarmoured DC cable without the fits with `julia --project=test benchmarks/cable_320kv_no_armour_dc_bipole.jl pscad`.
The 320 kV armoured DC cable uses separate sheath and armour bonds and a pulse window from its modal delays.
The 525 kV subsea armoured DC cable accepts the same arguments, with separate sheath and armour bonds and a 525 kV pulse window from its modal delays.
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

The run folder holds `segment.json`, `line_models.json`, `references.json`, `<scenario>_response_<index>.json`, `coarse_baseline.json`, `sampling_comparison.json` and `line_model_<index>_<quantity>_1_fitting_error.xlsx` workbooks.
Archives keep fit diagnostics, DC targets and achieved values, core 1 KPIs, wall times and transient RMS deviations scaled by the reference peak.
Growth ratios compare peaks in the final two 10 % windows.
A growth flag identifies a ratio that exceeds one and 110 % of its reference ratio.

Set `LINECABLEMODELS_PSCAD_CONFIG` to the station configuration path.
The section runs FD-Phase orders 10, 20 and 30 with each DC and passivity setting, and default FD-Mode.
It covers the single-point step at 10 Ω and the pulse.
It saves station identity, responses, comparisons, statuses with failure messages, and native projects, logs and outputs under `pscad_<scenario>_<index>/` in the run folder.
A station failure produces a status file, and execution continues.

`time_domain_tables.jl` reads the newest complete time-domain, PSCAD and
model-settings runs with the current archive filenames for the three DC cables.
Run `julia --project=test benchmarks/time_domain_tables.jl [output folder]` to write
eight LaTeX tables with captions and three significant digits (longtable, booktabs,
siunitx) in `.linecablemodels/benchmarks/tables/` by default.
The script prints its source folders and lists them in `time_domain_sources.txt`.
The tables cover fitting errors, DC behaviour, three transient scenarios, PSCAD,
sampling and model settings, with captions that identify omitted diagnostics.
The generator reads archives and workbooks without a study run.

`pscad_like_reference` runs the three nominal scenarios at 10 Ω with Pollaczek
earth impedance and ideal earth admittance. All other owned formulation selections,
windows and NLT controls match the Unified references. The run saves the Unified
waveforms in `references.json` and the alternative waveforms in
`pscad_like_references.json`.

Each alternative records receiving core 1 and sheath
voltage RMS differences, absolute peaks, signed peak differences and maximum
waveform differences against Unified. RMS ratios use each Unified terminal peak.

The same argument computes the single-point step and impulse at earth resistivities
of 100, 1000 and 10000 Ω·m. The `pscad_like_earth.json` archive keeps waveforms and peaks.
Both nominal formulations and every alternative resistivity point use the
`transient_sensitivity` passivity check over the scenario's full NLT frequency band.
Separate archives keep the Z/Y eigenvalue curves, band limits and first negative
sample. A negative eigenvalue prevents that transient calculation.

Add `pscad` to compute native line constants on the comparison grid:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl pscad_like_reference
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl pscad_like_reference pscad
```

The second command uses the configured station. Its `Formulation(:pscad)` keeps
the native default earth selections and the case's dielectric admittance laws.
Both owned Pollaczek + ideal and native line constants use the same 101
logarithmic samples from 0.1 Hz to 1 MHz. Each line-parameter archive records this
comparison grid's spacing, sample count and limits alongside the full frequency
coordinates. The owned and native Z/Y archives are `pscad_like_line_parameters.json`
and `pscad_native_line_parameters.json`. Owned earth coefficients remain separately
in `pscad_like_earth_parameters.json`, with the comparison grid recorded there too.

In station runs, PSCAD reported "ACLine segment manual solve: started" for its
stand-alone native line-constants solve with Pollaczek earth, but produced no
matrix file within the 7200 s run timeout at 1001 samples for the 320 kV armoured
and unarmoured cases or at 101 samples for the 320 kV armoured case, so the native
Z/Y comparison has no data.
The same solver with Pollaczek earth completed in about a minute inside every
observed FD-Phase transient run.

Generate CSV and LaTeX comparisons from completed folders:

```sh
julia --project=test benchmarks/pscad_like_tables.jl <run folder>... <output folder>
```

The generator compares the references and reads archived PSCAD transients listed
in the main tree's `local/report/tables/time_domain_sources.txt`. Each completed
PSCAD variant has receiving core 1 RMS deviation and signed absolute-peak error
against both references. RMS ratios use each reference peak. Failed
variants keep their status and empty metrics. Output folders follow the case
name and run timestamp.

The generator reports the native comparison only when
`pscad_native_line_parameters.json` exists. Its frequency tables give scaled deviations
for all, self, mutual, outer self and cross-cable outer mutual entries. The scale
is the largest absolute owned diagonal at each frequency. Individual entry tables
keep both complex matrices.

Cross-cable outer Z mutual entries isolate earth return.
Outer self entries also include conductor impedance. Native isolated
earth self coefficients are unavailable. Native conductor equations and dielectric
export representation also contribute to total Z/Y differences.

`dc_behaviour` runs the long-window DC fault study for the three DC study scripts.
Add `pscad` as a second argument to run the station matrix and the DC fault tests:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl dc_behaviour
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl dc_behaviour pscad
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl dc_behaviour pscad supplementary
julia --project=test benchmarks/cable_525kv_subsea_armoured_dc_bipole.jl dc_behaviour pscad supplementary <base folder>
```

Core 1 at the source end receives the same voltage step as the PSCAD study.
The bonds at that end each have a resistance of 10 Ω. The receiving core has a 0.01 Ω fault to ground.
The slowest L/R mode uses the first three baseline frequencies.
The duration is the larger of 20 times that mode and 50 times the longest modal travel time. The output
step is one tenth of the shortest modal travel time.

The NLT reference must
vary by less than 0.1 % over its last 10 %. An unsettled reference stops the
study and keeps its status for inspection.

Six ULM fits use orders 16/16, time rows, all three DC treatments and both
passivity settings on the baseline grid. Final values are means over the last
10 %. Signed differences use the NLT reference.
The current into the fault is `I_fault = V_receiving / R_fault`, with `R_fault = 0.01` Ω,
for NLT, ULM and PSCAD. The same definition applies to fitted DC circuit solutions.

The native receiving current points into the line for NLT and ULM, and out of the
line into the fault for PSCAD. Each status gives the raw measured mean, its orientation
and the sign that converts it to current into the fault.
It records the signed consistency residual and the maximum absolute residual over
the final window. The relative consistency error divides that maximum by the
maximum magnitude of the current derived from voltage.
These diagnostics describe the measurement agreement without a new acceptance test.

Long-window responses record the tail variation of the fault current
and both core voltages. This is the tail range divided by the magnitude of
its mean. A constant zero signal has variation zero. A signal with zero mean
and a nonzero range, or nonfinite samples, has infinite variation.
The settled flag requires current variation below 0.1 %.

ULM statuses give pole time constants for Yc, the H groups and the complete fit.
A time constant is the reciprocal of the smallest absolute real part of its poles.

A function without poles has time constant zero. A pole with zero real part
has infinite time constant. Added-pole fits also record the separate Yc and H
rates a0 in inverse seconds and their reciprocals in seconds.
Each fit archive keeps H(0), Yc(0), their physical targets and deviations.
This also applies to fits without DC correction.

The study records the baseline core self-conductance in S/m and the main
insulation's analytic conductance. A constant loss tangent specifies conductance
proportional to frequency. This conductance vanishes at DC.

The first 30 short-window FD-Phase variants compare lossless and lossy
exports, DC correction and passivity scans. They also compare five constant
shunt conductances, two correction methods and a lower start frequency for the fit.
The first five long-window variants compare disabled DC correction with both
correction methods and disabled or scan-only passivity at the physical shunt conductance.
These keep indices 1–35.

Short-window variants 36–39 use enforcement with both DC correction methods.
They use shunt conductance G_dc or 100 G_dc, and FS = 0.5 Hz.
Long-window variants 40–55 combine FS = 0.01 or 0.001 Hz, both correction methods,
and G_dc or 100 G_dc. Each combination uses scan only or enforcement.

The long window tests the guideline's 0.01 Hz start frequency for the fit
as well as 0.001 Hz. All 20 additions use enablf = 0 and orders 20/20.

The master-library definition `Line_FrePhase_Options` in `m7-master.pslx`
gives DCCOR = 1 as Functional Form and DCCOR = 2 as Add Pole/Residue.
CPASS = 0 disables passivity checks, 1 selects Scan only, and 2 selects Enforcement.
The field `shntcab` is Cable Shunt Conductance (obsolete), in mho/m,
with component default 1e-9. The recommended enablf = 1 uses loss tangent.
The constant-conductance variants explicitly select enablf = 0.

The same definition's category for passivity enforcement supplies these defaults:

- `ET_PE = 1e-10` (eigenvalue tolerance), `MER_PE = 2` (maximum percentage error)
  and `MIT_PE = 5` (maximum iterations).
- `FSP = 0.001 Hz`, `FEP = 1 MHz`, `NFP = 1000` and `FDIS = 3` (log + linear).

The component supplies these seven defaults when the exported XML and fit
arguments omit them. Each status records them in `effective_fit`
and cites the definition in `component_defaults_source`.

A full invocation exports 55 PSCAD projects and checks their XML settings.
The exporter writes six significant digits. Archives keep requested, exported
and effective values, with separate names for the correction and passivity modes.
The `pscad` argument starts station computations.

The `dc_behaviour pscad supplementary` arguments select indices 36–55
and the long and short NLT references. They skip the six ULM fits.
An optional base folder as the fourth argument also selects each index in 1–35
whose `dc_pscad_<index>/status.toml` file is absent or has an unfinished status.
The final PSCAD statuses are `completed`, `nonfinite_response`, `line_constants_error`
and `known_station_failure`. An `exported` or `started` status needs a fresh run.
The base folder must contain `dc_behaviour.toml` for the same cable.
Final status files keep their runs out of this recovery selection.

To run an explicit list of matrix indices, use:

```sh
julia --project=test benchmarks/cable_525kv_subsea_armoured_dc_bipole.jl dc_behaviour pscad runs 34 35 "<base folder>"
```

Indices must be unique integers from 1 to 55. This mode runs the selected PSCAD
variants and the NLT reference, and skips the six ULM fits. A selected index runs
even if its base status is final.

The script writes all selected results into a new run folder and records the
base folder and recovered indices in its summary. It reads the base folder
without changing its contents.
The same arguments work with all three DC cable scripts.

Each `dc_pscad_<index>/status.toml` keeps requested and exported settings,
loss tangents, line formulation, log paths, native failure messages and
available response metrics. The `dc_nlt` and `dc_ulm_<index>` folders keep
statuses, responses and fit archives. The `dc_behaviour.toml` file records
conductances, window choices, travel times and wall times.

A diverged PSCAD response records `status = "nonfinite_response"`,
`nonfinite_output_path`, `nonfinite_channel` and `first_nonfinite_time_seconds`.
The script keeps the raw outputs and continues with the next run.
Missing or unmeasured channels remain errors.

A full station matrix contains 34 short and 21 long windows.
The supplementary matrix contains four short and 16 long windows,
plus any base runs with absent or unfinished status files.
Native fits and compilation add station wall time.

The standalone generator groups run folders by cable and model.
For each cable, supply the base folder, its supplementary folders, and the Julia
rerun folder, in that order. The latest final attempt supplies each model's result.
An export-only record cannot replace a final station result. If no attempt is
final, the first supplied attempt gives the metadata and the table marks it incomplete.
The latest final NLT record supplies the common reference for that cable.

The generator writes CSV and LaTeX tables with longtable, booktabs and siunitx:

```sh
julia --project=test benchmarks/dc_behaviour_tables.jl <run folders>... <output folder>
```

The tables cover physical parameters, final values and errors, short-window
correction results, distinct exported settings, tail variations, pole time
constants and wall times. It discovers the model directories in each archive
and accepts full and supplementary folders. Final-window differences remain available
for all models.

The consolidated tables are `dc_pscad_short`, `dc_pscad_long`, `dc_ulm` and `dc_nlt`.
The short PSCAD table gives the voltage deviation.
It also gives the Tline errors and passivity diagnostics.
The long table gives voltage and current errors with settling, tail variation and divergence.
The ULM table adds fitted DC errors, pole time constants and reciprocals of added-pole rates.
CSV files keep the complete columns. LaTeX tables use shorter column selections.
The `dc_model_sources` table identifies the selected attempt and the other candidates.
The `dc_incomplete_pscad` table lists matrix entries that still need a final outcome.

The generator recomputes archived fault currents from stored receiving voltages
and the recorded fault resistance. Earlier archives use this study's fixed 0.01 Ω.
It recomputes the current and voltage errors against the selected NLT reference.

The `dc_current_consistency` table compares final-window means with the oriented
native current. Earlier archives stored minus the receiving-current channel for
all backends. The generator recovers that raw measured mean before it applies the
backend's physical orientation. The original archived current occupies a separate column for comparison.
This correction reads the archives without modifying them or repeating a station run.

The generator reads raw Tline logs from `line_constants_log_paths`.
If those paths are absent, it searches the model's archived native folders.
The case scripts only archive those logs.
The `dc_pscad_fits` tables give the last Hphase maximum and RMS errors before
Applying DC Correction. They also give the errors in Final Fit after correction.
The `dc_pscad_poles` tables give the final H pole counts and maximum ratio of
residue to pole. Yc errors and pole counts come from its own heading.

The log has no Yc error after correction, which remains `not found`.
The passivity tables give violation ranges, their smallest reported
eigenvalue and the scan or enforcement outcome. The patterns use native headings.
An absent heading or outcome gives `not found`.
The `dc_pscad_log_reports.toml` file keeps the extracted report text and log paths.
DCCOR and CPASS labels follow their master-library meanings.

Steady-state error columns need both that model and NLT to pass
the current variation criterion. The tables label unsettled models as unsettled
and leave their steady-state errors blank. Archives without these measurements
appear as `not_assessed` and have no steady-state error.
The file
`dc_sources.txt` lists the source folders.

Each ULM fit records its own DC circuit solution from `Yc(model, ComplexF64[0])`
and `H(model, ComplexF64[0])`, with the fault test's source and terminal impedances.
The study compares the DC fault current and both core voltages with the settled
NLT means. These DC errors apply independently of the transient settled flag.
The two-port assembly follows the private NLT construction in
`src/engine/transientanalysis/formulas/nlt.jl`, where H acts on current waves.

The same circuit also uses the physical DC targets retained by corrected fits.
Uncorrected fits use the independent physical DC limit. The archives record the
target solution and its signed differences from the NLT means.
The generator writes `dc_solutions` for all three quantities and `dc_target_check`
for the physical-target comparison. The `dc_fits` table gives each fit's current
errors, settled flag, slowest time constant and reciprocals of added-pole rates.

`pscad_options` compares fitting and simulation controls in the three DC study scripts.
A full invocation runs the Julia calculations. Add `pscad` to compute the native variants:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl pscad_options
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl pscad_options pscad
julia --project=test benchmarks/pscad_options_tables.jl <run folders>... <output folder>
```

To replace selected PSCAD runs in an existing folder, use:

```sh
julia --project=test benchmarks/cable_320kv_armoured_dc_bipole.jl pscad_options pscad runs 3 5 7 9 "<base folder>"
```

Indices are unique integers from 1 to 62 and select the same variants in both scenarios.
The base folder must contain `pscad_options.toml` for the same cable.
The `pscad_options pscad resume <base folder>` command completes an interrupted run of the same
cable and reuses completed PSCAD runs, Julia fits, transients and NLT references.
Append `check` after the folder to print the run plan without numerical work or station access.

The original indices 1–52 keep their settings. The tighter residue limits occupy
53–56 at order 20 and 57–60 at order 40, in ascending limit order.
Indices 61 and 62 enforce passivity with the scan band equal to the fitting band,
starting at 0.5 and 0.001 Hz, respectively. The script prints these indices for both scenarios.

This mode replaces only the selected PSCAD directories and statuses. It uses the archived
travel times and NLT references, and skips the Julia computations and reference calculations.
The generator accepts the same folder after these replacements.

The scripts print both scenarios' indices: `Numf=100` uses 2, 4, 6, 8.
`Numf=1000` uses 3, 5, 7, 9.

The core step uses single-point bonding at 10 Ω. The impulse uses the same bonds,
double exponential, duration and NLT settings as `time_domain`.
The step window is twenty times the shortest modal travel time, with 4000 intervals.
The impulse lasts at least 3 ms and three shortest travel times, sampled at 0.25 µs.

NLT uses `N=4096` for the step. The impulse uses the next power of two that reaches
2 MHz. Both use `Tfactor=1.25` and the default Hanning taper.
The time-step comparisons use NLT references at matching output times.

The FD-Phase baseline uses the lossy line formulation, `FS=0.5`, `FE=1e6`,
`Numf=100`, `YMaxP=AMaxP=20`, `YMaxE=AMaxE=0.2`, `MaxRPtol=2e6`,
`DCenab=0` and `CPASS=1`. The passivity scan uses the master-definition defaults
`FSP=0.001`, `FEP=1e6`, `NFP=1000` and `FDIS=3`.
The baseline interpolates travel times. Each comparison keeps the other controls fixed.

| Comparison | PSCAD settings | Julia settings |
|---|---|---|
| Fitting band and samples | `FS` = 0.5, 0.05, 0.01, 0.001 Hz, each with `Numf` = 100, 1000. | The lower limit = 1, 0.5, 0.1, 0.01, 0.001 Hz at 25 intervals per decade. Density = 25, 75 at the baseline band. |
| Residue/pole limit | `MaxRPtol` = 2e6, 1e9, 1e12 at 0.5 Hz. Appended limits = 1e2, 1e3, 1e4, 1e5 at orders 20 and 40 with error targets 0.2 and 0.02 percent, respectively. | `h.max_ratio` = 100, 2e6, infinity on the baseline and 1 Hz grids. |
| Order and error target | Both maximum orders = 5, 10, 15, 20, 30, 40, each with both error targets = 2, 0.2, 0.02 percent. | Both orders = 4, 8, 12, 16, 24. |
| Passivity | `CPASS=2`, `MER_PE` = 2, 0.2, 10, each with `ET_PE` = 1e-10, 1e-6. Appended runs set `FSP=FS` at 0.5 and 0.001 Hz. | Enforcement off and on. |
| Time step and interpolation | Shortest travel time divided by 2, 5, 10, 20, 50, each with `Interp1` = 0, 1. | The same steps, each with `interp=:single`, `:double`. |
| Transformation frequency | FD-Mode `FT` = 50, 500, 2000, 1e4, 1e5, 1e6 Hz. | J. Marti `transformation.frequency` at those exact samples. |

The ULM baseline uses the 227-frequency grid, 16 poles for both sections,
`h.max_ratio=Inf`, `dc=(;)` and `passivity=false`.
Maximum pole counts equal the requested orders, so adaptive order growth cannot
confound the order comparison. The fits omit time rows and the step check.
The comparisons use the ordinary time-domain interpolation and physical polarity defaults.

J. Marti uses its default Bode controls. Its grid includes the requested transformation frequency
to avoid substitution by a nearby sample.

Each case uses the package API directly. The `options_fit_<index>` folders keep
fit observations, pole and residue arrays, fitting errors and passivity reports.
The enforced ULM fit keeps its baseline report before enforcement and its own
report afterwards. Response folders store voltage and current observations.
The `pscad_options.toml` file records the scenario controls, travel times and wall times.
Failed Julia fits keep their exception type and message in their status file.

The study exports and checks 124 PSCAD variants locally, which include
two baseline runs. The check reads back each requested field and the integration
step from XML. An invocation without `pscad` records `export_checked` and opens
no station connection. Native runs keep their copied projects, logs and outputs.

`LineConstantsError` records `line_constants_error` with its native diagnostic.
`NonfiniteOutputError` records `nonfinite_response` with the channel, output path
and first nonfinite time. Other exceptions keep their type as `execution_error`.
These records preserve the cause of each failure.
At 40–75 seconds per native run, 124 runs take about 83–155 minutes per cable.

The standalone generator writes CSV and LaTeX tables under `<output folder>/<cable>/`.
It accepts several run folders per cable, with the latest final outcome taking
precedence. Export records cannot replace a final outcome.
The tables give core RMS and maximum absolute waveform deviations at the receiving end,
scaled by the NLT core peak. They also give core and individual sheath peaks and their signed errors.
Arrival uses 10 percent of the core peak.

Growth appears for both ends of each terminal, using the final two ten percent windows.
A flag requires a peak ratio above one and above 110 percent of the NLT ratio.
The model's final-window peak must also exceed NLT's by more than one percent
of the scenario scale, the NLT receiving core peak.
The growth table gives final-window amplitudes relative to that scale.
The growth-count table compares both rules on the same completed runs.

Failed and incomplete rows have empty waveform measurements.

Pole frequencies use `abs(pole)/(2π)`. For the impulse, the spectral limit contains
99 percent of the continuous source spectral energy, capped by the NLT upper frequency.
An ideal step has unbounded bandwidth. Its reported pole counts use the NLT upper
frequency as the resolved limit. The pole table keeps the definition and individual
ratios of residues to poles.

The maximum propagation ratio is dimensionless.
Admittance and impedance ratios use separate columns in S and Ω.

Native logs that report counts without pole locations leave
the count beyond the spectral limit explicitly unavailable.

The generator extracts Tline fitting errors, pole counts and maximum ratios of residues to poles.
It also records scan bands, violation ranges and enforcement outcomes with the original report text.
The passivity-range table labels each range below, inside or above the fitting band.
A range that crosses a band limit lists every region it overlaps.
A missing heading remains `not found`.

Enforcement logs give initial, final reported and worst reported eigenvalue minima.
Their final reported minimum describes an enforcement iterate, with its frequency.
If the log omits an outcome, the generator reports `iteration limit reached, violation remaining`
only when the main-iteration count equals `MIT_PE` and the last reported minimum is below `-ET_PE`.
Other runs without an outcome report `no outcome reported`.
The generator compares voltage and current at every terminal and time
with a completed `CPASS=1` run having the same fit and time grid.
A failed enforcement with an identical transient states
`enforcement failed; unenforced model simulated`.
Comparison scans identify their source run and its scan band, independently of the requested scan band.

The output root also contains `pscad_options_cross_cable.csv` and its LaTeX table.
They place the cables side by side, with each scenario's settings and a default row for every study.
The transformation comparison uses 2000 Hz as the default for FD-Mode and J. Marti.

Core RMS and signed core peak errors use the NLT core peak at the line's far end as scale.
Sheath error uses the largest absolute peak difference across sheaths at the far end.
The scale is the largest sheath peak in the NLT reference for that scenario.
The table gives the ratio in percent.
Variant identities keep repeated default settings distinct, and failed runs have blank measurements.

ULM passivity rows give reports before and after enforcement, fit-error changes
and the transient change against the unforced baseline.
The summary ranks the sampled settings by NLT deviation for each scenario.
Recommendations need the measured fit, passivity and transient evidence together.

Interpret the option comparisons through these measured effects:

| Control | Evidence for the setting |
|---|---|
| Band and samples | A lower start frequency adds slow dynamics. More samples resolve the same band more closely. Compare their separate ULM sweeps before attributing an improvement to either. |
| Residue/pole limit | Large propagation ratios show cancellation between rational terms. The bound trades fitting error against amplification of time discretization errors. Check both errors together. |
| Order and error target | Compare transient errors with the pole counts above the event's spectral limit. Extra poles can improve a frequency fit without improving the transient. |
| Passivity | Compare the violation reports before and after enforcement, the fitting-error change and the transient change. Each backend supplies its own passivity assessment. |
| Time step and interpolation | Compare arrival and waveform convergence as the step shrinks. Compare interpolation settings at the same step. |
| Transformation frequency | J. Marti's off-diagonal ratios measure the coupling discarded by its constant transformation. Compare the resulting transient over the event's spectrum. |

## Verification checks against transient outcomes

`verification_tables.jl` compares archived model checks with measured errors
against the exact NLT reference. It derives tables from saved results:

```sh
julia --project=test benchmarks/verification_tables.jl <run folders>... <output folder>
```

Supply the time-domain, model-settings and PSCAD folders listed in
`time_domain_sources.txt`, the consolidated DC behaviour folders, the option study
folders and the transient sensitivity folders. Supply each folder once.
Add the PSCAD-like reference folders, which contain `pscad_like_references.json`, for the earth-model comparison.
DC reruns follow argument order. A final result replaces an earlier result.
Different time steps and interpolation settings remain separate configurations.

The generator counts detections in each configuration class and in the pooled data.
The classes describe contrasts against each study's baseline.

The generator keeps each row's configuration class before it compares the earth models.
The comparison requires equal cable and scenario identifiers, with exactly equal window and step controls.
The archived reference time vectors must also be equal.
The generator checks native text timestamps as before and compares samples directly.
Other grids keep only the unified outcome, with an explicit reason.

A row that fails against unified NLT and passes against PSCAD-like NLT belongs to `earth model`.
A row that fails against both keeps its configuration class and records both outcomes.
Each failure definition makes its own reassignment.
`cause_class` uses the strict limits, and `cause_class_loose` uses the loose limits.
`configuration_cause_class` keeps the original class.

Short DC windows map to `step_single_point`.
Long DC windows have no reference with the same controls.

`verification_earth.csv` gives both outcomes and comparison availability for every PSCAD row.
`verification_class_changes.csv` records counts before and after the reassignment.
The detection tables and silent lists use the class for their own failure definition.

| Study | Configuration rule |
| --- | --- |
| Time domain, Julia | The baseline is the dense ULM 16/16 fit with H time rows, passivity enforcement, `dc=()` and single interpolation. Order, time rows, passivity, DC correction and coarse sampling change the fit class. Marti and FDCM change the model structure. FDCM uses a frequency-dependent transformation. |
| Time domain, PSCAD | The baseline is FD-Phase, orders 20/20, `DCenab=0`, `CPASS=1`, with the original band, window and interpolation. Order, DC and passivity changes belong to fit settings. FD-Mode belongs to model structure. |
| PSCAD options | Compare each backend with its archived variant 1. Divisor zero and native `Interp1=1` or Julia `interp=single` are the simulation defaults. Compare actual controls, so duplicate default values remain baseline. Transformation-frequency choices of FD-Mode and Marti belong to model structure. |
| Model settings | The dense ULM baseline is the same as the time-domain baseline, with `max_ratio=Inf`. Both archived ULM fits use the coarse grid and belong to fit settings. Double interpolation adds simulation settings. Marti uses the dense grid and belongs to model structure. |
| DC behaviour | The baseline uses a lossy line with loss tangent. Native settings are `FS=0.5`, orders 20/20, `DCenab=0`, `CPASS=0`. ULM uses 16/16 with time rows, `dc=()` and `passivity=false`. Long scenarios add the long-window DC class. Changed fitting controls add fit settings. Lossless dielectric or constant conductance changes add a separate line-settings class. |

The original scenario windows are the simulation defaults in each study.
Changes to several controls within one class remain in that class.
Changes across classes receive a separate mixed class, with every component retained.
Unreassigned mixed rows do not enter either component's count.

Coarse-grid Marti combines fit settings and model structure.
Coarse-grid ULM with double interpolation combines fit and simulation settings.
The summary states each rule, and `verification_mixed.csv` lists every mixed row.

The output contains CSV and LaTeX tables for configurations and scenarios,
check detection counts, and silent failures. Separate tables keep passivity ranges
and terminal measurements. The coverage table records excluded runs.
Physical-input NLT responses from transient sensitivity have no fitted model
or independent model error. Fits without transient outcomes also remain outside
the detection counts.

Reported fit-error flags use maximum entry errors above 0.2 and 2 percent.
A separate flag uses reported RMS error above 2 percent.
FD-Mode log errors keep their raw values because the native column has no unit label.
The percent-threshold checks remain unavailable for these rows.

DC-target flags use relative matrix deviations above 1 percent.
Another flag checks whether the spectral radius of H(0) exceeds one.
Propagation residue to pole ratios use thresholds of 1000 and 2000000.
These dimensionless comparisons exclude admittance and impedance ratios, which have units.

Empty fields identify unavailable checks. J. Marti passivity is unchecked.
A silent failure means no available check flags the response.
It does not mean every check passed.

The default failure definition includes divergence and these absolute errors:

- Core peak error above 2 percent.
- Scaled sheath peak error above 10 percent.
- DC endpoint error above 5 percent.

The second definition uses core and sheath limits of 5 and 20 percent.
Each sheath peak is the maximum absolute voltage at the receiving end over the scenario window.
The energized-sheath measure is the signed peak error divided by its own reference peak.
That sheath belongs to the cable with the driven core.

The scaled measure uses the worst absolute peak error across receiving sheaths.
Its denominator is the largest reference receiving-sheath peak of the scenario.
Both measures use percent. Failure criteria use the scaled measure.
A zero denominator leaves its measure unavailable.

The endpoint compares final-window voltage means for finite responses, whether settled or unsettled.
“Fine” means no threshold exceedance in retained outcomes.
Unavailable checks have separate counts for failed and fine responses.

Growth checks use the final two ten-percent windows at every voltage terminal.
The standalone check requires a peak ratio above 1.1.
The increase must also exceed one percent of the model receiving-core peak.
The divergence outcome uses the NLT comparison rule.
Its ratio must exceed one and 110 percent of the NLT ratio.
Its final-window peak must exceed NLT by more than one percent of the NLT receiving-core peak.

The NLT comparison requires the reference.
The generator reports it separately from the independent growth check.

Passivity ranges come from `PSCAD.read_pscad_fit_report` or the archived Julia fit record.
Locations distinguish frequencies below the fitting band, within it below the scenario band,
within the scenario band, and above either band.
A range that crosses a band limit keeps each location.
The scenario band uses the NLT limits.
For an impulse, its upper limit captures 99 percent of source spectral energy,
subject to the NLT upper limit.

Successful native enforcement supersedes its earlier negative samples.
Failed enforcement keeps the initial evidence for the unenforced model.
Eigenvalue normalizations differ between backends.

`verification_summary.txt` states these definitions and gives detection tables per class
and pooled, for both failure definitions.
`verification_classes.csv` gives row, failure and silent-failure counts per class.
`verification_power.csv` reports sensitivity as flagged-failed / all-failed
and specificity as quiet-fine / all-fine, with counts and fractions.
Unavailable checks stay in both denominators and contribute to neither numerator.
The fractions describe observed detection and clearance with the retained evidence.
Each denominator counts every row in its class, regardless of check availability.
A zero denominator has no rate.
`verification_silent.csv` lists each silent failure with its configuration and scenario.
It also lists errors and unavailable checks.
The combined check set reports results at each fit-error threshold and with RMS alone.
This comparison shows the effect of that threshold on silence.
`verification_sources.txt` records every input folder.
