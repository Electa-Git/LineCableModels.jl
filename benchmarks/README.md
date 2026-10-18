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
