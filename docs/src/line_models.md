# Line models

A fitted wideband line model represents characteristic admittance by a rational
matrix and propagation by groups of rational matrices with delays.
`LineModelProblem` holds a completed propagation segment.
`LineModelFormulation` selects the fitting method, and `compute` returns a
`LineModel` with coefficients and diagnostic records.

## Problem and computation

Start with a `line_problem` sampled at 31 logarithmically spaced frequencies from
1 Hz to 1 MHz. Compute its phase and modal parameters,
then bind a segment of length 1000 m:

```julia
using LineCableModels, UniversalLineModel

phase = compute(line_problem, LineParametersFormulation())
modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
segment = PropagationParameters(modal; line_length=1000.0)
problem = LineModelProblem(segment)
formulation = LineModelFormulation(:ulm;
    options=(yc=(npoles=8,), h=(npoles=8,)))
model = compute(problem, formulation)
```

`LineModelProblem(segment)` checks that the frequency vector is non-empty,
finite, positive and strictly increasing. Fitting uses the real frequency axis
only: the Laplace abscissa must be ``\sigma=0``. The segment length must be
positive. Series impedance, shunt admittance, propagation constants and both
modal transformation matrices must have finite entries.

## Formulations and options

The default selection, `:default`, resolves to `:ulm`. Load UniversalLineModel
for the `:ulm` and `:fdcm` fits. Both return a `WidebandModel` payload.

| Selection | Method |
|---|---|
| `:ulm` | Universal line model (1999), through `UniversalLineModel.fit_ulm`. |
| `:fdcm` | Frequency-dependent cable model (2016), through `UniversalLineModel.fit_fdcm`. |
| `:marti` | Constant-transformation modal line model (1982), returning a `ModalModel`. |

For `:ulm` and `:fdcm`, each section of the options defaults to an empty named tuple. Omitted
controls within `yc` and `h` use the defaults of the UniversalLineModel method
that receives them.

| Option | Default | Passed to UniversalLineModel |
|---|---|---|
| `yc` | `(;)` | `yc` in `fit_ulm` or `fit_fdcm`, then keyword arguments of `fit_yc` for characteristic admittance. `(npoles=8,)` sets eight poles. |
| `h` | `(;)` | `h` in `fit_ulm`, then keyword arguments of `fit_h` for propagation. For `:fdcm`, its entries pass as keyword arguments directly to `fit_fdcm`. |
| `passivity` | `(;)` | Passed to `fit_ulm` for `:ulm`. For `:fdcm`, a separate call to `enforce_passivity` follows the fit unless this is `false`. With `true` or `(;)`, enforcement uses its default controls. A named tuple supplies enforcement keywords. |
| `dc` | `(;)` | For `:ulm`, select `(method=:add_pole,)`, `(method=:add_pole, k=0.1)` or `(method=:functional_form,)`. Correction needs the physical `line`. The `:fdcm` fit rejects non-empty values. |
| `line` (`:ulm` only) | `Formulation()` | Line-parameter formulation for exact Laplace samples when `LineModelProblem(segment; line=...)` supplies the physical problem. |

The extension fits a `LineFrequencyData` adapter of the segment.
A supplied physical line provides exact Laplace samples for `:ulm`.
The extension measures passivity after either fit, including when
`passivity=false` disables enforcement.

## DC correction

The `:ulm` fit offers the Add pole and Functional form methods of
[De Silva, Gole and Wedepohl, IPST 2007](https://www.ipstconf.org/papers/Proc_IPST2007/07IPST024.pdf).
[PSCAD offers both methods](https://www.pscad.com/webhelp-v502-ol/Master_Library_Models/Transmission_Lines_Cables/Distributed_Line_Models/fd_phase_options.htm)
for its phase model. An empty `dc=(;)` keeps the ordinary fit.
The `:fdcm` and `:marti` fits refuse DC correction.

```julia
problem = LineModelProblem(segment; line=line_problem)
formulation = LineModelFormulation(:ulm; options=(dc=(method=:add_pole,),
    line=LineParametersFormulation(insulation_admittance=:lossy), passivity=false))
model = compute(problem, formulation)
model.dc.targets
model.dc.achieved
details(model).data.dc
```

The physical line and its `line` formulation follow the same consistency checks
as the exact step response. Finite insulation resistivity supplies positive DC
shunt conductance with `insulation_admittance=:lossy`.

LineCableModels computes the targets from the line's ``s\to0`` limit.
One scan evaluates ``f=(10^{-6},10^{-9},10^{-12})`` Hz at ``\sigma=0``.
For each pair of frequencies, the maximum entry difference in ``Z`` and ``Y``
is divided by the largest entry magnitude at the smaller frequency.
All three relative differences must be at most ``10^{-6}`` for each matrix.
The smallest-frequency matrices define

```math
D_Y=Z^{-1}\sqrt{ZY},\qquad D_H=\exp(-\sqrt{YZ}\,\ell).
```

LineCableModels symmetrizes ``D_Y`` by its transpose.
For each target, the largest imaginary magnitude divided by the largest real
magnitude must be below ``10^{-6}``. The real parts supply the DC targets.
``D_H`` follows the adapter's current-wave orientation and uses the segment length ``\ell``.
The computation throws when the parameters fail to converge or either target exceeds the imaginary tolerance.

Add pole appends a low-frequency pole and chooses its residue to reach the target.
With `k`, its decay rate is ``a_0=2\pi f_{\min}k``, where ``0<k<1``.
Without `k`, one scan supplies exact admittance and propagation at 31 logarithmic
frequencies from ``f_{\min}/1000`` through ``f_{\min}``.
These samples pass as `dc.reference=(s=2π*im*f, F=matrices)` to UniversalLineModel.
ULM selects the pole from the reference and fitted frequency bands.
Functional form identifies the rational model with a fixed DC term through
`dc=(method=:functional_form,)`.

`model.dc` retains `targets` and `achieved`, each with `admittance` and `propagation`.
The achieved matrices evaluate the final fitted model at zero.
`details(model).data.dc` retains convergence measures and ULM's admittance and propagation diagnostics.
For Add pole, these include `method`, `k`, `a0` in rad/s and `reference_range` in Hz.
Passivity enforcement can change the achieved propagation value after correction.

With zero shunt conductance, the targets ``Y_c(0)=0`` and ``H(0)=I`` leave the DC resistance and steady-state current undetermined.

## Time-domain rows and the step check

For `:ulm`, supply the physical line problem and select its line-parameter
formulation to evaluate the exact propagation step response on the Laplace plane:

```julia
problem = LineModelProblem(segment; line=line_problem)
formulation = LineModelFormulation(:ulm; options=(h=(td=true,),
    line=Formulation(earth_impedance=:unified, earth_admittance=:unified)))
model = compute(problem, formulation)
```

The line problem must have zero Laplace abscissa and the segment's frequencies.
The extension compares the supplied physical inputs and resolved formula
selections with the segment's recorded computation before fitting.
On a mismatch, the error specifies the first differing field.
Without these records, validation warns that it cannot check consistency.

The `h` controls pass unchanged to UniversalLineModel:

| Control | Default | Meaning |
|---|---|---|
| `td` | `false` | Add exact step-response samples to the residue fit alongside frequency samples. Poles and delays keep their frequency-domain values. Needs `LineModelProblem(segment; line=...)`. |
| `td_N` | `4096` | Laplace terms for the arrival window. Late-time scans use `max(1024, td_N ÷ 4)`. |
| `step_check` | `true` | Compare fitted and exact step responses whenever the physical line is supplied, including when `td=false`. |
| `step_tol` | `0.1` | Relative error threshold for ULM's step-check warning. |

ULM samples an arrival window around the fitted delays and, by default, late
times from twenty times the largest delay to `10/(2π * 1e-3)` seconds. Time rows
have relative weight `td_weight=0.03`. The default `td_window=true` includes
arrival rows, and `td_late=true` includes late rows. Set `step_check=false` to
disable the check. Exact samples still supply the time rows when `td=true`.

The step check and `td` rows use Laplace grids with highest frequency
`(2N−1)/(2T)`. Arrival grids extend far above the fitted band and can exceed the
line model's validity range, triggering the quasi-TEM warning above 100 MHz.
In that case the arrival-window error compares the fit with line parameters
evaluated outside their validity range. ULM's `td_window=false` keeps only
late-time rows when `td=true` and `td_late=true`.

For a fit with time rows or the step check, the extension computes one Laplace
scan for the arrival window. ULM requests another scan for each two decades of
late times. Each scan computes all its line-parameter samples in one `compute`
call. The adapter keeps samples during that fit. ULM performs its own eigenvalue
calculation, windowing and inversion. Every fit repeats these scans across
formulation combinations.

The result's `details(model).data.step_check` records whether the check ran.
When it ran, `step_err` and `step_err_late` record the largest element errors in
the arrival and late-time windows, relative to the largest exact sampled value.
Both errors are `NaN` when the check did not run, for both wideband fit methods.
Without a physical line, the fit uses its frequency samples and `td=true` throws.
The `:fdcm` and `:marti` options reject `line`.
If a problem supplies `line` for either fit, validation warns and fitting uses
the segment.

## J. Marti modal fit

The `:marti` fit uses a constant real current transformation ``T``. It takes the
real part of the segment's modal current transformation at the selected sample
and normalizes each column to unit norm. The voltage transformation is
``T_v = T^{-\mathsf{T}}``. At each frequency, the fit uses the diagonals of
``Z_m = T_v^{-1} Z T`` and ``Y_m = T^{-1} Y T_v`` to calculate the modal
characteristic impedance ``Z_{c,i} = \sqrt{z_{m,i}/y_{m,i}}`` and propagation
``A_i = \exp(-\gamma_i\ell)``, with
``\gamma_i = \sqrt{z_{m,i}y_{m,i}}`` chosen to have nonnegative real part.

Each propagation section retains a travel delay in seconds. Its initial value
comes from the highest-frequency propagation constant. A search over 300 delay
values within 10% of that estimate minimizes the RMS phase difference between
the delay-compensated propagation and the minimum phase implied by its magnitude.
The Bode gain and phase relation uses a discrete Hilbert transform on the
log-frequency grid with trapezoidal integration.

Bode asymptotes approximate each modal impedance and propagation magnitude with
slopes in multiples of 20 dB per decade. Corners supply starting real poles.
Variable projection then refines their logarithms by Levenberg–Marquardt least
squares. Every pole remains real and negative. For each trial, a linear least-squares
solve identifies a real constant and real residues against the complex impedance
or delay-compensated propagation, with weights equal to the inverse data magnitude.
Only steps that decrease the residual are accepted. Refinement stops below a
relative improvement of ``10^{-8}`` or after 50 iterations.

| Option | Default | Meaning |
|---|---|---|
| `transformation.frequency` | `nothing` | Frequency in Hz for the constant transformation. The nearest sample is used. `nothing` selects the last sample. |
| `bode.tolerance_db` | `0.5` | Magnitude tolerance in dB for the starting Bode asymptotes. |
| `bode.threshold` | `1e-2` | Propagation magnitude cutoff, strictly between zero and one. |

For each mode, the propagation fit and the delay estimate use the samples before the first one below this threshold.
The characteristic impedance keeps the full band.
`diagnostic.propagation_f_max` records the highest fitted propagation frequency of each mode in Hz.

These controls are stored in the `transformation` and `bode` named tuples.
Unknown option sections throw. The `:marti` sections also reject unknown keys.

```julia
model = compute(problem, LineModelFormulation(:marti;
    options=(transformation=(frequency=nothing,), bode=(tolerance_db=0.5,))))
diagnostic = details(model).data.marti
```

`diagnostic.offdiagonal_ratio` stores the largest off-diagonal magnitude divided
by the largest diagonal magnitude, over all frequencies, for ``Z_m`` and ``Y_m``
in that order. `diagnostic.phase_error_degrees` stores the largest absolute
propagation phase error per mode in degrees. The fitting-error matrices contain
per-mode magnitude errors on their diagonals and zeros elsewhere. Their RMS and
maximum values normalize by the largest data magnitude over frequency.
Passivity has `checked=false` and `passive=nothing`. Both DC fields are `nothing`.

## Result record

The returned `LineModel` has these fields:

| Field | Contents |
|---|---|
| `model` | The fitted payload. `WidebandModel` stores rational admittance and delayed propagation groups. `ModalModel` stores a real constant transformation and real-pole impedance and propagation sections. |
| `coordinates` | Phase labels in physical order. |
| `line_length` | Segment length in m. |
| `frequencies` | The fitted frequency samples in Hz. |
| `fitting_error` | `admittance` and `propagation` records, each with phase × phase matrices `rms`, `max` and `f_max`. Errors are relative fractions, and `f_max` gives the frequency in Hz of each maximum. |
| `passivity` | `checked`, `passive`, `bands`, `worst`, `f_worst` and `perturbation`. Bands are frequency pairs in Hz. The perturbation record is `nothing` when enforcement was disabled. |
| `dc` | `targets` and `achieved`, each with admittance and propagation matrices when correction is active, otherwise both `nothing`. |
| `details` | `ComputationDetails` with the formula identifier, options, source grid point and inputs. Timing is included when requested. |

`WidebandModel.admittance` contains `poles`, `residues` and `direct`. Each entry
of `WidebandModel.propagation` contains `delay`, `poles` and `residues`.
Residue arrays have axes phase × phase × pole.
`ModalModel` has fields `transformation`, `transformation_frequency`, `impedance`
and `propagation`. Its sections contain `constant`, `poles`, `residues` and `delay`.

For each matrix entry, the measured fitting errors normalize the absolute
difference by the largest absolute source sample of that entry over frequency.
`rms` is the RMS of these fractions, and `max` is their maximum.
Where a source entry is zero, an exact fit has zero error and any other fit has `Inf`.

## Selectors, observations and reports

Selectors expose poles in 1/s, residues, propagation delays in seconds,
relative errors and passivity diagnostics:

```julia
poles(model)
residues(model)
delays(model)
fitting_error(model)
passivity(model)
```

For a wideband payload, `poles` and `residues` return the admittance group
first, followed by the propagation groups. `delays` returns one value per
propagation group. `fitting_error(model)` returns both error records, while
`fitting_error(model, :admittance)` and `fitting_error(model, :propagation)`
select the corresponding RMS matrices.
`passivity(model)` returns the retained passivity record.

```julia
admittance_error = observe(model, fitting_error, :admittance)
propagation_error = observe(model, fitting_error, :propagation)
pole_groups = observe(model, poles)
propagation_delays = observe(model, delays)

observed = observables(model, (fitting_error,))
display(report(observed))
```

The `fitting_error` observation expands to the admittance and propagation RMS
matrices, with detached values and phase labels. `observables(model)` includes
poles, delays and both error matrices. Reports accept these observations.

## Line-model sweeps

A selection `Grid` forms a formulation `Gridspace`. This example compares `:ulm`
and `:fdcm` with eight poles for each rational approximation and passivity
enforcement disabled:

```julia
formulations = LineModelFormulation(Grid((:ulm, :fdcm));
    options=(yc=(npoles=8,), h=(npoles=8,), passivity=false))

# segments is a completed ParametricResult of PropagationParameters.
problems = Gridspace{LineModelProblem}(segments)
run = compute(problems, formulations)
model = run[1, 2]
selected_formula = formula_id(run.axes.formulations[2].formula)
selected_options = formulation_options(run.axes.formulations[2].formula).data
```

`run` is a `ParametricResult` in problems × formulations order. The problem index
varies fastest, so point `(i, j)` has linear index `i + (j - 1) * length(segments)`.
The formulation axis identifies the selection and numerical options at each point.

Grids inside option records form the same formulation space. This example
varies the admittance pole count, propagation pole count and passivity enforcement
over eight fits:

```julia
formulations = LineModelFormulation(:ulm;
    options=(yc=(npoles=Grid((6, 10)),), h=(npoles=Grid((8, 12)),),
        passivity=Grid((false, true))))
run = compute(problems, formulations)
```

The number of poles for admittance varies fastest. The number for propagation
varies next. The passivity setting varies last. Each entry on the formulation
axis keeps its selected scalar options.
