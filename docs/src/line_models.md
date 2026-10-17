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
| `dc` | `(;)` | A non-empty value throws before either fit. No direct-current targets pass to UniversalLineModel. |

The extension calls both fits with `geom=nothing` and a `LineFrequencyData`
adapter of the segment. It measures passivity after either fit, including when
`passivity=false` disables enforcement.

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
| `dc` | `targets` and `achieved`, both `nothing` for the current fits. |
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
