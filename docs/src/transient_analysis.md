# Transient analysis

`TransientAnalysisProblem` holds an uncomputed `LineParametersProblem` or a fitted
`LineModel`, assignments at both ends and a time window. `TransientFormulation(:nlt)`
selects the numerical Laplace transform of the exact line two-port, and `:td`
selects time-domain simulation of a fitted model. `compute` returns a
`TransientResponse` with voltage in V, current in A and time in s.

## Terminals and waveforms

Each end is a tuple of assignments or `Matched()`. Each assignment specifies a
terminal label such as `"cable:2:core"`. Assignments can select cores, sheaths or armour.
Use `terminal_labels(line.system)` to list valid labels for terminal assignments.
Each terminal can have one assignment per end. Unlisted terminals are open.
`TerminalImpedance` connects a terminal to ground through a series resistance R
in Ω and inductance L in H. Zero values for both ground the terminal solidly.

```julia
sending = (Step("cable:1:core", 1.0; t0=0.0),
    TerminalImpedance("cable:2:core"; R=1.0),
    TerminalImpedance("cable:3:core"; R=1.0))
receiving = (TerminalImpedance("cable:2:core"; R=10.0, L=0.0),
    TerminalImpedance("cable:3:core"; R=1.0))
```

Waveforms are ideal voltage sources between their terminal and ground:

```julia
Step("cable:1:core", 1.0; t0=2e-6)                 # amplitude [V], delay [s]
Cosine("cable:1:core", 1.0, 50.0; phase=0.0)       # amplitude [V], frequency [Hz], phase [rad]
DoubleExponential("cable:1:core", 1.0, 1e4, 1e6)    # amplitude [V], decay rates [1/s]
```

The cosine is ``A\cos(2\pi f t+\varphi)``. The double exponential is
``A(e^{-\alpha t}-e^{-\beta t})``, with ``\beta>\alpha>0``.
With `:nlt`, `Matched()` terminates a whole end with the characteristic admittance
matrix ``Y_c(s)`` from the line. It replaces the tuple at that end.

## Problem and computation

Given a `line_problem` with the required geometry and physical length, declare
the time window and select the line formulas through the transient options:

```julia
problem = TransientAnalysisProblem(line_problem; sending, receiving,
    window=(tmax=2e-4, step=1e-6))
formulation = TransientFormulation(:nlt; options=(
    N=4096, Tfactor=1.25, taper=:hanning,
    line=Formulation(earth_impedance=:unified, earth_admittance=:unified)))
response = compute(problem, formulation)

matched_problem = TransientAnalysisProblem(line_problem;
    sending, receiving=Matched(), window=problem.window)
```

`tmax` and `step` use s. Both must be positive and `step < tmax`. Samples start
at zero, with `round(Int, tmax / step) + 1` points separated by `step`.
Response arrays have axes end × terminal × time (2 × n × samples).
The sending end comes first, followed by the receiving end. Currents enter the line at both ends.
`response.coordinates` retains terminal labels in physical order.

The `:nlt` options control the scan and inversion:

| Option | Default | Meaning |
|---|---|---|
| `N` | `4096` | Number of Laplace samples, an `Int` of at least two. |
| `Tfactor` | `1.25` | Sets the transform interval to `Tfactor * tmax`. Must be finite and greater than one. |
| `taper` | `:hanning` | Spectral taper before inversion: `:hanning`, `:blackman` or `:lanczos`. |
| `line` | `Formulation()` | `LineParametersFormulation` used to compute the Laplace-grid scan. |

`compute` first builds the Laplace grid, with ``T=\mathrm{Tfactor}\,t_{\max}``:

```math
f_k = \frac{2k-1}{2T},\qquad
\sigma = \frac{\ln(N^2)}{T},\qquad
s_k = \sigma + j2\pi f_k,\qquad k=1,\ldots,N.
```

It computes the line parameters at those frequencies and that real Laplace
abscissa. The original frequency samples in `line_problem` are replaced for
this scan. At each grid point, it constructs the exact line two-port, applies
the terminal assignments and solves for the terminal voltages and currents.
Terminals removed by the selected line formulation cannot receive assignments.
The tapered spectra are inverted at each requested time:

```math
x(t) \approx \frac{2e^{\sigma t}}{T}
\operatorname{Re}\!\left[
\sum_{k=1}^{N} w_k X(s_k)e^{j2\pi f_k t}
\right],
```

where ``w_k`` is the selected taper. `details(response).data.nlt` retains the
effective grid and taper controls and an aliasing estimate. These numerical
records support inspection. Physical accuracy needs separate assessment.

## Observations and archives

`observe` selects the stored voltage or current. End and terminal selectors
accept indices or labels:

```julia
voltages = observe(response, voltage)  # end × terminal × time [V]
receiving_voltage = observe(response, voltage, "receiving", "cable:2:core", :)
sending_current = observe(response, current, 1, 1, :)

observed = ObservedResult(response, (voltage, current))
display(report(observed))
save(observed, "transient.json")
save(observed, "transient.jls")
```

`ObservedResult` retains detached values, terminal identities and time
coordinates. Reports use time rows with a column for each end and terminal.
Saved observations retain these quantities and coordinates for later reporting.

## Response KPIs and comparison

Each KPI accepts a response and `voltage` or `current`. Its arrays are indexed
by end × terminal (2 × n). The final value is the mean over the last 10 % of
time samples, rounded up to at least one sample.

```julia
extrema = peak(response, voltage)
extrema.value  # largest absolute sample [V]
extrema.time   # first time of that peak [s]
arrival = arrival_time(response, voltage; fraction=0.1)  # [s]
decrement = damping(response, voltage)                   # dimensionless
offset = dc_offset(response, voltage; expected=0)        # [V]
```

`arrival_time` returns the first sampled time when the absolute signal equals
or exceeds the requested fraction of its peak. The fraction must be in (0, 1]. A zero
peak gives `NaN`. Referencing the peak allows a pulse that settles at zero to
have an arrival time.

`damping` returns ``\ln(a_1/a_2)``, where ``a_1`` and ``a_2`` are the first
two strict local maxima of the signal minus its final value, strictly after
arrival at fraction 0.1. Fewer than two maxima gives `NaN`. For
``Ae^{-at}\cos(\omega t)``, the decrement approaches ``a\,2\pi/\omega`` as
the final-value estimate approaches zero and time resolution improves.
`dc_offset` subtracts `expected` from the final value. `expected` defaults to
zero and uses the same units as the selected signal.

`compare(reference, result, selector)` returns `RMSError`. For each end and
terminal it computes RMS error over time, normalized by the reference RMS:

```math
E_{\mathrm{abs}} = \sqrt{\frac{1}{N_t}\sum_k |x_k-r_k|^2},
\qquad E_{\mathrm{rel}} =
\frac{E_{\mathrm{abs}}}{\sqrt{\frac{1}{N_t}\sum_k |r_k|^2}}.
```

```julia
import LineCableModels.Engine: compare

error = compare(reference, response, voltage)
error.absolute  # end × terminal [V]
error.relative  # dimensionless fractions
details(error).data
```

The operands must have the same terminal coordinates and time axes equal
within relative tolerance 1e-12, with zero absolute tolerance. Relative error
is `missing` where the reference RMS is zero. Details retain the quantity,
`:reference_rms` normalization, sample count and time window.

## Time-domain simulation of a fitted line model

Load UniversalLineModel, fit a line model, then select `:td`:

```julia
using UniversalLineModel

phase = compute(line_problem, Formulation())
modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation())
segment = PropagationParameters(modal; line_length=line_length(line_problem))
model = compute(LineModelProblem(segment), LineModelFormulation(:ulm))
problem = TransientAnalysisProblem(model; sending, receiving,
    window=(tmax=2e-4, step=1e-7))
response = compute(problem, TransientFormulation(:td;
    options=(interp=:single, enforce_sign=false, initial_solve=false)))
```

The time step must be at most the smallest fitted propagation delay.
The `:td` options pass to the UniversalLineModel network and time loop:

| Option | Default | Meaning |
|---|---|---|
| `interp` | `:single` | Delayed-history interpolation. `:single` uses one linear interpolation per step. `:double` splits the delayed step into two segments at the grid sample inside it when the delay is not a whole number of steps. |
| `enforce_sign` | `false` | When true, sets the companion conductance diagonal to absolute values and its off-diagonal entries to negative absolute values. Passed to `add_line!`. |
| `initial_solve` | `false` | When true, solves the companion network at time zero with zero history sources. Otherwise the initial values are zero. Passed to `simulate`. |

`Matched()` terminates a whole end with the fitted characteristic admittance
``Y_c(s)`` of the line model in `:td`, and with the exact ``Y_c(s)`` in `:nlt`.
Comparisons between them include the admittance fitting error at the
matched end, along with the line fitting and numerical inversion errors.
A `ModalModel` also throws with `:td`.

## Choosing a formulation

`:nlt` computes the exact two-port of an uncomputed `LineParametersProblem` on
its own Laplace grid, then inverts the spectra numerically. `:td` simulates a
fitted `LineModel` with a `WidebandModel` payload through UniversalLineModel.
`Matched()` is available with both formulations.

Use the same terminal assignments, line length and time window to compare the
two responses. After the fit and `:td` computation above, compute the `:nlt`
reference on that window:

```julia
import LineCableModels.Engine: compare

reference_problem = TransientAnalysisProblem(line_problem; sending, receiving,
    window=problem.window)
reference = compute(reference_problem, TransientFormulation(:nlt;
    options=(line=Formulation(),)))
voltage_error = compare(reference, response, voltage)
current_error = compare(reference, response, current)
```

The reference uses the same line formulation as the scan used for the fit.
The returned errors compare the stored samples over this window.

## Sweeps

Numeric `Grid` values in terminal assignments lift the transient problem into
a `Gridspace`. The existing traversal computes one response per point:

```julia
receiving = (TerminalImpedance("cable:2:core"; R=Grid((1.0, 10.0))),
    TerminalImpedance("cable:3:core"; R=1.0))
problems = TransientAnalysisProblem(line_problem; sending, receiving,
    window=problem.window)
responses = compute(problems, formulation)
peaks = [peak(result, voltage) for result in responses]
```
