# Line models

A fitted line model represents characteristic admittance by a rational matrix and
propagation by groups of rational matrices with delays. `LineModelProblem` holds a
completed propagation segment. `LineModelFormulation` selects `:ulm`, `:fdcm` or
`:marti`. The default selection is `:ulm`. Load `UniversalLineModel` to fit `:ulm`
and `:fdcm`. The `:marti` selection currently throws an error.

## Problem and computation

Start with a `line_problem` sampled at 31 logarithmically spaced frequencies from
1 Hz to 1 MHz, as in the extension example. Compute its phase and modal parameters,
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

The returned `LineModel` retains the fitted coefficients and diagnostic records.
Selectors expose poles, residues, propagation delays in seconds, relative errors
and passivity diagnostics:

```julia
poles(model)
residues(model)
delays(model)
fitting_error(model)
passivity(model)
```

## Line-model sweeps

A selection `Grid` forms a formulation `Gridspace`. This example compares `:ulm`
and `:fdcm` with eight poles for each rational approximation and passivity
enforcement disabled:

```julia
formulations = LineModelFormulation(Grid([:ulm, :fdcm]);
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

Grids nested inside formulation option records currently remain inside those
records. Where such grids should lift remains an open design question.
