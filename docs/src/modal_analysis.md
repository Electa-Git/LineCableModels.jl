# Modal analysis and finite segments

`ModalAnalysisProblem` consumes one completed phase-domain `LineParameters`.
`ModalAnalysisFormulation` selects the decomposition. The returned
`LineParameters` contains modal coefficients, modal-to-phase voltage and
current bases `Tv` and `Ti`, and retained propagation roots `gamma`.

```julia
phase = compute(line_problem, line_formulation)
modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
restored = transform(PhaseDomain, modal)

Tv(modal)       # phase × mode × frequency
Ti(modal)       # phase × mode × frequency
gamma(modal)    # mode × frequency
Zc(modal)       # mode × frequency
Yc(modal)       # mode × frequency
Zc(modal, PhaseDomain)  # ordered phase × phase × frequency
```

The bases map modal coordinates to phase coordinates. At each frequency the
coefficients are `Tv \ (Zphase * Ti)` and `Ti \ (Yphase * Tv)`; inverse conversion
uses the retained bases. The characteristic quantities use the diagonal modal
approximation and the retained root branch. Phase characteristic matrices
preserve row and column order and need not be symmetric.

The convenience call performs the same two public actions in sequence:

```julia
modal = compute(line_problem, line_formulation; modal=:default)
modal = compute(line_problem, line_formulation;
    modal=formula(:default; options=(iteration=(convergence=1e-8,),)),
    modal_options=(offdiagonal_tolerance=1e-6,))
```

`modal_options` controls the modal action only. A nonempty value requires
`modal`. The upstream solver, including its selected backend and batched path,
completes before modal analysis begins. Upstream `on_result`, progress, and
timing belong to that upstream call; modal `on_result`, progress, and timing
belong to the downstream call. Enabling timing does not change the source
gridpoint identity. A missed iteration or coupling target is recorded and
warned about; it does not discard finite results. The iteration convergence
setting is a numerical target, not a physical error bound. Inspect
`details(modal).data.modal.diagnostics` for fallback and missed frequency
indices, Z/Y residual coupling, returned eigen-residuals, iteration counts,
and convergence flags. These are numerical diagnostics, not a scientific
acceptance decision. The scalar characteristic and response construction
uses the selected modal approximation when coupling remains.

The built-in formula selects paired voltage and current eigenvectors and
retains a forward root with nonnegative real part, or nonnegative imaginary
part for a purely imaginary root. The modal workspace supplies impedance and
admittance slices, coordinate conversion scratch, the admittance–impedance
product, eigenpair history and voltage-vector scratch. The built-in formula
tracks eigenvalues of its normalized, shifted eigenproblem; multiplying
`eigenvalue + 1` by the scale gives an eigenvalue of the physical product,
distinct from a propagation root. The selected formula allocates its
normalized eigenproblem, Hungarian assignment work and Levenberg–Marquardt
least-squares work through `initialize_buffers`. Its arithmetic remains in
`decompose!(::Val{:chrysochos2014}, ...)`. `Tv`
and `Ti` map modal coordinates to phase coordinates. Total source coefficients
are normalized by their declared source length before the finite segment is
bound.

## Finite segments and source length

`PropagationParameters` binds a modal scan to one physical segment. A source
with total coefficients retains its positive source normalization length;
the finite segment may have another length. For an external scan without a
known source length, supply the target length explicitly.

```julia
segment = PropagationParameters(modal; line_length=150.0)
next_segment = PropagationParameters(segment; line_length=300.0)
length_space = PropagationParameters(modal;
    line_length=Grid((150.0, 300.0))) # lazy finite segment Gridspace

H(segment)                              # mode × frequency
alpha(segment)                          # attenuation constant, Np/m
beta(segment)                           # phase constant, rad/m
velocity(segment)                       # phase velocity, m/s
H(segment, PhaseDomain; field=:voltage) # ordered phase voltage response
H(segment, PhaseDomain; field=:current) # ordered phase current response
```

The modal forward factors are `exp.(-gamma(segment) .* line_length(segment))`.
`alpha` and `beta` are the real and imaginary parts of the per-unit-length
`gamma`; they do not depend on the selected segment length. `velocity` is
`2πf / beta` at each retained mode and frequency. Undefined values remain
unavailable in observations rather than being clipped or assigned a finite
surrogate. The default display units for the first two are Np/km and rad/km.
The phase voltage and current responses use `Tv` and `Ti`, respectively.
Zero length gives identity propagation. Rebinding does not change the embedded
source, and a changed length receives a new gridpoint identity. Indexing a
modal or finite result by frequency slices coefficients, bases, roots, and
frequency-aligned diagnostics together without a new solve.

## Collections and observations

Completed phase result spaces transport through
`Gridspace{ModalAnalysisProblem}(phase_results)`. A single modal formulation
keeps one downstream result per source point; a formulation Gridspace forms
the ordinary product or zip selection. Completed modal results transport
through `Gridspace{PropagationParameters}(modal_results)` for finite lengths.
Scalar quantities compose with broadcasting, for example `gamma.(modal_results)`
and `H.(segment_results)`; each tensor remains one value.

`ObservedResult` retains detached quantities and the existing four sections:
gridpoint, quantities, errors, and timings. Requests can mix phase coefficients
with `gamma`, `Zc`, `Yc`, `Tv`, and `Ti`. A finite segment also offers `H`,
`alpha`, `beta`, and `velocity`. Bare complex requests acquire complete real and
imaginary pairs; a raw plotting convenience also completes a single component
while displaying only that component. An observed-input plot selects retained
components without calculating an absent representation.
Bind phase representations explicitly:

```julia
phase_Zc = Base.Fix2(Zc, (domain=PhaseDomain,))
voltage_H = Base.Fix2(H, (domain=PhaseDomain, field=:voltage))
current_H = Base.Fix2(H, (domain=PhaseDomain, field=:current))
observed = ObservedResult(segment,
    (gamma, Zc, H, voltage_H, current_H, Tv, Ti, velocity))
another = PropagationParameters(segment; line_length=300.0)
observed_pair = ObservedResult.((segment, another), Ref((H, gamma)))
plot(observed; ydata=(gamma, velocity)) # all modes overlaid for one gridpoint
plot(observed_pair; ydata=(gamma,), overlay=:coordinates)
```

Mode vectors have a mode axis and a frequency axis. `Tv` and `Ti` have phase
rows, mode columns, and frequency samples; retained coordinates keep both
label sets. Phase matrices have ordered phase rows and columns. `abs`, `angle`,
`real`, and `imag` can be requested through the existing request tuples.
For example, `@observe alpha[:, :]` and `@observe Tv[:, :, :]` keep original
indices; `@observe (H, abs)[:, :]` requests the explicit magnitude component.
Complete magnitude/angle pairs are acquired separately from the default
Cartesian components. Each component has its own physical name, unit, table,
and plot family. A vector component table has frequency rows and one column
per selected mode, named `Mode 1`, `Mode 2`, and so on. Full transformation
tables include every phase × mode entry, including off-diagonals.

To view how each conductor participates in each mode, overlay the rows of the
transformation matrices:

```julia
transform_matrices = observables(segment,
    ((Tv, abs), (Tv, angle), (Ti, abs), (Ti, angle)))
display(report(transform_matrices))
plot(transform_matrices; overlay=:rows, layout=(3,3))
```

Each component has mode panels containing conductor curves. Panel titles name
modes; legends name conductors. Selected modes fill each page left to right and
then top to bottom. For example, 18 modes produce two pages per component with
`layout=(3,3)`. Every selected coefficient remains plotted.

For a collection, acquire the same requests at each point and pass the retained
results to the same plotting call:

```julia
requests = ((Tv, abs), (Tv, angle), (Ti, abs), (Ti, angle))
transform_space = observables.(segment_results, Ref(requests))
plot(transform_space; overlay=:rows, layout=(3,3))
```

Each gridpoint starts its own figures with the same block layout. Curves and
panels from different points never share a figure under `overlay=:rows`, even
for equivalent results or an explicit reference. Figure titles identify the
point through its compact description. Conductor styles stay consistent across
modes, pages and points. Existing automatic matrix plots and explicit
`overlay=:coordinates` remain available.

Unavailable nonlinear values retain masks and reasons. Selection, tables,
plots, unit conversion, and archive loading use stored observations without
reconstructing a segment or acquiring another representation. Save and load
observations with the existing `.json` or `.jls` archive functions; loaded
selectors identify retained data and are never evaluated.
Reacquire old modal observations to obtain the corrected component identities
and vector coordinates; existing non-modal archives remain supported.

`H` binds to a finite segment only. There is no distance option on `H` for a
modal scan. A phase `H` selector binds both `domain=PhaseDomain` and
`field=:voltage` or `:current`. A phase `Zc` or `Yc` selector binds
`domain=PhaseDomain`. Exact request unit overrides take precedence over
broader selector overrides. Retained rows record mode indices, ordered phase
labels where applicable, source gridpoint ancestry, and relevant upstream and
modal formulation assumptions; unknown assumptions remain unknown.

## Manual broadband study

The nine-cable coaxial research study is a direct, editable Julia
script under `test/manual/calculations/`. It is excluded from TestItemRunner
discovery and ordinary tests. Parse or inspect it without executing the study.
With LineCableModels and GLMakie in the active REPL project, include it from the
repository root. The existing manual plotting environment includes GLMakie:

```sh
julia --project=test/manual/plotting
```

Then in that REPL:

```julia
include("test/manual/calculations/modal_analysis.jl")
modal
line
plots
diagnostic_figure
```

The script keeps each calculation and live GLMakie figure in a named REPL
variable. Save any retained observation explicitly after inspection if needed. Modal
order, residual coupling, and conditioning are research diagnostics, not
acceptance thresholds.
