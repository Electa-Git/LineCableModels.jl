# Modal analysis and line segments

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
coefficients are `Tv \ (Zphase * Ti)` and `Ti \ (Yphase * Tv)`. Inverse conversion
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

The native, FEM and PSCAD line-parameter calls and the parametric studies accept the
`modal` keyword. The positional form takes a `ModalAnalysisFormulation` and composes any
formulation whose `compute(problem, formulation; options)` returns `LineParameters`. A
user-owned formulation composes too:

```julia
modal = compute(line_problem, user_formulation, ModalAnalysisFormulation(:default))
```

A problem `Gridspace`, one of its grid points and a `ParametricProblem` compose in the
same way. Cable constants are not line parameters and do not compose with modal analysis.

`modal_options` controls the modal action only. A nonempty value requires
`modal`. The upstream solver, including its selected backend and batched path,
completes before modal analysis begins. Upstream `on_result`, progress, and
timing belong to that upstream call. Modal `on_result`, progress, and timing
belong to the downstream call. Timed results preserve the source gridpoint
identity. A missed iteration or coupling target produces a warning and a
diagnostic record. Finite results are returned. The iteration convergence
setting is a numerical target. It does not bound the physical error. Inspect
`details(modal).data.modal.diagnostics` for fallback and missed frequency
indices, Z/Y residual coupling, returned eigen-residuals, iteration counts,
and convergence flags. These diagnostics describe the numerical solve. Scientific
acceptance requires a separate decision. The scalar characteristic and response construction
uses the selected modal approximation when coupling remains.

The built-in formula selects paired voltage and current eigenvectors and
retains a forward root with nonnegative real part, or nonnegative imaginary
part for a purely imaginary root. The modal workspace supplies impedance and
admittance slices, coordinate conversion buffers, the admittance-impedance
product, eigenpair history and voltage-vector buffers. The built-in formula
tracks eigenvalues of its normalized, shifted eigenproblem. Multiplying
`eigenvalue + 1` by the scale gives an eigenvalue of the physical product,
distinct from a propagation root. The selected formula builds the buffers of its
normalized eigenproblem, Hungarian assignment and Levenberg–Marquardt least squares
in its `initialize_buffers` method. Its arithmetic remains in
`decompose!(::Val{:chrysochos2014}, ...)`. `Tv`
and `Ti` map modal coordinates to phase coordinates. Total source coefficients
are normalized by their declared source length before the line segment is
bound.

## Alternative complex LM implementation

Select `:vieira2026` for the complex Levenberg–Marquardt implementation adapted
from the supplied Vieira tracker, based on Chrysochos et al. (2014):

```julia
selected = ModalAnalysisFormulation(:vieira2026)
modal = compute(ModalAnalysisProblem(phase), selected)

# The convenience call performs the same two computations.
modal = compute(line_problem, line_formulation; modal=:vieira2026)
```

The method consumes the completed matrices at the stored frequencies. It uses
a two-sample predictor, SVD minimum-norm steps, greedy eigenvalue assignment and
clustered eigenvector refinement. When tracking misses its matching or predictor
target, it uses a correlation-matched eigensolution at that same frequency.
There is no interpolation, adaptive frequency insertion or physical reevaluation.
The default formula remains `:chrysochos2014`.

Numerical controls belong to the selected formula:

```julia
selected = ModalAnalysisFormulation(:vieira2026; options=(
    iteration=(convergence=1e-13, max_iterations=100),
    tracking=(predictor_tolerance=0.25, eigenvalue_tolerance=1e-5,
        cluster_tolerance=1e-3, order_by_velocity=true),
))
```

`convergence` is the Euclidean norm target for the complex eigenpair residual,
including its bilinear normalization constraint. It is not the differently
interpreted convergence control of `:chrysochos2014`, nor an accuracy guarantee.
Matching compares eigenvalue distances against `eigenvalue_tolerance` times the
larger of one and the largest normalized eigenvalue magnitude. Clusters connect
normalized eigenvalues whose gap is below `cluster_tolerance` times the larger
of their shifted magnitudes and `1e-3`. `predictor_tolerance` bounds relative
eigenvalue and isolated-vector changes from the predictor. These decisions select
the numerical tracking route. Finite results are retained with warning diagnostics.

Tracking history uses `transpose(t)*t = 1`. Published current bases have unit
Euclidean column norms. Voltage bases are the normalized paired `Z*Ti` columns.
With `order_by_velocity=true`, decreasing phase constant at the final frequency
sets one mode order for the entire scan, including its diagnostics. Disable this
option to retain the seed's tracked order. The first frequency is obtained by
direct eigendecomposition and has no LM iteration count or convergence flag.
Subsequent flags describe LM convergence even when a conventional eigensolution
is ultimately returned. The reported eigen-residual describes that returned pair.

This adaptation uses `s=j2πf` and the package's vacuum permittivity
`8.8541878128e-12 F/m`. The supplied `EigTrack.jl` used `8.854187817e-12 F/m` and
also supported complex-frequency reevaluation. The formula builds its complex
least-squares, predictor and assignment buffers in its `initialize_buffers` method and
reuses the common modal buffers. The detailed method attribution is attached to its
formula description.

## Newton–Raphson implementation

Select `:wedepohl1996` for the Newton–Raphson eigenpair continuation adapted
from `eig_newton` in the supplied UniversalLineModel implementation of
Wedepohl, Nguyen and Irwin (1996), DOI: 10.1109/59.535695:

```julia
selected = ModalAnalysisFormulation(:wedepohl1996; options=(
    iteration=(convergence=1e-9, max_iterations=60),))
modal = compute(ModalAnalysisProblem(phase), selected)
```

The seed modes are ordered by decreasing attenuation. Subsequent samples use
the previous eigenpairs in Newton corrections of `YZ / opnorm(YZ, 2)` (the matrix
spectral norm), with the bilinear constraint `transpose(t)*t = 1`. The
convergence target bounds the largest absolute correction in this normalized
problem. A singular iteration, missed target or duplicate eigenpair causes
same-frequency direct eigendecomposition and greedy correlation matching.
These events remain warning diagnostics. No physical samples are reevaluated.

## Shared basis rotation and sign continuity

All modal formulations, including user-defined formulations, use
`options=(rotate=true,)` by default. After decomposition and before computing
modal Z/Y, each Ti column is multiplied by a unit complex factor that minimizes
its squared imaginary norm. The same factor multiplies its Tv column, preserving
the voltage and current pairing. This does not independently minimize the imaginary
part of Tv. The remaining 180-degree ambiguity is resolved by requiring a
nonnegative real overlap with the previous frequency's Ti column. Mode order
and the decomposition's tracking history are unchanged.

```julia
modal = compute(ModalAnalysisProblem(phase), selected; options=(rotate=false,))
modal = compute(line_problem, line_formulation; modal=:wedepohl1996,
    modal_options=(rotate=false,))
```

`rotate=false` retains the formula's complex column phases apart from fixing
180-degree flips. It does not disable the formula's internal normalization or
tracking. The choice is retained at `details(modal).data.modal.rotate` and in
observation assumptions. Rotated and unrotated observations remain distinct
when grouped.

Column norms, diagonal modal Z/Y/Zc/Yc, propagation constants and reconstructed
phase quantities (including H) remain unchanged under the paired operation. Residual
off-diagonal modal entries retain their magnitudes but can change phase with the
coordinate choice. No values are clipped or set to zero by this operation.

## Finite segments and source length

`PropagationParameters` binds a modal scan to one physical segment. A source
with total coefficients retains its positive source normalization length.
The line segment may have another length. For an external scan without a
known source length, supply the target length explicitly.

```julia
segment = PropagationParameters(modal; line_length=150.0)
next_segment = PropagationParameters(segment; line_length=300.0)
length_space = PropagationParameters(modal;
    line_length=Grid((150.0, 300.0))) # lazy line segment Gridspace

H(segment)                              # mode × frequency
alpha(segment)                          # attenuation constant, Np/m
beta(segment)                           # phase constant, rad/m
velocity(segment)                       # phase velocity, m/s
H(segment, PhaseDomain; field=:voltage) # ordered phase voltage response
H(segment, PhaseDomain; field=:current) # ordered phase current response
```

The modal forward factors are `exp.(-gamma(segment) .* line_length(segment))`.
`alpha` and `beta` are the real and imaginary parts of the per-unit-length
`gamma`, independent of segment length. `velocity` is
`2πf / beta` at each retained mode and frequency. Undefined values remain
unavailable in observations, with their masks and reasons retained. The default
display units for the first two are Np/km and rad/km.
The phase voltage and current responses use `Tv` and `Ti`, respectively.
Zero length gives identity propagation. Rebinding preserves the embedded source.
A changed length receives a new gridpoint identity. Indexing a
modal or finite result by frequency slices coefficients, bases, roots, and
frequency-aligned diagnostics together without a new solve.

## Collections and observations

Completed phase result spaces transport through
`Gridspace{ModalAnalysisProblem}(phase_results)`. One modal formulation
keeps one downstream result per source point. A formulation Gridspace forms
the ordinary product or zip selection. Completed modal results transport
through `Gridspace{PropagationParameters}(modal_results)` for finite lengths.
Scalar quantities compose with broadcasting, such as `gamma.(modal_results)`
and `H.(segment_results)`. Each tensor remains one value.

`ObservedResult` retains detached quantities and the existing four sections:
gridpoint, quantities, errors, and timings. Requests can mix phase coefficients
with `gamma`, `Zc`, `Yc`, `Tv`, and `Ti`. A line segment also offers `H`,
`alpha`, `beta`, and `velocity`. Bare complex requests acquire complete real and
imaginary pairs. A raw plotting convenience also completes one component
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
rows, mode columns, and frequency samples. Retained coordinates keep both
label sets. Phase matrices have ordered phase rows and columns. `abs`, `angle`,
`real`, and `imag` can be requested through the existing request tuples.
`@observe alpha[:, :]` and `@observe Tv[:, :, :]` keep original
indices. `@observe (H, abs)[:, :]` requests the explicit magnitude component.
Complete magnitude and angle pairs are acquired separately from the default
Cartesian components. Each component has its own physical name, unit, table,
and plot family. A vector component table has frequency rows and one column
per selected mode, named `Mode i`, where `i` is the mode index. Full transformation
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
modes. Legends name conductors. Selected modes fill each page left to right and
then top to bottom. 18 modes produce two pages per component with
`layout=(3,3)`. Every selected coefficient remains plotted.

For a collection, acquire the same requests at each point and pass the retained
results to the same plotting call:

```julia
requests = ((Tv, abs), (Tv, angle), (Ti, abs), (Ti, angle))
transform_space = observables.(segment_results, Ref(requests))
plot(transform_space; overlay=:rows, layout=(3,3))
```

With `overlay=:rows`, each gridpoint starts its own figures with the same block
layout, including equivalent results and an explicit reference. Figure titles
identify the point through its compact description. Conductor styles stay
consistent across modes, pages, and points.

Unavailable nonlinear values retain masks and reasons. Selection, tables,
plots, unit conversion, and archive loading use stored observations. Save and
load observations with the `.json` or `.jls` archive functions. Loaded selectors
identify the retained quantities and coordinates.

`H` uses the length bound to a line segment. A phase `H` selector binds both
`domain=PhaseDomain` and
`field=:voltage` or `:current`. A phase `Zc` or `Yc` selector binds
`domain=PhaseDomain`. Exact request unit overrides take precedence over
broader selector overrides. The retained rows contain indices of modes, ordered phase
labels where applicable, source gridpoint ancestry, and relevant upstream and
modal formulation assumptions. Unknown assumptions remain unknown.

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
