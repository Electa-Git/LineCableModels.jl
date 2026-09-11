# MATLAB adaptive trapezoids: reusable sampling ideas

Reviewed 2026-09-10 while CIM performance changes were paused for this review.
The source is `M_avg_mut.m` and `M_avg_self.m` in the user's
`LineCableLab-bkp/Submarine_1Layers_2Cond_in_Sea` directory. Their
`sommerfeld_trapz_auto`, diagnostic helper and evaluation helper are identical.
The routine starts at line 439 in the mutual file and line 418 in the self file.

The adaptive option is `integration='adapt-trapz'`; ordinary `trapz` is a
different implementation that compares grids of N and 2N points.

## Useful ideas

The adaptive routine selects an initial cutoff Λ from exponential damping H,
using L = max(20, −log(0.01 rtol)) and Λ = L/H. It sometimes enlarges Λ using
the larger material root magnitude. The finite mapping is
λ = scale t/(1−t), with scale = Λ/10 and tmax = 10/11.

Two distinct initial-resolution conditions control the cosine and large-λ
geometric envelope:

\[
|y|\Delta\lambda\leq\pi/10,\qquad H\Delta\lambda\leq 0.25.
\]

The maximum mapping derivative converts these conditions into a global sample
count. Rounding to 2^p+1 points makes subsequent doubling of subintervals nested.
Diagnostics separately report phase demand, envelope demand, cutoff, grid size,
iterations and estimated error. These are useful distinctions to preserve.

## Limitations relevant to the Julia implementation

- Refinement doubles the entire grid and reevaluates old nodes. The nesting
  permits sample reuse, but the MATLAB evaluator does not exploit it.
- With the chosen mapping, tmax times the maximum mapping derivative is 11Λ.
  Applying that maximum globally imposes the largest mapped spacing constraint
  on every cell. Local resolution budgets can avoid this over-allocation.
- Grid differences exclude the omitted integral after Λ. The initial
  exponential cutoff is a heuristic, not a bound on the complete integrated
  tail, including amplitudes, denominator behavior and cancellation.
- Only the larger material root affects initial coverage. Smaller roots,
  denominator scales and narrow near-contour branch/pole features need their
  own metadata; global phase and geometric decay alone cannot detect them.
- The difference divided by three estimates composite-trapezoid error in the
  smooth second-order regime. It is not a general convergence guarantee and
  must not be transplanted into the double-exponential backend's estimator.
- Nonfinite samples are silently replaced by zero. The maximum-grid condition
  warns and returns a result even if the tolerance remains unmet; an initial
  count above maxN is evaluated before that condition is checked.

## Numerical check and provenance

MATLAB startup failed with `Unable to load ApplicationService for command
client-v1`; the process was interrupted. No successful MATLAB execution is
claimed. A Python transcription of the adaptive sampler was run against two
analytic integrals, both with H=2, rtol=1e-8 and atol=1e-12:

| Integrand | Analytic integral | Returned value | Estimated error | Final N |
|---|---:|---:|---:|---:|
| exp(−2λ) cos(λ) | 0.4 | 0.4000000013928689 | 1.4300e-9 | 4097 |
| exp(−((λ−0.371)/1e-7)²) exp(−2λ) | 8.4397276229e-8 | 0 | 0 | 2049 |

The second is a synthetic coverage test, not an earth-matrix result. Both
nested grids miss the narrow feature, so agreement alone reports convergence.
The current Julia feature-coverage tests already guard against this class of
failure when the feature is declared.

The transcription and output are retained in
`.linecablemodels/qa/matlab-adaptive-trapz/` as `audit.py` and `python.log`.

## Incorporation into the pending CIM work

The reusable addition is an explicit, local resolution budget for phase and
envelope variation, alongside the existing physical features and separate
tail/error budget. Evaluate these conditions in the actual contour coordinate,
including radial-root derivatives and Bessel/path terms where applicable.
Retain nominal sampling decisions while preserving uncertain physical values.

For trapz, these conditions can guide panel subdivision and diagnostics; the
qualified DE package continues to own its rule and local refinement. For CIM,
they guide fitting regions and physical-contour validation. Pencil samples
must remain uniformly spaced in the coordinate of the exponential expansion:
uniform mapped-t samples cannot simply be fed to the existing pencil.

This complements Rallis's compact fitting regions. It does not establish a
universal fixed point count, replace fit certification, or justify restoring
the many overlapping SVD windows that dominate the current CIM cost.
