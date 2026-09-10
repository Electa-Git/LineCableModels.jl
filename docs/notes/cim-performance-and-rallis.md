# CIM cost and the Rallis thesis

The implementation follow-up is recorded in
[`spectral-backend-optimization.md`](spectral-backend-optimization.md).

Measured on 2026-09-09. This investigation profiles the current complete earth
closure; it does not change the physical formulation or production integration
implementation.

## Measured cause

The accepted three-wire geometry at 1 MHz uses radius 0.0425 m, depth 1 m,
horizontal positions 0/1/2 m, earth resistivity 0.1 Ω m and Γ = 0. The complete
Ze/Pe/Ye preparation uses `rtol=1e-6`, Julia 1.12.7 and one BLAS thread.
Compilation is excluded and the physical workspace is reused.

An uninstrumented run took **1.665 s and allocated 1,198,725,432 bytes**
(minimum of three warmed calls). The instrumented call took 1.680 s and
allocated 1,201,781,032 bytes. All three resulting matrices agreed with the
uninstrumented results under Julia's default `≈` comparison. The instrumented
copy is created in memory from the production function; production source is
unchanged.

| Stage | Calls | Time |
|---|---:|---:|
| Local matrix-pencil SVD | 1,214 | 599.6 ms |
| Global amplitude SVD | 18 | 577.6 ms |
| Pivoted QR to select image candidates | 9 | 129.4 ms |
| Quadrature of the absolute weighted fit residual | 18 | 138.4 ms |
| Basis construction/selection copies | 27 | 87.1 ms |
| Hankel matrix construction | 2,428 | 26.9 ms |
| Held-out sample checks | 18 | 24.7 ms |
| Independent quadrature of the original integral | 18 | 1.15 ms |
| Analytic image sums | 18 | 1.05 ms |
| True-tail quadrature | 18 | 0.38 ms |

The table omits smaller stages. SVD and QR account for about **78%** of the
instrumented elapsed time. The two SVD stages alone allocate 819 MB. Local
Hankel matrices are 64×65; the largest amplitude factorization is 1282×192.
All 18 fits succeed at CIM refinement zero: the large cost is already present
in the initial attempt, rather than being caused by repeated failures.

There are nine ordered interactions and two integrated kernel families. Final
representations contain 158–192 images. Replaying the same 18 accepted radial
image sums in an allocation-light scalar loop took **0.692 ms and 288 bytes**
(minimum of 20 warmed calls). This replay excludes construction, verification,
error propagation and matrix assembly; it is evidence about image evaluation
cost, not a benchmark of an implemented cached solver.

The prior complete production benchmark measured direct quadrature at 1.89 ms
for this three-wire case. Consequently, making image evaluation cheap cannot
by itself make the current end-to-end CIM implementation competitive.

## Why the implementation creates so much work

`src/engine/integration.jl:cim_estimate` clears the image buffers at entry and
constructs a new representation for each integral. Those buffers are scratch
storage, not a reusable fit cache. For each adjacent feature interval it runs a
pencil; it also runs overlapping pencils from zero to successive interval
endpoints. `src/engine/spectralsampling.jl:cim_windows` bridges separated scales
by factors of four. In this case that produces 1,286 sampled windows across
the 18 fits, of which 1,214 require SVD.

Candidates from these windows feed a large global weighted amplitude fit.
When necessary a pivoted QR caps the candidate count at 192. The amplitude
problem is then solved with another dense SVD. This successfully covers
separated scales, but it multiplies factorization work and retains many images.
An independently integrated original kernel and an integrated absolute fit
residual certify every result. The latter is substantially more expensive
because each residual evaluation sums the image expansion.

For this equal-radius, equal-depth case, identical scalar integrations repeat
at separations 0, 1 and 2 m. More strongly, each same-earth radial kernel in
`EarthRadialSpectrum` depends on the medium state and a scalar normalization;
height and separation appear in the analytic weight. Thus two underlying
kernel functions generate the 18 independent fits here. The existing fitting
metric, sample coverage and acceptance budget are geometry-dependent, so
sharing a fit requires validating that shared domain, not merely retaining
whichever pair was fitted first. Reusing identical scalar interactions does
not require imposing symmetry on the manuscript's ordered matrix closure.

## What Rallis does differently

Source: Konstantinos Rallis, *Electromagnetic study of underground conductors*,
[official thesis record](https://phdtheses.ekt.gr/eadd/handle/10442/34633), and
the supplied `PhD_Ralis.pdf`. Page references below are printed page numbers;
the corresponding PDF page number is one greater.

In Chapter 4, equations (4.5)–(4.9), pp. 70–71, the buried-conductor treatment
fits the material factor

\[
\frac{1}{\lambda+\sqrt{\lambda^2+k^2}}
\approx\sum_n c_n\exp\!\left(s_n\sqrt{\lambda^2+k^2}\right).
\]

It samples uniformly along a finite straight segment in the transformed
variable γ = √(λ²+k²), starting at k. The demonstrated single-level settings
are 100 samples, 14 terms and an endpoint determined by T₀ = 10|k|.
On p. 71 he reports that 8–10 terms suffice in most examples and explicitly
explains that the coefficients are independent of conductor geometry: once
constructed, they can be reused for different depths and separations.
Frequency enters k, so this observation alone does not justify reuse across
arbitrary frequency-dependent medium states.

The two-level method, pp. 71–72 and 83, uses different sampling steps in two
regions: dense sampling where the function changes rapidly, and coarser
sampling in the smoother region. This is much more compact than building
pencils on every local interval plus cumulative overlapping intervals.

There are two important limits to the comparison:

- On p. 75, he reports loss of single-level accuracy at high frequencies and
  large separations. Increasing the sample count increases cost. Single-level
  GPOF is faster than numerical Pollaczek integration in his tests; two-level
  GPOF has comparable cost. The benefit of reuse is reiterated on p. 84.
- The mixed air/earth fit in equation (4.11), p. 75, includes
  exp(−h₂γ)/(λ+γ), so burial depth remains inside that fitted kernel. Even this
  thesis does not provide one geometry-independent fit for every interaction.

His buried transform yields K₁ image expressions, whereas our existing radial
decomposition uses a different residual and K₀ image expressions. The fitting
architecture is transferable; his coefficients and fixed sample counts are
not a validated replacement for the complete manuscript's Ze/Pe/Ye kernels,
prescribed Γ, unequal permeabilities and voltage-path terms.

## Recommended implementation direction

1. Separate a typed image representation from the per-integral scratch
   workspace. Key it by the complete material/frequency/Γ state, kernel family,
   branch/contour and normalization. Specify the geometry range and error
   budget for which it has been verified. Reuse identical scalar interactions
   immediately; share material-only fits where the exact decomposition allows.
2. Keep the physical feature description reusable across backends, but let CIM
   use it to choose a small adaptive set of fitting regions. A feature point
   need not trigger both a local SVD and a cumulative SVD. Test compact
   one-/two-region fits first; add regions or rank only where the residual
   requires them. Exact extraction of constant/asymptotic terms is another
   candidate for reducing rank without changing the physical formula.
3. Certify the representation over its declared geometry range during
   construction. Cache hits must have a valid propagated error estimate;
   outside-range uses trigger validation/refinement. Keep independent
   quadrature as an audit and construction check rather than repeat both full
   original and residual integrals for every already-certified use. A finite
   set of representative geometries alone is not a proof of an envelope bound.
4. Compress the accepted expansion and reuse factorization buffers where
   practical. Benchmark fit construction, certification, reuse and complete
   matrix assembly separately, retaining the frozen physical accuracy suite.

The desired speedup comes from fewer factorizations, fewer images and repeated
use of each verified fit. Removing quadrature checks alone would leave the
dominant cost. A single previously unseen interaction can still favor direct
quadrature, even with a well-designed CIM backend.

## Reproduction evidence

Local artifacts are in `.linecablemodels/qa/cim-performance/`: `profile.jl`,
`summary.toml`, `stages.csv` and `run.log`. The profiling script writes its
results to `/tmp/lcm-cim-profile`; create that directory before reproducing.
It creates an instrumented function copy and a process-local dispatch override,
so run it in a fresh Julia process from the repository root:

```sh
mkdir -p /tmp/lcm-cim-profile
OPENBLAS_NUM_THREADS=1 JULIA_DEPOT_PATH=/tmp/lcm-julia-depot:/home/amartins/.julia julia --startup-file=no --compiled-modules=existing --project=. .linecablemodels/qa/cim-performance/profile.jl
```

Profiled `src/engine/integration.jl` SHA-256:
`f456c9708fc600e2d24b8c421a846462643465de52a7a4641094193a7b2bbed7`.
The script checks that this source remains unchanged. Timings include minor
instrumentation overhead and ordinary run-to-run variation; they identify the
dominant work rather than predict an optimized implementation's speedup.
