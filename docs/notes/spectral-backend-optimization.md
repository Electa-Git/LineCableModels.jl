# Spectral backend optimization — 2026-09-10

This implements the performance follow-up to
[`cim-performance-and-rallis.md`](cim-performance-and-rallis.md) and the sampling
review in [`matlab-adaptive-trapz-review.md`](matlab-adaptive-trapz-review.md).
The manuscript kernels, current closure, voltage references and author
registrations are unchanged.

## Reusable complex images

`src/engine/compleximages.jl` defines typed workspace-owned image fits and
certificates. Built-in kernel identities copy s, Γ, material admittivities,
permeabilities and both medium roots, plus the ordered media, kernel family
and every geometric parameter occurring inside that kernel. They do not
retain aliases to mutable material/current-map arrays or rely on hash equality.
The radial fitting origin and physical scalar precision must also match.

Same-earth radial kernels admit a separate scalar exponential normalization;
the image values and their error estimates scale together. Geometry-dependent
mixed/path kernels retain their geometry in the identity. Arbitrary callbacks
without an identity method are not cached. Both uncertain kernels and uncertain
analytic weights are explicitly rejected by CIM; trapz retains their correlated
physical values.

Each representation stores either point certificates or geometry-envelope
certificates. For a radial residual F−F_N, the latter estimates

\[
E=\int_0^\infty |F(u(t))-F_N(u(t))|
\frac{\exp[-h_\min\Re u(t)+y_\max|\Im\lambda(t)|]}{|u(t)|}
\left|\frac{d\lambda}{dt}\right|dt.
\]

The same estimate covers h ≥ h_min and y ≤ y_max on that contour, since
Re(u) ≥ 0 and |cos(yλ)| ≤ exp(y|Im λ|). Cosine and Bessel/cosine weights use
the corresponding exponential envelope, with y+r for the latter. Kernel
identity, contour admissibility and convergence of the image transforms remain
prerequisites. A finite list of sampled geometries is not used as an envelope
certificate. These are numerical estimates, not rigorous analytic bounds:
quadrature error and rounding are included in the reported estimates.

Construction still checks the original integral and the complete weighted
residual. A new geometry can require a new full-range residual certificate or
a new fit. A valid cache hit evaluates the analytic images and propagates the
stored error; it performs no quadrature. The cache holds at most 64 recently
constructed fits, each with at most 32 certificates after extensions.
Each integration request still evaluates one kernel prototype to check the
physical scalar type; this is included in the evaluation counter on a cache hit.

The complete earth assembler visits lower-height, larger-separation interactions
first and uses a common admissible angle cap for same-earth radial terms. This
lets demanding weights seed fits reusable by nearby interactions. It does not
impose symmetry on the ordered K/H/L matrices. In the equal-depth three-wire
1 MHz example, two shared fits serve all 18 scalar integration requests.

Construction starts with at most eight fitting regions and a 96-image initial
cap, bounded by the caller's `max_terms`. Verification failures expand the
region set before increasing uniform pencil resolution. Mandatory material
features also enter amplitude fitting. These are initial work budgets, not
fixed physical cutoffs or accuracy assumptions. The established full coverage
and supplemental exponent construction remain available within refinement
limits. Sample/Hankel storage is reused; Hankel SVD operates in place. Held-out
checks and analytic image evaluation use scalar loops rather than allocating
temporary exponential arrays at every sample.

## Local transformed trapezoids

The DE package remains the numerical rule. Phase and envelope changes of the
analytic weight guide finite panel seeding in the actual contour coordinate,
including the radial root transformation. Existing feature breakpoints remain
mandatory. The seeding horizon never truncates the integral: the remaining
interval to infinity is integrated and verified.

Independent panel quadrature now supplies both normalization and verification.
The redundant separate whole-integral reference is removed. Reference errors
are tightened when their sum exceeds the allocated budget, and local absolute
discrepancies cannot cancel between panels. A panel's reference storage is then
reused for its verified DE result. Successful panels remain intact while the
largest-error panel is subdivided; the package reuses nested DE samples within
each panel. Sampling and QuadGK segment buffers are reused.

The reference and refinement functions explicitly specialize on their callable
type. Without that declaration Julia's forwarding-function specialization
heuristic introduced dynamic calls and per-panel allocations despite inferred
test probes appearing concrete. Resolution diagnostics live outside the
homogeneous vector-buffer tuple; mixing them into that tuple similarly broke
the existing allocation behavior. Quadrature retains its 784-byte complete
preparation benchmark.

## Validation and reproducible measurements

The regression suite covers frozen two-/three-wire matrices at 0.1 through
1 MHz, mixed layers and prescribed Γ, deep/interface/finite/scalar voltage
references, physical scalar types, public uncertainty derivatives, formula
ownership and author registrations. New tests cover fit identity invalidation,
geometry-envelope reuse, normalization, inferred result types, mutable
uncached callbacks, oscillatory cancellation and uncertain analytic weights.

`test/gauntlet/spectral_backend_cost.jl` measures complete Ze/Pe/Ye assembly at
1 MHz with one BLAS thread, after compilation. `construction` empties the image
cache before each call while retaining scratch buffers; `reuse` retains the
certified images and still reassembles the matrices. Its JSON records source
hashes, allocations, kernel evaluations, fit/pencil/certificate counts and image
counts, and checks the resulting matrices against the frozen fixture.

The final measured results and test log are retained under
`.linecablemodels/qa/spectral-backend-cost/`. Construction and reuse must be
reported separately: a first CIM fit can remain more expensive than direct
quadrature, even when repeated image evaluation is faster. Trapz's local
sampling controls and reduced allocation do not imply that it always outpaces
the previous rule or direct quadrature.

### Final complete-matrix measurements

The combined run passed **5,150 selected unit checks** and **432 mixed-reference
checks**. This includes 3,888 frozen matrix checks over both geometries, eight
frequencies and all three backends; 85 public uncertainty/sharing checks; and
182 scalar-type/air-root checks. The mixed test includes unequal permeability,
prescribed Γ and deep/interface/finite/scalar voltage references. The benchmark
also asserts agreement with the frozen 1 MHz Ze/Pe/Ye entries. All seven recorded
production source hashes were unchanged through validation and measurement.

Times below are minima of three warmed calls, in milliseconds. Construction
clears only the CIM fit cache; reusable scratch remains allocated. Reuse retains
the accepted representations and still assembles the complete matrices.

| Wires | Backend | Construction time | Reuse time | Construction allocation | Reuse allocation |
|---|---|---:|---:|---:|---:|
| 2 | quad | 0.626 ms | 0.617 ms | 784 B | 784 B |
| 2 | trapz | 2.475 ms | 2.464 ms | 42,448 B | 42,448 B |
| 2 | cim | 22.558 ms | 0.095 ms | 21,260,768 B | 3,504 B |
| 3 | quad | 1.737 ms | 1.783 ms | 784 B | 784 B |
| 3 | trapz | 6.623 ms | 6.440 ms | 84,016 B | 84,016 B |
| 3 | cim | 157.481 ms | 0.467 ms | 102,931,232 B | 6,704 B |

The previous complete-matrix CIM measurements were 772 ms / 528 MB for two
wires and 1,778 ms / 1,199 MB for three wires. Fresh construction is therefore
approximately **34× and 11× faster**, respectively, with about **25× and 12×
less allocation**. These historical timings were taken in a separate run and
remain subject to machine-load variation. Quadrature is still substantially
cheaper for first construction; repeated certified CIM reuse is faster here.

Two fits serve each geometry: 47 and 41 images for two wires; 192 and 45 for
three wires. Three-wire construction now uses 85 local pencils, two accepted
fits and two envelope certificates, versus the previously profiled 1,214
pencils and 18 separate fits. Reuse needs 18 prototype evaluations and analytic
image sums, with no pencils or quadrature. The remaining 192-image response
shows why compact fitting is an adaptive starting budget, not a universal
small image count.

### Same-process trapz comparison

The comparison alternates copies of the previous and current trapz adapters
inside one process, with the same current physical kernels and workspace. Both
adapters are called directly through the same dispatch wrapper. The earlier
`trapz-compare-final.log` used `invoke` only on the current adapter. The direct-call repeat removes
that dispatch asymmetry; the earlier log is retained as superseded evidence.

Each alternating batch records the minimum of five warmed calls. The table
uses the faster of two batches per implementation; both produced matching
Ze/Pe/Ye matrices.

| Wires | Previous time | New time | Previous allocation | New allocation |
|---|---:|---:|---:|---:|
| 2 | 2.254 ms | 2.435 ms | 73,440 B | 51,152 B |
| 3 | 6.116 ms | 6.407 ms | 193,840 B | 103,600 B |

This controlled comparison shows about **30% and 47% less allocation**, with
**8% and 5% higher elapsed time**, respectively. Kernel evaluations changed from
15,664 to 15,920 and from 42,148 to 42,070. Trapz's improvement here is local
sampling control and allocation; there is no measured runtime speedup. The
comparison wrapper has additional allocation compared with the standalone
complete-matrix benchmark above, so compare each table internally. The current
source hash was checked before and after the comparison. Raw results are in
`trapz-compare-direct.log`; machine load still affects timings.
