# Gmsh/GetDP optimization acceptance

All four retained 18 kV campaigns pass the full numerical comparison: 404
frequency slices over 0.1 Hz–1 MHz, comprising 98,172 primitive Z/P and final Y
entries. The refactor preserves the reference equations and results. This is
an equivalence check on identical meshes; it does not establish mesh convergence
for every physical observable.

## Numerical results

Errors below are the maximum relative Frobenius norm of any frequency slice,
against the corresponding retained campaign. Reduced Z/P and primitive Y were
also checked through the existing Engine reduction/inversion path. Every
comparison passed the `1e-8` relative acceptance threshold.

| Geometry and material profile | Retained run | Z error | P error | Final Y error | Solve and comparison time |
|---|---|---:|---:|---:|---:|
| Homogenized, lossless default | `run-YiqIRh` | 3.37e-13 | 6.07e-12 | 3.17e-12 | 108.7 s |
| Homogenized, Ametani2004 | `run-aZVWtq` | 8.30e-13 | 9.25e-12 | 2.65e-9 | 129.7 s |
| Detailed, Ametani2004 | `run-Pbj2uz` | 1.88e-12 | 8.98e-10 | 1.97e-10 | 553.9 s |
| Detailed, lossless default | `run-orkY1H` | 1.30e-12 | 1.08e-12 | 1.24e-12 | 412.3 s |

The freshly resolved material data match each historical `model_data.pro`
byte for byte. All 101 frequency-specific mesh plans match, and the actual
retained meshes are used read-only. Historical raw checksums were verified
unchanged after the comparisons. The validation uses the production batching,
process coordinator, column checks, aggregate parser and Engine reductions.

The largest P difference occurs at frequency index 14 in the detailed lossy
case. Repeating the unchanged single-column solver with the same one-thread
settings produces **exactly the same 162 Z/P entries as the optimized batch**,
including that difference from the historical result. It is therefore not
caused by factor reuse. The largest Y sensitivity occurs in the homogenized
lossy case, whose reduced P condition number reaches 4.80e8; its largest Y
difference is 2.65e-9 relative, or 0.000000265%.

Independent small-coax checks retain their analytical acceptance limits:
capacitance within 2% and lossy conductance within 3%. The frozen Python
reference comparisons also pass their existing primitive/reduced checks.

## Performance and retained evidence

The full comparison used four frequency workers and one BLAS/OpenMP thread
per GetDP process. The four runs required 1,204.5 s (20.1 min) in total. The
retained campaign timings sum to 22,086.8 s (6.14 h), an observed ratio of
18.3. These are historical campaign and current frozen-mesh timings, with
different background load; they are not a controlled claim about mesh
generation speed. The new timings include coordinator/comparison overhead
and exclude mesh generation. Fresh meshing is separately exercised by the
public small-case integration tests.

PETSc reports exactly **404 symbolic factorizations, 404 numeric factorizations
and 3,636 solves**. Every fresh frequency factors once and solves nine right-hand
sides. The historical orchestration required 3,636 numerical process launches
and factorizations for the same four scans.

The [machine-readable evidence](gmsh-performance-validation.json) contains
per-campaign mesh hashes, numerical errors, inversion diagnostics, source/run
locations, summed worker phase timings and PETSc event counts. Worker phase
times are sums across concurrent processes and must not be added to wall time.
The earlier [controlled probes](gmsh-performance-evidence.json) measured about
5.8× from batching and 1.98× from two concurrent frequency batches.

Reproduce the full comparison from the repository root:

```sh
JULIA_LOAD_PATH=@:.:@stdlib julia --project=test/gauntlet test/gauntlet/fem_performance_validation.jl
```

Set `LINECABLEMODELS_GETDP` only to override the package-owned GetDP artifact.
The driver creates a fresh `.linecablemodels/fem/validation-*` root;
it never resumes or overwrites a historical campaign. The companion
`test/gauntlet/fem_worker_scaling.jl` checks one, two, four and eight workers on
eight retained detailed-cable meshes, requiring exact serial agreement and
one factorization per frequency.

The eight-frequency scaling check produced identical Z/P values at every
worker count (1,296 primitive entries per execution):

| Workers | Wall time | Speed relative to one worker | Peak sampled aggregate solver RSS |
|---:|---:|---:|---:|
| 1 | 108.6 s | 1.00× | 0.92 GiB |
| 2 | 53.2 s | 2.04× | 1.83 GiB |
| 4 | 28.5 s | 3.81× | 3.65 GiB |
| 8 | 17.9 s | 6.08× | 7.13 GiB |

RSS is sampled from each active child process on Linux. The default remains
two workers to bound memory across unknown cable sizes. Four workers are a
validated choice for these 18 kV cases; eight improve throughput further when
roughly 8 GiB of solver memory plus coordinator/system headroom are available.

## Operational validation

The headless extension checks pass 905 assertions covering meshing,
enclosures, parsing, failures, recovery, maps, analytical and legacy numerical
references. The existing real GUI integration passes 61 assertions. Core
formulation/options/display checks pass 654 assertions, including the optional
uncertainty display test in an environment exposing its Distributions dependency.

Recovery tests cover corrupt checkpoints, complete attempt adoption, retries
starting at terminal 2, serial/parallel agreement, changing worker count,
refusing changed solver thread settings, cancellation, process reaping, and
run ownership. A malformed historical attempt manifest is ignored without
discarding other valid columns. Additional checks verify cross-process lock
exclusion and an actionable capability error for an unsupported solver.

Four isolated real UI scenarios pass: closure before mesh generation, closure
before solving, cancellation during active solving, and successful display of
all 36 expected field maps for a two-frequency/two-terminal case. Caller-owned
models, parameters and views are restored. The event pump checks availability
before `fltk.wait`, which can otherwise initialize a GUI again after closure.

Validation platform: Linux, Julia 1.12.7, Gmsh API 4.15.0 and GetDP
3.6.0-git-1cf7fa06 with complex PETSc 3.14.4 and MUMPS. The Windows filesystem
locking branch has not been exercised here. Gmsh calls remain with the session
owner; external code must not mutate the same Gmsh session during a computation.
