```@meta
EditURL = "../literate/gauntlet.jl"
```

# Gauntlet: validating against PSCAD

Gauntlet is a separate validation application for LineCableModels. It turns
independently harvested PSCAD evidence into typed cases, runs the native
`compute!` path, and retains enough diagnostics to locate a discrepancy by
frequency and matrix entry. The root package does not depend on Gauntlet.

The tracked smoke corpus works offline and contains 16 fixed SHA-derived
records. A full corpus uses the same objects and functions.

````@example gauntlet
using Gauntlet
using LineCableModels
using DataFrames
````

## Open and inspect a corpus

`Corpus` is indexed lazily: opening it reads only the small index. A case is
materialized from its normalized JLD2 record when requested.

````@example gauntlet
corpus = Corpus(:smoke)
(length = length(corpus), first_ids = first(keys(corpus), 3))
````

Select a sparse coaxial reference. Its `problem` and `formulation` are native
LineCableModels values reconstructed through the normal public modeling API.

````@example gauntlet
case_id = "29b9150df4fe52667d41f6e8d423b5bbaace24ee148a3c523c55070840b2ce03"
case = corpus[case_id]
case
````

Scientific context is explicit. Assumptions explain deliberate
PSCAD-to-LineCableModels approximations; provenance identifies the source,
definition, campaign, hashes, and PSCAD version.

````@example gauntlet
assumptions(case)
source = provenance(case)
(source.version, source.definition, source.case_sha256)
````

## Run a case and inspect matrix diagnostics

`gauntlet(case)` normally runs every check available for that record. Here we
request Z and Y explicitly so the guide remains quick while still calling the
same native `compute!` entry point used by the full suite.

````@example gauntlet
trial = gauntlet(
    case;
    checks = (MatrixCheck{:Z}(), MatrixCheck{:Y}())
);
nothing #hide
````

The case is classified `Approximate`, so native/reference discrepancies are
diagnostic rather than gating failures. The numerical record is still
complete: maximum absolute and scaled relative error, a Frobenius residual at
every frequency, the worst frequency and port pair, reciprocity, symmetry,
and positive-real diagnostics.

````@example gauntlet
matrix_diagnostics = DataFrame([(
                                    check = string(typeof(comparison.check)),
                                    verdict = string(nameof(typeof(comparison.verdict))),
                                    max_relative = comparison.metrics.max_rel,
                                    worst_frequency_hz = comparison.metrics.worst_frequency,
                                    worst_entry = comparison.metrics.worst_entry
                                ) for comparison in trial.comparisons])
````

Native and reference values remain ordinary `LineParameters`; the regular
accessors provide the full frequency response at one matrix position.

````@example gauntlet
actual_z11 = Z(trial.actual, 1, 1)
reference_z11 = Z(reference(case).phase, 1, 1)
(actual_z11, reference_z11)
````

## Inspect modal evidence and evaluate a PSCAD vector fit

A detailed smoke record contains the calculated modal transform, propagation
quantities, characteristic admittance, and fitted channels over frequency.
Mode comparison first performs one-to-one correlation assignment, phase
alignment, and degenerate-subspace alignment; incidental order and complex
phase therefore do not become false discrepancies.

````@example gauntlet
detailed_id = "07a34a318668e01782dfbaff5f814dc15ce5b7e6ec74c5e88782bca39feefd5f"
detailed = corpus[detailed_id]
modal = reference(detailed).modes
(
    frequencies = length(modal.frequency),
    transform_size = size(modal.transform),
    propagation_size = size(modal.propagation)
)
````

Gauntlet reads `.clo` and `.tlo` poles, residues, constants, and delays but
does not fit a new time-domain model. `evaluate` calculates the imported model
at any frequencies inside its recorded range.

````@example gauntlet
fit = reference(detailed).fit
fit_frequencies = exp.(range(
    log(first(fit.frequency_range)),
    log(last(fit.frequency_range));
    length = 3
))
fit_values = evaluate(fit, fit_frequencies)
(
    characteristic = size(fit_values.characteristic),
    propagation = size(fit_values.propagation)
)
````

## Run a suite and write reports

A suite fixes case identifiers, checks, fidelity policy, tolerances, and the
performance selection. Rejected PSCAD cases are preserved as searchable
source evidence; they never masquerade as Julia solver failures.

````@example gauntlet
rejected = only(value for value in corpus if value.fidelity isa Rejected)
suite = Suite(:guide; corpus, ids = [rejected.id])
report = gauntlet(suite)
DataFrame(report)
````

Reports are Tables-compatible and have explicit Markdown, JSON, and CSV
writers. The documentation uses a temporary file; normal applications choose
a persistent report directory.

````@example gauntlet
report_path = tempname() * ".md"
write_report(report_path, report)
first(split(read(report_path, String), '\n'), 8)
````

The smoke corpus currently contains no native reconstruction promoted to
`Exact`: formulation equivalence has not yet been demonstrated. This is a
scientific status, not a missing result. Approximate cases retain all
diagnostics, while exact cases alone gate frozen numerical tolerances.

````@example gauntlet
(
    exact = count(value -> value.fidelity isa Exact, corpus),
    approximate = count(value -> value.fidelity isa Approximate, corpus),
    rejected = count(value -> value.fidelity isa Rejected, corpus)
)
````

## Read performance evidence

The Tutorial 3 EMT baseline is a typed `Perf` record captured after two
warm-ups with ten BenchmarkTools samples and `evals=1`. Bytes and allocations
have a five-percent regression gate; wall time remains an informational trend.

````@example gauntlet
baseline = performance_baseline("tutorial3-emt")
(
    minimum_ms = 1e3 * baseline.minimum_seconds,
    median_ms = 1e3 * baseline.median_seconds,
    bytes = baseline.bytes,
    allocations = baseline.allocations,
    frequencies = baseline.frequencies
)
````

## Use a local full artifact

The full suite is explicit and network-free once an artifact has been
extracted. It uses the same grammar as the smoke examples:

```julia
full = Corpus("/path/to/pscad-normalized-v1")
full_case = full["canonical-case-sha256"]
full_trial = gauntlet(full_case)
full_report = gauntlet(Suite(:full; corpus = full))
write_report("gauntlet-full.json", full_report)
```

The command-line equivalent is:

```bash
julia --project=gauntlet gauntlet/bin/run.jl full \
    --artifact-dir /path/to/pscad-normalized-v1
```
