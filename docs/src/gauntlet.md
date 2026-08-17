```@meta
EditURL = "../literate/gauntlet.jl"
```

# Gauntlet datasets and datasource adapters

Gauntlet is a separate, datasource-agnostic validation application for
LineCableModels. It turns independently produced evidence into typed cases,
runs the native `compute!` path, and retains enough diagnostics to locate a
discrepancy by frequency and matrix entry. The root package does not depend
on Gauntlet.

PSCAD is the first datasource implemented end to end. Its tracked smoke
dataset works offline and contains 16 fixed SHA-derived records. A full PSCAD
dataset and datasets from future datasources use the same owned objects and
verbs.

````@example gauntlet
using Gauntlet
using LineCableModels
using DataFrames
````

## Datasource dispatch

A datasource owns external parsing, reference decoding, native case loading,
and ingestion. Symbols resolve through `Val` specializations,
following the same style as the package's external-format exports.

````@example gauntlet
registered = (pscad = datasource(:pscad), fem = datasource(:fem))
````

`:pscad` supports `ingest`, `load`, and `decode`. `:fem` is deliberately
registered so the vocabulary is stable, but those actions currently throw
`FEM datasource is not implemented`. A new implementation extends the same
short verbs for its own record types:

```julia
struct SolverData <: Datasource end
Gauntlet.datasource(::Val{:solverdata}) = SolverData()
function Gauntlet.ingest(::SolverData, source, destination)
    # Normalize the external dataset.
end
function Gauntlet.load(::SolverData, record)
    # Reconstruct one native `Case` from a normalized record.
end
function Gauntlet.decode(::SolverData, record)
    # Decode its normalized evidence into a `Reference`.
end
```

`Dataset(path)` reads the required `datasource` key from `index.toml` and
routes case loading to that datasource. The dataset layer never interprets
`.cli`, `.tli`, `.clo`, `.tlo`, finite-element files, or any other external
format itself.

## Open and inspect a dataset

`Dataset` is indexed lazily: opening it reads only the small index. A case is
materialized from its datasource-owned normalized record when requested.

````@example gauntlet
dataset = Dataset(:smoke)
(length = length(dataset), first_ids = first(keys(dataset), 3))
````

Select a sparse coaxial reference. Its `problem` and `formulation` are native
LineCableModels values reconstructed through the normal public modeling API.

````@example gauntlet
case_id = "29b9150df4fe52667d41f6e8d423b5bbaace24ee148a3c523c55070840b2ce03"
case = dataset[case_id]
case
````

Scientific context is explicit. Assumptions explain deliberate
PSCAD-to-LineCableModels approximations; provenance identifies the datasource,
definition, campaign, hashes, and PSCAD version.

````@example gauntlet
assumptions(case)
origin = provenance(case)
(origin.datasource, origin.version, origin.definition, origin.case_sha256)
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
detailed = dataset[detailed_id]
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
rejected = only(value for value in dataset if value.fidelity isa Rejected)
suite = Suite(:guide; dataset, ids = [rejected.id])
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

The smoke dataset currently contains no native reconstruction promoted to
`Exact`: formulation equivalence has not yet been demonstrated. This is a
scientific status, not a missing result. Approximate cases retain all
diagnostics, while exact cases alone gate frozen numerical tolerances.

````@example gauntlet
(
    exact = count(value -> value.fidelity isa Exact, dataset),
    approximate = count(value -> value.fidelity isa Approximate, dataset),
    rejected = count(value -> value.fidelity isa Rejected, dataset)
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
full = Dataset("/path/to/pscad-normalized-v1")
full_case = full["canonical-case-sha256"]
full_trial = gauntlet(full_case)
full_report = gauntlet(Suite(:full; dataset = full))
write_report("gauntlet-full.json", full_report)
```

The command-line equivalent is:

```bash
julia --project=gauntlet gauntlet/bin/run.jl full \
    --artifact-dir /path/to/pscad-normalized-v1
```
