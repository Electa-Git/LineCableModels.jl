# Scientific consumers

The registered `ichqp-showcase` deck and `cable-study` workbench compose the same
components. They do not import numerical packages. Worker execution and private
terminal availability remain subject to runtime acceptance gates; adding a UI
driver does not enable an unsafe fallback.

## Ownership

| Source | Responsibility |
| --- | --- |
| `src/scientific/StudyCases.jl` | Complete passive specimen inputs and bounded result projections; ordinary dispatch on two concrete case types |
| `src/scientific/ScientificViews.jl` | Existing typed fields, `ScientificJob`, persistent SVG plot and `DataTable` |
| `src/scientific/CableGeometry.jl` | Local proportional radial illustration, independent of solver inputs |
| `src/scientific/StudyRuntime.jl` | Existing worker selectors, explicit preparation and diagnostics |
| `assets/scientific-views.css` | Shared structure; all colours come from the existing brand palette |
| `src/applications/Showcase.jl` | Session-local factories for the deck's five owned live routes |
| `src/workbenches/CableStudy.jl` | Workbench `initialize` / `compose` / `handle!`; no copied scientific or runtime controls |

The numerical implementation remains in `worker/src/operations/`, loaded by its
separate profile environment. UI defaults are explicit curated specimen inputs,
not a second set of numerical defaults. `study_case_validation.jl` checks that
authoritative worker normalization preserves those inputs and their hashes.
Profile registration and input validation import neither numerical engine.
The optional `load_profile!` dispatch hook loads only the profile's fixed
dependencies inside its executor, on explicit preparation/execution. The child
emits `loading_environment` progress; the protected runtime exposes coarse phase,
fraction and elapsed time, not arbitrary child stage/message text. Bootstrap and
readiness remain separate. Package loading uses the declared preparation deadline
rather than the command reader's bootstrap bound.

## Using the consumers

1. Build the publication with `./playground/lcm playground build`.
2. Start the configured runtime gateway using the existing instructions in
   [runtime/CONFIGURATION.md](runtime/CONFIGURATION.md). Do not start another
   listener on an occupied publisher port.
3. From Presentations, launch the scientific showcase; from Workbenches, launch
   CableStudy. Each launch owns an independent UI process and run namespace.
4. Open **Workers and preparation** (the second content slide in the deck).
   Assign the `parameters` / `line-parameters` and `power-flow` / `power-flow`
   roles separately. Explicitly prepare each and inspect its reported evidence.
5. Open a scientific view, edit inputs and choose **Run calculation**. Changing
   fields, slides or display quantities never submits a job. Invalid/incomplete
   fields disable submission. A local edit also disables Run until its canonical
   input draft is acknowledged; Run cannot overtake a Bonito field round trip.
   The shared job control owns this fence, cancellation, receipts,
   last-good provenance and stale-result fencing.
6. Optionally assign the `terminal` / `julia-terminal` role and explicitly
   Connect from the terminal view. It is not the solver's namespace or cache.
   A missing isolation capability is an unavailable feature, not native fallback.

Public/static viewing starts no owned frames. The `bonito` shortcode's
`requires-run="true"` attribute requires a `public-url` fallback. Without a valid
`lcm-run` context, the deck displays that link immediately. Overview/PDF retain
the existing static placeholders, not screenshots or concealed calculations.

## Scientific meaning

Line parameters use two identical coaxial cables in homogeneous earth, with
grounded sheaths. The fixed material/construction inputs are explicitly present
in `StudyCases.inputs(LineParameters())`. Display quantities are Re/Im Z in Ω/m
and Re/Im Y in S/m, for self (1,1) and mutual (1,2) matrix entries. The sample
table uses the same plotted unit. Preparation executes representative frequency
endpoints; it does not claim the requested sweep is already calculated.

The OHL/UGC case uses `ohl_ugc_transition_v1`, earth resistivity 100 Ω m, and the
existing solved/linearized reference network. Evaluations alter passive corridor
lengths while retaining active-device linearization. Curves show nominal and
±length variation at B5, in dB re 1 Ω. These are deterministic sensitivity cases,
not confidence intervals or a new hosting-capacity solver. UGC share excludes
zero-length endpoints to match the worker's regularization explicitly.

The lightweight geometry illustration is intentionally independent of numerical
input state. It changes only radial proportions. This curated version does not
provide persistent projects or import/export files; those are not prerequisites
for repeating the small, fully visible specimen inputs.

## Focused checks

From the repository root:

```sh
julia --startup-file=no --project=playground playground/test/runtests.jl
julia --startup-file=no --project=playground/worker/profiles/line-parameters playground/runtime/test/study_case_validation.jl line-parameters
julia --startup-file=no --project=playground/worker/profiles/power-flow playground/runtime/test/study_case_validation.jl power-flow
LCM_RUNTIME_BROWSER_SUITE=scientific bash playground/runtime/test/run-bonito.sh
bash playground/runtime/test/run-broker.sh scientific
```

The browser suite uses actual registered UI hosts and the rendered Reveal deck,
with its own temporary database, Chrome profile and cleanup. It checks UI
composition, validity, isolation, themes, X-ray, persistent nodes, terminal
mounting and static PDF behavior. It does **not** claim a scientific execution
or quota-enforced terminal-isolation pass. The full runtime goal additionally
requires end-to-end numerical, failure/recovery and deployment acceptance.

The separate `run-broker.sh scientific` gate drives the actual selectors,
preparation and Run controls in registered consumers through protected HTTP,
TLS NATS and private artifact storage. Its numerical children use the existing
explicit finite native test driver, not a production isolation bypass. It owns
and removes its test containers, browser and UI/worker processes. Retained
`scientific-results.json` contains inputs, results, provenance and measured times.
The scenario explicitly releases/reassigns a role, cancels cold preparation,
prepares a replacement and calculates again before testing worker loss. Warm
scientific calculations may finish between status polls; running-job cancellation
is tested separately with the finite-child TLS fixture. A failed scenario retains
`scientific-partial-results.json` for numerical diagnosis, never as a complete pass.
The successful gate automatically runs the direct-engine reference comparison
in its separate numerical project. To repeat that comparison on retained data:

```sh
julia --startup-file=no --compiled-modules=existing --project=playground/worker/profiles/line-parameters playground/runtime/test/study_result_validation.jl /tmp/EXACT_TEST_DIRECTORY/scientific-results.json
```

These commands define the acceptance checks; consult `RUNTIME_PLATFORM_PROGRESS.md`
for their actual pass/fail evidence. A configured `assigned_execution` capability
does not mean a worker is assigned, ready or certified for restricted execution.
