# Planned, changed, delivered

This page records how LineCableModels 0.1.0 came to be, with what was planned on 2026-10-23,
what changed on the way and what the release delivers. Each date is the date of the
decision.

## The starting point

Before the release there were three branches. The `main` branch could no longer be
installed: the package it depended on for GetDP was archived. The 0.2.0 branch had grown
into a large tree of formulations, studies, reports and plots, much of it written by
automated agents, and its architecture had eroded. The finite-element
investigation branch had sound physics: the unified earth-return formula, the coaxial
computation and a finite-element backend validated against it.

On 2026-10-23 the release started from a green field. The first registered release would be
0.1.0, built on an orphan branch from an empty tree. Its scope was a data model and an
engine only, with a clean architecture. The designs were two bare wires, two insulated
wires, a coaxial cable with a core of circular strands, and the same cable with rectangular
strands. The finite-element backend came in because it was complete and sound. Plots,
observables, studies, exports and libraries stayed out. The previous release's code had
stated the intent as a prototype and had not held as production code. Its concepts remained
design intent, and its code was read but never merged.

Before any physics was written, two things were fixed. The first was the guards, so that no
cleanup would ever be needed. The second was the register of names, against the growth of
names that had marked the earlier releases. Every guard blocks at zero from its
first commit, with neither a baseline nor a count.

## What was planned

The release was planned as pieces, then as work packages with milestones.

1. The bootstrap, the guards, the register, the writing checks and the frozen references of
   the four cable systems, computed on the investigation branch at a known commit.
2. The materials, the earth and the data model, with the cross-sections, the declaration
   nodes, the design, the system, the blueprint and the builders of the surface syntax.
3. The engine, specified whole before its first line. It was implemented in the subset the
   four references need, and ported stage by stage against frozen intermediate matrices, so
   that each step had a number to hit.
4. The finite-element backend as a package extension, ported from the investigation branch
   with its geometry, meshes and solver processes, against the same references.
5. Uncertainty propagation with Measurements, added to the plan on 2026-10-25 as a condition
   of the release.
6. The release, with this documentation, the registration workflows and the audits.

It ran with three roles, an architect who specified and audited, a coordinator who ran the
steps and an executor under strict bounds. The [workflow](workflow.md) page describes them.

## What changed on the way

Each entry gives the change, the date and the reason. Together they are the architecture's
own history, and the [rules](rules.md) page lists what never returns.

### The data model

- **The declaration was redesigned before its first step** (2026-10-24). The previous
  release nested groups inside groups and conflated the electrical mapping with the tree.
  The redesign assumed that the whole concept, as it stood, was wrong. The result was one
  job per node. `Terminal` defines the electrical mapping. `Repeat` and `Lay` replace the
  group. `Enclosure` is a container with a fill and contents. The description is flat
  regions with tables.
- **Terminal names are identifications** (2026-10-24). A proposal typed the terminals and
  reserved names for grounding. The name of a terminal identifies it and does nothing else.
  The circuit role belongs to the connections, with an index per terminal name, a repeated
  index for a bundle and zero for elimination, every name present.
- **Connections are named tuples, and dictionaries were refused** (2026-10-24). A dictionary
  drops a duplicate key in silence and invites a default. Copies that repeat a name all map
  to one index, which is bonding by a repeated index.
- **Touching cables keep their declared placements** (2026-10-24). The previous release
  moved cables a micrometer apart because its mesh needed it, and its earth formula had
  inherited the guard, so a constraint of the mesh had passed for one of the physics. The
  system rejects overlap only. The finite-element geometry learned to mesh the pocket
  between touching cables.
- **`Shell` became `Offset`** (2026-10-24). A layer declared by its thickness is an operator
  on the outline inside it. It is not a physical element, and an operator does not take the
  name of one. `CircularPattern` takes a radius or an `Offset`.
- **`stranded` takes an explicit `boundary`** (2026-10-25), the finished outline of the
  formation, as the previous release had it.
- **`BentStrip` became `AnnularSector`** (2026-10-27). The name describes the form, as
  `Circle` and `Annulus` do. `DifferenceShape` stayed: it is the cross-section of every
  fill.
- **A fill's holes are the outlines of its outermost contents** (2026-10-27). Nested
  contents no longer produce nested holes.
- **A design's number type promotes over its materials** (2026-10-28), so that an uncertain
  material widens the type. The resolved shapes of a design are kept in that type.
- **The equivalent tube's bounds are taken over the free regions and the fills' containers**
  (2026-10-28), never over a fill's contents. A bound chosen by a last bit of rounding at a
  tie has no derivative.

### The engine

- **One formula object per family** (2026-10-25). The architect first proposed author-named
  types such as `Ametani()`. The existing primitive stayed: `Family.Formula{ID}` with the
  identifier in the type and a `Val` lift in the surface syntax. One author proposes
  several formulas, and the family namespace is the disambiguation.
- **One earth family, two slots** (2026-10-25). A wrapper that paired an impedance formula
  with a potential formula was withdrawn, because it made a property of one formula into an
  architectural concept. A formula gives one quantity or both by the methods it defines.
  If both slots take a formula of one identifier, the two formulas must be equal.
- **The integrals inside a formula** (2026-10-25). The earth formulas went through two
  overdesigned drafts, a `Quadrature` family and then a set of required functions read off
  the unified implementation. The structure became `SpectralIntegral`, `integrate` selected
  by a tag, and a formula that defines its own kernels.
- **No kernel types** (2026-10-26). A prototype built callable structs for the spectral
  terms. They mixed the kernels with `SpectralIntegral` and gave the engine something
  particular to one formula. Each integrand became one anonymous function written where the
  term is.
- **The levels of tools** (2026-10-26). A name at the engine or family level is declared by
  every formula of that level. The ray of a spectral integral became
  the engine's `ray!`. The folds and the contour function went.
- **`root` is Julia's `sqrt`** (2026-10-26), after 8.6 million calls showed identical
  results.
- **One Γ for impedance and admittance** (2026-10-26). Γ is the common propagation constant
  of the primary source that drives both quantities. It became an option of the formulation,
  read by the formulas that model the longitudinal dependence.
- **Problem, formulation and backend, and their names** (2026-10-26). A reading had bound
  each formulation type to one backend. That binding inflated the names and the structures,
  so it was refused. The formulation is the methodology and the backend the implementation:
  `CoaxialFormulation` and `LCMBackend`.
- **What goes in a buffer** (2026-10-26). A buffer stores a value computed once and reused,
  or one step's hand-off to the next. Inline mathematics is preferred wherever possible, so
  that the code does not conceal the physics it writes.
- **The result is built after the loop** (2026-10-26), from arrays the loop fills, so that
  the result's constructor validates the computed values once.
- **The allocation criterion is the loop's own bytes** (2026-10-26): the pivot vectors of
  the four factorizations per frequency, measured exactly.
- **The engine's interfaces before the finite-element port** (2026-10-26). The earth methods
  took the per-frequency arguments, so that a formula without shared state computes in its
  quantity method. The reduction plan became one function, `reductions`, that both backends
  call. The trace was typed per backend. The slots took their families' quantity names.
- **`Float32` promotes to `Float64` at the problem** (2026-10-25). A `Float32` input is
  accepted without a chase for its precision. `BigFloat` is unsupported.

### The finite-element backend

- **The drawing converged in a package module** (2026-10-27). The Gmsh extension and a
  future Makie extension draw with one minimal set of methods for the primitives, and each
  extension specializes them for its own target. The `CAD` module defines the curves, the
  contours and `draw`.
- **A native name or none** (2026-10-27). A proposed `launch` verb was refused, because the
  name of an act keeps one meaning in every context where it is used. Processes start
  with Base's `run`, inline. The pool became the one method of `mesh`.
- **The model data writer is a method of Base's `write`**, and the scan of the frequencies
  is `scan` with Base's `read` for one frequency (2026-10-27 and 2026-10-28). A name that
  bundled two acts under one meaning was refused.
- **A cable tangent to the ground surface is refused** (2026-10-26), as the unified formula
  refuses it, because no source meshes it. In-run mesh sharing was dropped as inert. MUMPS
  is the one solver. Failures after the run started raise `error` with the run directory and
  the log tail.
- **The coaxial criterion** (2026-10-28). The port's strand coordinates differ from the
  source's by one or two units in the last place, which moves the meshes at round-off. The
  port's code was not changed to reproduce the last bits of the source. The wire cases match
  the frozen rows within 1e-6. The coaxial cases are held at the 1 % accuracy level against
  the frozen rows, with the achieved maxima recorded as known status.
- **One process pool as a method of Base's `run`, and one walk as methods of Base's `map`**
  (2026-10-28). The backend computes every problem in deterministic `Float64`, the only
  numbers that Gmsh and GetDP accept.
- **The GetDP binary** (2026-10-28). The identifiers went through `:installed`, `:auto` and
  `:host` before they came to name the action that delivers the binary:
  `GetdpBinary{:discover}` and `GetdpBinary{:download}`. The Windows
  binary was repacked unchanged into a reproducible tarball.
- **The finite-element job runs on pull requests that are not drafts**, the whole suite with
  the reference cases (2026-10-28).

### Uncertainty propagation

- **The extension is methods only** (2026-10-28): `nominal` and four Bessel methods on
  SpecialFunctions' own names, declared in the extension's imports, with no wrapper name.
- **The finite-element backend drops the uncertainty** (2026-10-28) and returns a `Float64`
  result of the nominal problem. A Monte Carlo evaluation is a later method of `compute`.
- **The criteria were measured** (2026-10-28). The nominal parts agree within 1e-14 with the
  deterministic engine. The sensitivities agree within 1e-5 with central differences, at the
  step where truncation and noise balance.

### The tests

- **The test tree was audited for gamed tests** (2026-10-28): tests that reach no feature
  code and assertion counts that do not measure a feature. About 150 assertions that checked
  no feature were removed or folded.
- **Coverage at 95 %, from tests of feature results** (2026-10-28). The number is a project
  status on pull requests over the merged report of the four jobs. The audit called three
  branches dead. They were kept as live code. One is reached through an imaginary
  Γ and is tested by the continuity of the two forms at their switch. The other two lie
  beyond every admitted input and are recorded as uncovered.
- **The references became six** (2026-10-29). The overhead and the mixed bare-wire pairs,
  frozen from the investigation branch, put the air-side kernels through `compute`.

## What was delivered

The release is 49 commits between 2026-10-23 and 2026-10-29.

| Part | Content |
|---|---|
| the data model | the cross-sections, `Offset`, the declaration nodes, `CableDesign`, `CableSystem`, `CableBlueprint`, the builders |
| the engine | the eight families with one formula each, `CoaxialFormulation`, `FEM`, the two backends, the numerical input, the buffers, the loop, the reductions, the spectral tools, the result |
| the Gmsh extension | the model, the geometry, the meshes, the GetDP scan, the result, the run directory, the discovered or downloaded binary |
| the Measurements extension | `nominal` and the four Bessel methods |
| the references | six cable systems with their analytical rows, finite-element rows, intermediate matrices and geometry files, each from a cited commit of the investigation branch |
| the guards | twenty-six guards at zero, with planted negative controls |
| the register | every name of the package with its meaning, 350 lines |

The source is about 6 000 lines under `src/` and 3 000 under `ext/`, with 6 500 lines of
tests. The suite runs 3 389 tests without the finite-element cases and 4 615 with them. The
measured coverage with the finite-element job is 96.7 %.

The finite-element results against the analytical ones, as the references record them:

| Case | `e_Z` | `e_Y` |
|---|---|---|
| two bare wires | 0.57 % | 0.80 % |
| two insulated wires | 0.41 % | 0.045 % |
| coaxial, circular strands | 4.9 % | 10.5 % |
| coaxial, rectangular strands | 5.2 % | 10.4 % |
| overhead bare wires | 0.61 % | 0.82 % |
| mixed bare wires | 0.61 % | 0.62 % |

The coaxial differences are explained and recorded in `test/references/README.md`: the
analytical model treats the armor wires as a tube, and a two-dimensional field model has no
lay correction.

## Deferred, and work in progress

These parts of the previous release are not in 0.1.0. Each will be refactored into this
structure as its own change and never merged from the old code.

- The corpus of earth-return formulas of Carson, Pollaczek, Wise, Saad, Wedepohl and Wilcox,
  Gary, Lucca, Ametani 2009 and Xue. Each is one `EarthReturn.Formula{ID}`.
- The lossy dielectric law, pipe-type cables and ducts, sector and Milliken conductors,
  tapes, rectangular armor wires and compaction.
- Layered, frequency-dependent and equivalent homogeneous earth.
- Cable constants, modal decomposition and finite-length propagation.
- Grids and parametric studies, Monte Carlo evaluation, observables and plotting through
  Makie, reports and tables.
- Exports to ATPDraw, PSCAD and TRALIN, material and cable libraries, saving and loading.
- In the finite-element backend, the resume of a run, mesh reuse across runs, the ONELAB
  export, field maps and the GMRES solver.
