# Legacy Pluto showcase: recovered presentation intent

This document preserves the presentation intent of the retired Pluto showcase without
preserving its implementation. It is a rebuilding brief for new showcase pages, not a
line-by-line transcription of a notebook.

## Recovery scope

The history was traced backward from `main` through `showcase/showcase1.jl`. Two
materially different presentations were found:

| Version | Historical reference | Visible presentation identity | Distinguishing purpose |
| --- | --- | --- | --- |
| A | `6bc9af9b`, committed 15 September 2025 | *Improved models and tools for cable systems*, dated 10 April 2025 | An end-to-end tour of the LineCableModels toolbox and its data-exchange workflow |
| B | `e2edce7d`, committed 16 September 2025 | *Uncertainty of Frequency Dependent Impedance Parameters for Transmission Assets*, Second Annual Belgian Energy Transition Workshop, dated 16 September 2025 | A focused engineering study connecting uncertain cable inputs to frequency-dependent sequence impedances |

The intervening revision at `196bb556` is a transition between these stories: it adds
the uncertainty-oriented introduction while still retaining most of the original
toolbox tour. Its useful material is represented in Version B below rather than treated
as a third coherent presentation.

The later revisions `98fc847f`, `d3a2ebf8`, and the `main` revision `2356399e` retain
Version B's narrative. Their changes are formatting or API maintenance, not a third
content version. The Markdown content is identical across those four revisions.

The historical file `showcase/ backup 1.jl` was also checked. It was committed as an
empty file and contains no recoverable presentation material.

Excluded deliberately:

- Pluto setup, Binder setup, CSS, custom navigation, and DOM manipulation
- Julia construction code and obsolete package APIs
- local plotting helpers and notebook-only output formatting
- images and generated artifacts

Every visual below is represented by a placeholder describing what the rebuilt page
should communicate.

## Shared opening narrative

Both versions begin from the same institutional and system-level context. This can be
rebuilt once and reused by either presentation.

### 1. Cover and attribution

Present the talk title, presenter, KU Leuven Etch/EnergyVille affiliation, and date.
Keep partner branding secondary to the title.

> [Image placeholder: Etch, EnergyVille, and KU Leuven branding]

### 2. Etch and the future grid

Introduce Etch as the Energy Transmission Competence Hub. Frame its mission around
future-proof electricity networks with large-scale HVDC and underground-cable
integration.

The four challenges named in the notebooks were:

- more underground cables;
- protection of cable-based systems;
- control interactions;
- resilient HVDC grids.

> [Image placeholder: future-grid/HVDC system graphic]

### 3. Why undergrounding matters

Balance the motivation and difficulty:

- underground and HVDC assets can improve environmental impact, reliability, and
  public acceptance;
- cost, technical complexity, accessories, installation environment, and system
  interaction prevent this from being a simple substitution for overhead lines;
- research must connect detailed asset models to planning, operation, diagnostics,
  protection, and control.

> [Image placeholder: two contrasting underground-cable installation photographs]

The research-roadmap wording changed between versions. Version A emphasizes modeling
and design, uncertainty/optimization, diagnostics, protection/control, and cables in
multi-GW grids. Version B sharpens this to cables and joints in complex environments,
multi-terminal HVDC systems, operational/diagnostic tools, and cable hosting capacity.

## Version A — toolbox and workflow tour

### Central intent

Demonstrate that LineCableModels can carry a real cable specification through one
continuous engineering workflow:

1. describe materials and physical layers;
2. assemble and inspect the cable design;
3. calculate base electrical quantities;
4. store and exchange the design;
5. place three cables in a physical system;
6. export that system for use in an EMT tool.

The desired audience takeaway is breadth and continuity: the same domain description
supports geometry, parameter calculation, reuse, system assembly, and downstream data
exchange.

### Suggested page sequence

#### A1. What the toolbox covers

Summarize the promised capability rather than listing package types:

- coaxial cables with arbitrary conductor and insulation layering;
- solid, tubular, and stranded conductors, screens, armoring, and semiconductive
  layers;
- DC and AC line/cable constants with temperature, stranding, and helical corrections;
- propagation and earth-return effects across frequency;
- analytical and finite-element formulations plus EMT-tool interfaces;
- reusable materials and cable-design libraries;
- uncertainty propagation as a cross-cutting capability.

> [Image placeholder: LineCableModels identity beside a layered cable cross-section]

#### A2. Reference cable specification

Use the single-core aluminum 1000 mm², 18/30 kV NA2XS(FL)2Y cable as the running case.
Show a layer schedule beside a product photograph. The layer schedule should move from
the stranded conductor outward through semiconductive layers, main insulation, wire
screen, copper and water-blocking tapes, aluminum moisture barrier, and PE jacket.

The historical nominal checkpoints were 0.0291 Ω/km DC resistance, 0.39 µF/km
capacitance, and 0.3 mH/km trifoil inductance. Treat them as comparison values, not as
hard-coded truths for a future data model.

> [Image placeholder: cable photograph beside a layer/thickness/diameter table]

#### A3. Materials and physical construction

Explain that the material library supplies reusable electromagnetic properties, then
build the design in visually meaningful stages:

1. a 61-wire aluminum stranded core and its semiconductive/XLPE insulation stack;
2. copper wire screens, equalizing tape, and water-blocking layer;
3. aluminum moisture barrier and outer PE protection.

After each stage, show the resulting cross-section. The point is not to expose old
constructors; it is to make the layer hierarchy and cumulative geometry legible.

> [Generated-figure placeholder: core and main-insulation cross-section]

> [Generated-figure placeholder: cross-section after the wire screen is added]

> [Generated-figure placeholder: finished cable cross-section with all layers]

#### A4. Base cable quantities and validation

Show three levels of result:

- a compact comparison with datasheet R, L, and C values;
- equivalent electromagnetic properties grouped by cable component;
- a detailed engineering report for inspection or debugging.

The visual story is “specification in, inspectable engineering quantities out.” Tables
or metric cards are more important here than source code.

> [Output placeholder: datasheet-versus-calculated RLC summary]

> [Output placeholder: detailed component parameter table]

#### A5. Reuse and data exchange

Demonstrate that completed cable designs and material definitions can be stored in
libraries, listed, and serialized for transfer between projects. The notebook used JSON
for readable exchange and noted a binary option for sensitive specifications.

This should be presented as a workflow and governance feature, not as a serialization
API tutorial.

> [Output placeholder: small library/catalog view with one stored cable design]

#### A6. Three-phase cable system

Place three copies of the cable in trifoil at approximately one metre burial depth and
associate them with a frequency-dependent earth model spanning 1 Hz to 1 MHz.

Explain phase mapping explicitly:

- the three cores map to phases A, B, and C;
- screens and jackets map to ground and are eliminated from the final phase model;
- shared phase assignments would represent bundled conductors.

> [Generated-figure placeholder: three single-core cables in buried trifoil formation]

#### A7. PSCAD handoff

Close the workflow by showing the cable-system preview and the resulting PSCAD-facing
artifact. The important message is that the domain model is useful beyond the Julia
session.

> [Image placeholder: PSCAD import/export result for the three-phase cable system]

### Version A interaction model

The original story was progressive rather than strongly parametric. A rebuild should
prefer stable, staged views of one design. Useful interaction would be layer inspection,
switching between construction stages, or selecting the result-table level. It does not
need to reproduce notebook cells or execute each construction step in front of the
audience.

### Version A closing message

- conductor and material modeling materially affects cable parameters and propagation;
- earth, soil, and modal formulations need continued expansion;
- more cable models and uncertainty studies are needed;
- general multilayer formulations are important for semiconductive materials;
- finite-element validation is part of the research direction;
- interoperability with engineering tools is a practical requirement.

## Version B — uncertainty-to-impedance application study

### Central intent

Explain how uncertainty in geometry, material properties, field conditions, and modeling
choices propagates into the frequency-dependent parameters used for power-system and
EMT studies.

The desired audience takeaway is causal: uncertain physical inputs change base cable
quantities, those quantities feed a three-phase frequency-domain model, and the final
sequence impedances carry visible uncertainty.

### Suggested page sequence

#### B1. Where parameter uncertainty comes from

Group the sources into internal and external origins:

- geometrical and material tolerances;
- field-dependent values such as resistivity and the as-built conductor layout;
- nearby interferences;
- choices made in the parameter and EMT modeling procedure.

Use skin/proximity behavior and earth return as the two physical examples. They explain
why frequency-dependent impedance cannot be reduced to a single nameplate value.

> [Image placeholder: conductor current-density/skin-effect illustration]

> [Image placeholder: earth-return field or current-path illustration]

#### B2. Uncertainty propagation primer

Introduce every physical quantity as a nominal value with an associated uncertainty.
Keep the mathematical message compact:

- absolute uncertainties combine quadratically for sums and differences;
- relative uncertainties combine quadratically for products and quotients;
- a general function uses first-order sensitivity through its partial derivatives.

Call out the counter-intuitive point retained in the notebook: subtracting nominal
values does not subtract their uncertainties.

> [Diagram placeholder: uncertain cable dimensions flowing through a calculation]

#### B3. Reference application

Reuse the same 1000 mm², 18/30 kV aluminum cable and its layer schedule. This anchors
the study in a recognizable design before any parameters are varied.

> [Image placeholder: reference cable photograph beside its layer schedule]

#### B4. Interactive base RLC quantities

Expose the physical controls that defined the historical experiment:

- number of conductor-wire layers;
- wire diameter and its percentage uncertainty;
- inner-semiconductor thickness and its percentage uncertainty;
- main-insulation thickness and its percentage uncertainty;
- outer-semiconductor thickness and its percentage uncertainty.

Update one cable cross-section and three headline outputs from those inputs:

- conductor cross-sectional area;
- core resistance per kilometre;
- main-insulation capacitance per kilometre.

The presentation should emphasize cause and effect. Keep domain construction behind a
small page-local boundary so the current data-model API can replace the old one without
changing the authored narrative.

> [Interactive-figure placeholder: cable cross-section responding to geometry inputs]

> [Metric placeholder: area, resistance, and capacitance with uncertainty]

#### B5. Equivalent EMT representation

Show the detailed cable next to its reduced/equivalent representation and provide a
compact base-parameter table. Explain which physical detail is collapsed and which
electrical behavior the equivalent is expected to retain.

> [Comparison placeholder: detailed cable and equivalent EMT cable cross-sections]

> [Output placeholder: equivalent base-parameter table]

#### B6. Frequency-domain system setup

Create a three-cable trifoil system at approximately one metre burial depth. Associate
it with a frequency-dependent earth model from 1 Hz to 1 MHz and show the system
geometry before presenting electrical results.

> [Generated-figure placeholder: trifoil cable system within the earth model]

#### B7. Sequence-component impedance result

Carry the uncertain cable inputs through the line-parameter calculation and a Fortescue
transformation. Present zero-, positive-, and negative-sequence behavior over logarithmic
frequency.

The historical result view contained:

- resistance and inductance curves versus frequency;
- uncertainty bars or bands, deliberately enlarged in the notebook to remain visible;
- a tabular R, L, C, and G representation for each sequence.

A rebuild should label any visual uncertainty amplification explicitly. It must not make
an enlarged display band look like the physical uncertainty itself.

> [Plot placeholder: zero-, positive-, and negative-sequence resistance versus frequency with uncertainty]

> [Plot placeholder: zero-, positive-, and negative-sequence inductance versus frequency with uncertainty]

> [Output placeholder: sequence-component RLCG tables]

### Version B interaction model

The interaction is the argument of the presentation, not decoration. Changing a nominal
dimension should update the geometry, base quantities, and final frequency response.
Changing only an uncertainty should leave the nominal curve centered appropriately while
changing its interval. Expensive frequency-domain work may be separated from continuous
geometry feedback, provided the page makes the cached/recomputed state explicit.

The old notebook duplicated its geometry controls for the base-quantity and sequence-
impedance sections. That is an implementation artifact, not a content requirement. A
new page can use one shared parameter state or two intentionally labeled scenarios.

### Version B closing message

- correct conductor and material models are necessary for credible propagation studies;
- earth-return and modal choices are part of the uncertainty budget;
- additional cable models and systematic validation remain necessary;
- multilayer semiconductor formulations and FEM comparison are priority extensions;
- uncertainty should be carried into the engineering outputs rather than discarded at
  the input boundary.

## Recommended reconstruction boundaries

The following content units fit the current discoverable-page showcase without reviving
the Pluto notebook structure:

| Page unit | Reuse | Primary content |
| --- | --- | --- |
| Institutional context | Shared | Etch mission, future-grid challenges, undergrounding motivation |
| Reference cable | Shared | Product identity, layer schedule, nominal checkpoints |
| Layered cable workflow | Version A | staged physical construction and previews |
| Engineering quantities | Version A | validation tables and detailed results |
| Libraries and exchange | Version A | reusable designs/materials and external handoff |
| Trifoil and PSCAD | Version A | system geometry, phase mapping, PSCAD artifact |
| Sources of uncertainty | Version B | tolerance, field, interference, and modeling origins |
| Uncertainty primer | Version B | propagation rules and interpretation warning |
| Parametric base quantities | Version B | geometry controls, cross-section, R and C feedback |
| EMT equivalent | Version B | detailed-versus-reduced comparison |
| Sequence impedance | Version B | trifoil, earth model, frequency curves, intervals, RLCG table |

Keep narrative prose, page layout, and image placeholders independent of package object
construction. Each page may own a small adapter that creates its current domain object
and extracts the results it displays. This preserves the recovered story while allowing
the LineCableModels tree and plotting APIs to continue evolving.

