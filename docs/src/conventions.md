# Development conventions

LineCableModels follows the SciML formatter style. Run

```julia
using JuliaFormatter
format(".")
```

before committing.

## Ownership-centred recursive module layout

Source code is organised by **owner first** and **responsibility second**. The
directory names the module, domain concept, or action that defines a decision.
The filename names the operation implemented for that owner.
For example, Engine owns line-parameter result semantics, so their array and
accessor methods live under `engine/lineparameters/`. PlotBuilder owns the
renderer-independent plotting grammar, while the optional Makie extension owns
Makie figures and controls under `ext/`.

Apply these rules when adding or moving code:

1. Put a declaration or method with the owner whose API or numerical meaning
   determines when it changes.
   Do not group unrelated owners into `helpers`, `utils`, `handlers`, or other
   mechanism-first buckets.
2. Keep a function at the package root only when at least two independent
   sibling modules require the same meaning and neither sibling defines it.
   Calculation functions shared by Engine,
   ParametricBuilder, UQ, and external backends lives in `Grammar`. Engine-only
   names such as `Z`, `frequencies`, and `description` remain in `Engine`.
3. Name files after a precise responsibility: `interfaces.jl` declares owned
   hooks, `types.jl` holds owner-wide types, `base.jl` contains only Base/Core
   protocols, and an action file such as `compute.jl`, `validate.jl`, or
   `render.jl` keeps one visible call sequence.
4. Promote one file to a directory when the concept has several coherent
   responsibilities. A directory does not imply a Julia submodule. Introduce a
   submodule only for a separate namespace, dependency set, or independently
   stated interface.
5. Keep each `Module.jl` as an index: module documentation, explicit imports,
   exports, includes in dependency order, and child reexports. Constructors,
   algorithms, validation, and presentation methods belong in
   the indexed files.
6. Declare shared interfaces at the lowest legitimate owner. Children import
   and extend those functions explicitly. They must not mutate parent modules
   with `eval` or registration machinery. Sibling communication follows a
   directed import graph.
7. Keep optional-package translations in Julia package extensions. Core
   `src/` may define package-neutral requests and encoded values, but it must
   not import the optional package. Makie rendering and XLSX workbook writing
   therefore live under `ext/`.
8. Keep development-only environments and manual diagnostics under `dev/`,
   not `src/` or the automated test tree.
9. Assign file formats by purpose. ReportBuilder owns human-facing tables and
   workbook descriptions; ImportExport owns lossless persistence and
   external-tool interchange. The line-parameter workbook model therefore
   lives in `reportbuilder/xlsx.jl`, its XLSX.jl writer lives in the package
   extension, and ATP, PSCAD, and TRALIN translations remain in
   `importexport/`.

The following sanitised snapshot is the reference shape. The snapshot omits leaf files
that do not clarify the layout, but every shown path exists in the maintained
tree:

```text
src/
├── LineCableModels.jl                 # package index
├── interfaces.jl                      # genuine multi-owner root contracts
├── grammar/
│   ├── Grammar.jl                     # shared calculation-type index
│   ├── types.jl
│   ├── interfaces.jl
│   └── uncertainty.jl
├── units/
│   ├── Units.jl               # unit/quantity owner
│   ├── quantities.jl
│   ├── definitions.jl
│   └── scaling.jl
├── validation/
│   ├── Validation.jl                  # validation index
│   ├── interfaces.jl
│   ├── rules.jl
│   └── validate.jl                    # fixed validation action
├── datamodel/
│   ├── DataModel.jl
│   ├── packing.jl                     # cable packing constraints
│   ├── cabledesign/
│   │   ├── cabledesign.jl
│   │   ├── cableconstants.jl
│   │   └── dataframe.jl               # authored-design tables only
│   └── preview/                       # renderer-independent preview definitions
├── engine/
│   ├── Engine.jl
│   ├── interfaces.jl
│   ├── compute.jl                     # analytical choreography
│   ├── earthreturn.jl                 # named physical algorithm
│   └── lineparameters/
│       ├── lineparameters.jl
│       ├── quantities.jl              # Engine ↔ Units methods
│       └── plotdefinition.jl          # renderer-independent plot definition
├── plotbuilder/
│   ├── PlotBuilder.jl
│   ├── interfaces.jl
│   ├── composition.jl
│   ├── render.jl                      # fixed recipe action
│   └── backends.jl                    # extension activation only
├── uq/
│   ├── UQ.jl
│   ├── linearerror.jl
│   └── montecarlo/
│       ├── compute.jl
│       └── plot.jl
├── reportbuilder/
│   ├── ReportBuilder.jl               # report owner and index
│   ├── grammar.jl                     # fixed report action
│   ├── tables.jl                      # completed-result tables
│   ├── montecarlo.jl                  # UQ summary projection
│   └── xlsx.jl                        # human-facing workbook
└── importexport/
    ├── ImportExport.jl
    ├── interfaces.jl
    └── pscad/
        ├── pscad.jl
        ├── project.jl
        └── import.jl

ext/
├── LineCableModelsMakieExt/
│   ├── LineCableModelsMakieExt.jl     # Makie extension index
│   └── UIComponents.jl                # Makie-owned translation and controls
└── LineCableModelsXLSXExt.jl          # XLSX workbook writer

dev/
└── plotting/                          # manual plotting environment
```

Structural regression tests verify this layout. If a change needs a path
that contradicts the snapshot, first identify the new owner and responsibility
and update both the convention and its layout test. Do not
retain a misleading path through an alias or compatibility namespace.

## Public symbol selectors and `Val` dispatch

Public format, backend, preset, and ranking selectors accept ordinary
`Symbol`s. The generic's owner defines a local façade and retains `Val` as the
dispatch hook:

```julia
function owned_action(selector::Symbol, args...; kwargs...)
    return owned_action(Val(selector), args...; kwargs...)
end

owned_action(::Val{:example}, args...; kwargs...) = ...
```

The façade belongs beside the generic or the owned `Val` implementation. Julia
has no method that can rewrite calls to every unrelated generic, so do not add
a package-wide forwarding helper, macro, registry, or symbolic switch.

Apply this rule only to selectors that are part of a public call. Internal
dispatch axes such as `:self`, `:mutual`, plot dimensions, composition modes,
serialization tags, and `Val(OwnerType)` option dispatch remain `Val`-only.
Current public examples include:

```julia
Formulation(:analytical)
export_data(:atp, system, earth)
LineParameters(:tralin, path)
estimate[:match]
PlotBuilder.layout_preset(:single, 1)
```

Unknown selectors continue to reach the owner's unmatched `Val` dispatch.
Each owner decides whether that means an ordinary `MethodError` or an explicit
`ArgumentError`.

Versions follow [Semantic Versioning](https://semver.org/). Public behaviour is
kept compatible within a minor release. Deprecations must include a migration
path before removal.

Commit subjects use scoped Conventional Commits, start with a lowercase
description, and remain within 72 characters. For example:

```text
fix(engine): reject unsupported formulation options
```

Changes must include tests at the closest relevant scope. Core tests must not load
optional packages. CairoMakie is verified in a separate check.
Examples in docstrings should be executable and self-contained. Examples that
require fixtures, user interfaces, or external executables belong in developer
documentation instead.

Physical quantities state their SI units. A docstring for an implemented
equation places that equation in the method description and defines its
symbols. Use `# Notes` only for assumptions, limitations, or supplementary
mathematics.

Docstrings use `DocStringExtensions` abbreviations as the codebase-wide default
for generated signatures, type declarations, fields, and module inventories.
The complete templates, unit notation, mathematical requirements, and example
rules are documented in [Docstrings](@ref).
