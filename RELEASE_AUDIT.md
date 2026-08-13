# v0.2.0 release audit

Audit baseline: `main` at
`88f5de54202f5a97c7556b4ceac662b997f13de2`, before release-preparation
edits. The baseline test evidence was 1,750 package assertions and 19 FEM
integration assertions.

## P0 — release blockers

- **Unregistered hard dependency:** GetDP.jl was referenced through
  `[sources]`, which is not acceptable for General registration. Resolution:
  remove the dependency and vendor the compatible frontend privately from
  commit `b1d91b0d8974ea642b772462edcf6e26299fdf0a`, tree
  `2a2fe2782b1c17a235a63d3639b99ae012dab1f5`.
- **Optional systems loaded by core:** Makie, Gmsh, and GetDP were hard
  dependencies. Resolution: package extensions and explicit activation.
- **FEM coupled to routine coverage:** the external solver integration was
  discoverable by the normal test runner. Resolution: a manual release-only,
  coverage-disabled gate with a pinned GetDP archive.
- **Publication ambiguity:** two workflows could create releases. Resolution:
  TagBot is the sole tag and GitHub-release authority.
- **Registry metadata:** package version and installation guidance targeted a
  development checkout. Resolution: version 0.2.0 and clean-install checks.
- **Credential scan:** no private credentials were found. The only tracked
  token-shaped literal is the existing public Codecov badge identifier, not an
  upload credential. The repository history was not rewritten.
- **Third-party licensing:** the vendored Julia frontend retains its upstream
  BSD license and provenance. `THIRD_PARTY_NOTICE.md` also identifies the
  registered Gmsh wrapper and the licenses of the separately supplied Gmsh and
  GetDP programs. No GetDP binary or downloader is distributed.

## P1 — correctness and maintenance

- Aqua was skipped and reported three method ambiguities. The dispatch
  intersections are now explicit, and Aqua is unconditional.
- Three plotting exports had no definitions. They were corrected while moving
  implementation code to extensions.
- Binder, the TODO-to-issue scraper, the duplicate release template, and the
  accidental standalone cable-builder tree were maintained dead weight and
  were removed.
- The documentation build used `warnonly`, rewrote source files during an
  opaque build phase, and opened a browser locally. It was split into explicit
  instantiate, doctest, build, and deploy phases.
- Development instructions tied mathematical documentation to the `calc_`
  prefix and required unsustainable cross-reference lists. The maintained
  convention now follows physical meaning and contains no “See also” sections.
- GetDP executable resolution could invoke a downloader. The compatibility
  frontend now checks `GETDP_EXECUTABLE`, then `PATH`, and otherwise errors.

## P2 — deliberately deferred

- Wholesale public renaming, including the historical `calc_` prefix.
- Broad module or type redesign outside the dependency-extension boundary.
- Speculative performance changes and cleanup of nontrivial historical code.
- Simplification of the future FEM API beyond the compatibility facade.

## Branch evidence

The abandoned `updates-runner` tips recorded before deletion were:

- local: `47b53243b8353985665fb08e56389814d0da2410`
- origin: `56bc4ee0b3d987fc56c0b6e752180b6fa7158c0d`

Local and GitHub branch listings after cleanup contain only `main`,
`refactor/cablebuilder`, and `gh-pages`.

## Verification evidence

The following locally executable gates passed on Julia 1.12.6:

| Gate | Result |
| --- | --- |
| Routine package tests, including unconditional Aqua | 1,734/1,734 assertions; no plotting/FEM/formatter dependency resolved |
| CairoMakie extension tests | 1,744/1,744 assertions |
| Release-only FEM integration | 19/19 assertions with GetDP 3.5.0 |
| GetDP archive verification | SHA-256 `d3c28fa18f20d6147b4c7367d4dd802e9f7ddb58c608688bbb71919dbca8041d` |
| Independent Aqua run | all checks passed, including ambiguities, exports, stale dependencies, compat, piracy, and persistent tasks |
| Documentation doctests | passed |
| Strict documentation build | passed without `warnonly` |
| SciML formatting | JuliaFormatter 2.12.4 returned `true` with `overwrite=false` |
| Commit policy | scoped, lowercase, 72-character gitlint policy validated with gitlint 0.19.1 |
| Optional dependency boundary | core loaded no Makie backend, Gmsh, or GetDP; missing-backend errors were actionable |
| Registry dependency audit | all 32 dependencies and weak dependencies have registered or standard-library UUIDs |

The baseline and release trees were also compared structurally: every baseline
test expression remains represented, the 19 FEM assertions moved intact to the
manual integration harness, and two plotting lifecycle assertions were added.
The numerical totals of the routine runner changed because FEM discovery was
removed and Aqua was made unconditional, not because covered behavior was
discarded.

The prerelease-Julia job is defined in CI and must run on GitHub because only
Julia 1.12 is installed locally. The immutable release commit and its clean
commit-based installation result are recorded in the release handoff; a commit
cannot embed its own SHA.
