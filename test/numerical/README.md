# Numerical snapshots — deferred

Numerical snapshot testing will not be activated before the first stable
publication. The user will then choose a small number of Gauntlet artifacts.
There is no snapshot-approval task for this prerelease. See the governing
[testing policy](../../docs/src/developers.md#testing-policy).

The pasted snapshot concept describes a possible later continuity check; its
suggestion to begin before release is superseded. No reference is to be selected,
generated, refreshed or added to CI now. No new framework or dependency is needed.

## Existing inactive provision

`references.jl`, `runtests.jl`, `approved.toml` and `Artifacts.toml` are existing
inactive replay machinery. The list is empty and CI does not invoke its command.
The legacy command errors on that empty list; this is not a failed calculation
or a missing prerequisite for the first stable release. Its protocol tests use
temporary records and do not select any numerical baseline.

The current reader supports deterministic coaxial per-metre phase matrices and
uses per-entry RMS across frequencies. It is not a finished general snapshot
system: do not advertise FEM/UQ replay or elementwise continuity that it does not
implement. Leave it dormant rather than expanding it speculatively.

## Minimal later design

Once the user selects artifacts after stable publication:

- Reuse existing Gauntlet records, JLD2, the runner and `Test`. Adapt only what the
  selected records require.
- Compare a backend with its own retained output across revisions. Preserve input
  and formulation settings, frequency/terminal identity, units/basis, relevant
  numerical quantities and execution provenance.
- Check real/imaginary components at individual entries/frequencies with explicit
  continuity tolerances; RMS may summarize differences. Report the differing
  entry, old/new values and permitted difference.
- Keep the selected reference fixed. Updating it requires an explicit user
  decision; CI must not generate missing references or refresh them after passing.
- Keep private workspaces, struct layouts, temporary paths and incidental output
  out of the snapshot. A source revision is provenance, not a reason to bypass
  the comparison.

These are behavior snapshots, not scientific approvals. No previous-push service,
reference promotion workflow, new artifact catalogue or scientific manager is
part of the present work.
