# Reviewed numerical references

This gate is prepared but not enabled in CI: no reference has been approved.
It does not run Gauntlet, PSCAD, GetDP, a case importer or a sampler. It checks
the current coaxial engine against explicitly reviewed, pinned Gauntlet outputs.
The existing independent physics fixtures remain in `test/fixtures/reference`.

New deterministic campaign artifacts carry the serialized scalar problem actually
computed, including the normalized frequency samples, and the full
physical formulation declaration. Replay uses those stored inputs through the
package's existing deserializer and `Formulation` constructor, not today's case
catalogue. Frequency samples, matrix ordering, basis and reduction options are
not silently changed. A recorded `:default` exercises the current contextual
default; a changed result requires investigation and explicit review, not an
automatic reference refresh. Explicit author choices retain their own tests.

After reviewing an artifact's input data, backend comparisons, physical
assumptions and numerical results:

1. Pin its published archive in this directory's standard Julia `Artifacts.toml`.
2. Add a `[[references]]` entry to `approved.toml` with `id`, `artifact`, `file`,
   `sha256`, `review`, `Z_atol`, `Y_atol` and `rtol`. `file` is relative to the
   artifact root; `sha256` selects the exact reviewed JLD2 bytes. `review` records
   the review decision or its link. This manifest is the separate approval;
   the immutable source artifact can remain marked `:unreviewed`.
3. Choose tolerances explicitly. The existing `compare` API calculates RMS
   error across the full frequency vector for every Z/Y matrix entry. Every
   entry must satisfy its absolute **or** relative tolerance. `Z_atol` is in
   Ω/m, `Y_atol` in S/m, and `rtol` is a fraction, not percent. No tolerance is
   inferred from the candidate result or from another backend's RMS difference.
   Replay passes `atol=0` to `compare`: the scientific numerical-zero policy
   must not hide drift from a reviewed CI reference. Only these explicit
   manifest tolerances govern acceptance.
   The command reports the absolute and relative errors for every matrix entry,
   identifying a failing row and column rather than only a matrix-wide maximum.
4. Run `julia --project=test/gauntlet --startup-file=no test/numerical/runtests.jl`.
   Enable this command in CI only when the bindings and review are complete.

An empty approval manifest fails explicitly when this command is invoked. The
gate never generates, approves, binds or overwrites reference data. File hashes
protect the reviewed bytes; implementation fingerprints are provenance, not a
reason to invalidate the baseline when the implementation under test changes.
Historical artifacts without replay inputs still work in reports, but are not
silently reconstructed from current case files.

This initial gate supports deterministic, per-metre phase matrices from owned
coaxial calculations. External results inform their review; they are not
assumed numerically interchangeable. UQ moments and other problems require an
explicit reference-replay method before inclusion.

`@inferred compute` checks the replayed scalar computation. Wall-clock timing
is deliberately absent: the core suite owns hot-path allocation checks, while
controlled performance comparisons remain a separately recorded manual check.
