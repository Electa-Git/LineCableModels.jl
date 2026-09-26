# Three bare wires: retained baseline

The physical fixtures live in `test/support/scenarios.jl` and select no formula.
Run `dev/run_three_bare_wires_baseline.jl` in an environment with this package
and Gmsh available; the script does not activate or install an environment.

Each additive `capture-*` directory retains five complete nine-frequency scans
for current unified and FEM (`quasi_tem`), or the actual failure for each scan.
CSV files contain every ordered entry of phase-domain Z [Ω/m] and Y [S/m] with
round-trip Float64 precision. The manifest records inputs, source provenance,
timings, failures and checksums. FEM artifacts remain at the recorded persistent
run directories under `.linecablemodels/fem`; preserve those directories with
the capture. Nothing here silently approves a reference or sets tolerances.

The baseline was reviewed before production cleanup. The human authorized
non-degradation of the observed agreement, allowing improvement and insignificant
Y changes; this does not assert exact FEM agreement. Other formulas can reuse the
physical fixtures without inheriting the unified comparison requirement.

`capture-20260917T102438-Nt9nic` is the immutable pre-refactor baseline.
`validation-20260917T111424-R4HWRE` retains all post-refactor entries, before/after
errors and fresh-workspace timings. Its comparison never reruns FEM or replaces
the captured reference.

The versioned checkpoint contains the original numerical payloads and capture
metadata, plus the Engine-convergence starting scan
`validation-20260919T144543-aqKVmx` and final scan
`validation-20260919T170632-eFIkJx`. Historical source copies, dirty patches,
execution logs, native FEM artifacts and other intermediate validation captures
remain local evidence; they are not required to run the maintained fixture tests.
Historical serialized problems retain their original schema and are not current
execution inputs. Tests construct the physical fixture from `scenarios.jl`.
