# References

## Approved references

The first approved independent references sit beside their test, as the next section
asks. The directory `test/integration/fem_references/` contains one plain-text file per
case of the FEM investigation: `two_bare_wires.csv`, `two_insulated_wires.csv` and
`solid_1000mm2_single.csv`. Each file gives the FEM and analytical Z and Y at the 11 FEM
frequencies of its case. The header of each file states the source run folder. The user
approved these files on 2026-10-10.

`test/integration/fem_references.jl` checks every frequency of each case against them.
The analytical result stays within 1 % of FEM by the metric of the user, for Z and for Y.
That metric divides the largest entry difference by the largest diagonal entry. The
result also matches the recorded analytical values to a relative tolerance of 1e-12.

## Retired inherited references

The inherited numerical payloads have all been removed. No file in this directory is
an approved initial-release reference. Necessary scenarios and independent expectation
recipes belong beside their maintained test families.

Do not regenerate rejected expectations after a failure.
