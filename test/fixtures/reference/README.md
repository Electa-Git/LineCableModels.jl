# References

## Approved references

The first approved independent references sit beside their test, as the next section
asks. The directory `test/integration/fem_references/` contains one plain-text file per
case of the FEM investigation: `two_bare_wires.csv`, `two_insulated_wires.csv` and
`solid_1000mm2_single.csv`. Each file gives the FEM and analytical Z and Y at the 11 FEM
frequencies of its case. The header of each file states the source run folder and the
checkout that produced it. The user approved these files on 2026-10-10.

The FEM investigation regenerated them on 2026-10-12 with the production meshing of that
date. Their analytical values did not change, and their FEM values moved by at most
0.009 % of the largest diagonal entry.

The same directory holds cable cases copied on 2026-10-12. The header of each states what
FEM confirms:

- `cable_525kv_land_no_armour_dc_bipole.csv`: Y at every frequency, and Z up to 1.58 MHz.
  At 10 MHz the core self-impedance differs by 2.78 %.
- `cable_320kv_no_armour_dc_bipole.csv` and `cable_380kv_no_armour_ac_flat.csv`: Y at every
  frequency. Z differs from 4 Hz because the analytical stranded core includes the lay
  correction, which the 2D FEM model does not have.
- `cable_132kv_cigre_tb880_case0_630cu_trefoil.csv`: Y up to 39.8 kHz. Z differs from
  25 Hz because of the proximity effect between the closely spaced cables (spacing 2.2 times
  the outer radius), which the analytical model does not represent.

`test/integration/fem_references.jl` checks every frequency of each case against them.
The analytical result stays within 1 % of FEM by the metric of the user, for Z and for Y,
wherever FEM confirms the case.
That metric divides the largest entry difference by the largest diagonal entry. The
result also matches the recorded analytical values to a relative tolerance of 1e-12.

## Retired inherited references

The inherited numerical payloads have all been removed. No file in this directory is
an approved initial-release reference. Necessary scenarios and independent expectation
recipes belong beside their maintained test families.

Do not regenerate rejected expectations after a failure.
