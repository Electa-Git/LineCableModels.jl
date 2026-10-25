# Reference results

Per-unit-length line parameters of the four cable systems of this release, computed with the
Unified earth formulation and with FEM. The tests compare the package against them.

The case scripts in `test/manual/cases/` of the branch `release/v0.2.0-fem-investigation` of
this repository computed them, each at the commit that the table gives. Commit `2bff171d` adds
the two coaxial scripts to `a7640beb` and changes nothing else. Commit `989af96c` reduces every
terminal to its equivalent conductor by the current share of its regions, as this release does.
This reduction changes only the core of `coaxial_circular_strands`. At `989af96c`, the other
three runs repeat their analytical values within 4e-15 of the largest diagonal entry, and their
FEM values within 1e-8. Run each script from the repository root as
`julia --project=test test/manual/cases/<case>.jl`.

| Case | System | Commit | Run |
|---|---|---|---|
| `two_bare_wires` | two copper wires of radius 42.5 mm, 1 m apart, buried 1 m deep in earth of 0.1 Ω·m (εr 1) | `a7640beb` | `20261016-162340` |
| `two_insulated_wires` | the same wires with 1 mm of insulation (εr 2.3) | `a7640beb` | `20261016-162539` |
| `coaxial_circular_strands` | one coaxial cable with a core of circular strands, buried 1 m deep in earth of 100 Ω·m (εr 10) | `989af96c` | `20261025-000358` |
| `coaxial_rectangular_strands` | the same cable with a core of rectangular strands | `2bff171d` | `20261023-215619` |

The coaxial cable, from the center out:

- The core of `coaxial_circular_strands` is 37 copper wires of 3 mm: a straight center wire and
  courses of 6, 12 and 18 wires at lay ratio 12.
- The core of `coaxial_rectangular_strands` is a copper center wire of 3 mm and courses of 6, 12
  and 18 copper strands of 3.0 × 3.0 mm at lay ratio 12. Each strand is bent into an annular
  sector of the same area.
- Both cores lie under 1 mm of semicon, 8 mm of XLPE and another 1 mm of semicon.
- A copper sheath of 1 mm and a PE bedding of 3 mm follow.
- The armor is 40 steel wires of 4 mm at lay ratio 10, with relative permeability 36.6, under
  a PE jacket of 4 mm.

## Files

Each case folder holds `matrices.csv`, as its script wrote it. Its comment lines give the
conductors in row and column order, and the units. Each row is one entry of one matrix:
`source` (`analytical` or `fem`), `frequency_hz`, `quantity` (`Z` in Ω/m or `Y` in S/m), `row`,
`column`, `real`, `imaginary`. The analytical rows hold 101 frequencies from 0.1 Hz to 10 MHz,
evenly spaced on a logarithmic scale. The FEM rows hold 11 of them, one per 0.8 decade.

Running a script again gives analytical values identical to the bit and FEM values within about
1e-8 of the largest diagonal entry. FEM meshes depend on incidental state such as the length of
the run directory's path.

## Status

FEM against the analytical results, per FEM frequency, as e_Z = max|ΔZ_ij| / max|Z_ii| and
e_Y = max|ΔY_ij| / max|Y_ii|. The FEM accuracy level is 1 %.

- `two_bare_wires`: e_Z at most 0.57 % and e_Y at most 0.80 %, under 1 %.
- `two_insulated_wires`: e_Z at most 0.41 % and e_Y at most 0.045 %, under 1 %.
- `coaxial_circular_strands` and `coaxial_rectangular_strands`: e_Y from 8.5 % to 10.5 %, e_Z up
  to 5.2 %. Differences between the two models account for them. Variant runs of the
  rectangular-strand case on `a7640beb` remove one difference at a time and show each part. The
  last variant was run for both cases.
  - e_Y, 8.5 % to 10.5 % at every frequency, in Y33: the analytical model treats the armor
    wires as a tube. The wires meet the 3 mm bedding only at their crowns, with fill (εr 1) in
    the wedges between them. With the armor as a steel tube, e_Y is at most 0.094 %.
  - e_Z below about 25 Hz, up to 4.5 % in Z33: the armor lay correction, which the analytical
    equivalent conductor carries and a two-dimensional FEM model does not. With straight armor
    wires, e_Z at 0.1 Hz is 0.49 %.
  - e_Z above about 40 kHz: the core lay correction, in the same way. In Z11 at 10 MHz, e_Z
    reaches 5.2 % for the rectangular strands and 4.9 % for the circular strands. With a
    straight core, at lay ratio 1e6, e_Z at 10 MHz is 0.55 %.
  - e_Z from 25 Hz to 40 kHz, 1.1 % to 1.9 % with both lays straight: the magnetic armor wires
    against the tube.
  - With the armor as a tube and a straight core, e_Z is at most 0.42 % and e_Y at most
    0.094 % for the rectangular strands, and 0.63 % and 0.13 % for the circular strands.
