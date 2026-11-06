# Reference results

Per-unit-length line parameters of the six cable systems of this release, computed with the
earth-return formula `:martins2026` and with FEM. The tests compare the package against them.

The case scripts in `test/manual/cases/` computed them, each at the commit of this repository
that the table gives. Commit `2bff171d` adds
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
| `overhead_bare_wires` | two bare copper wires of radius 42.5 mm, 1 m apart, 1 m above earth of 0.1 Ω·m (εr 1) | `e41bba13` | `20261031-100133` |
| `mixed_bare_wires` | the same wires at (0, 1) m and (1, −1) m, one overhead and one buried | `e41bba13` | `20261031-100324` |

The `fem` rows of `coaxial_rectangular_strands` come from the run `20261106-131133`, at the
commit that meshes each of its rectangular strands with transfinite strips. That run computes
the case as `test/references/problems.jl` declares it, with the options of the FEM test of the
suite. The `analytical` rows are those of `2bff171d`.

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

Each case folder also holds `physical.geo`, the geometry file of the FEM run of its
script. At commit `989af96c`, the four original scripts ran unchanged with
gmsh_jll 4.15.2. Their FEM rows match `matrices.csv` within 5.4e-9 of the largest
diagonal entry.

Each `physical.geo` is the retained file of its case's run of 2026-10-25. It holds the lines
that the extension has written since then, without the lines `FEMConductorWidth` and
`FEMConductorSector`, which the extension no longer writes. The file of
`coaxial_rectangular_strands` is the retained file of its run `20261106-131133`. The geometry
tests compare against these files.

Running a script again gives analytical values identical to the bit and FEM values within about
1e-8 of the largest diagonal entry. FEM meshes depend on incidental state such as the length of
the run directory's path.

## Intermediate matrices

Each case folder also holds `stages.csv`, and the folder of `two_insulated_wires` holds
`stages-gamma.csv`. They give the four matrices from which the analytical computation builds `Z`
and `Y`, at the same 101 frequencies:

- `Zin`, per conductor: the series impedance inside each cable from its conductors to its outer
  surface, in Ω/m.
- `Pin`, per conductor: the potential coefficients inside each cable in m/F.
- `Zg`, per cable: the earth-return impedance of the cables from their outer surfaces in Ω/m.
- `Pg`, per cable: the earth-return potential coefficients of the cables from their outer
  surfaces in m/F.

Copies of the four original case scripts computed them at `989af96c`. Copies of the overhead
and mixed scripts computed them at `9621e504`. The runs appear below. Each copy keeps the
script's system and formulation, leaves out FEM, computes with `options = (trace = true,)`
and writes the four matrices of `result.trace`. The comment lines give the
conductors and the cables in row and column order, and the units. Each row is one entry of
one matrix: `source` (`analytical`), `frequency_hz`, `matrix`, `row`, `column`, `real`,
`imaginary`.

`Z` is `Zin` plus `Zg`, whose entry for two cables adds to every pair of their conductors, and `Y`
is `jω (Pin + Pg)⁻¹`, with `Pg` added in the same way. Rebuilt so, each `stages.csv` gives the
analytical rows of its `matrices.csv` within 4e-15 of the largest diagonal entry.

`stages.csv` uses `Γ = 0`, as `matrices.csv` does. `stages-gamma.csv` prescribes the longitudinal
wavenumber `Γ_k = j·1.01·2π·f_k/c` on both earth matrices, with `c = 299792458` m/s. It changes
`Zg` and `Pg` by up to 7.5e-4 of their values and leaves `Zin` and `Pin` unchanged. A second run of
each copy gives its data rows again, byte for byte.

| File | Run |
|---|---|
| `two_bare_wires/stages.csv` | `20261025-151919` |
| `two_insulated_wires/stages.csv` | `20261025-151919` |
| `two_insulated_wires/stages-gamma.csv` | `20261025-152032` |
| `coaxial_circular_strands/stages.csv` | `20261025-151920` |
| `coaxial_rectangular_strands/stages.csv` | `20261025-152029` |
| `overhead_bare_wires/stages.csv` | `20261028-235644` |
| `mixed_bare_wires/stages.csv` | `20261028-235756` |

## Status

FEM against the analytical results, per FEM frequency, as e_Z = max|ΔZ_ij| / max|Z_ii| and
e_Y = max|ΔY_ij| / max|Y_ii|. The FEM accuracy level is 1 %.

The port's run gives these maxima against the analytical rows:

| Case | e_Z | e_Y |
|---|---|---|
| `two_bare_wires` | 0.568 % | 0.798 % |
| `two_insulated_wires` | 0.407 % | 0.0447 % |
| `coaxial_circular_strands` | 4.92 % | 10.5 % |
| `coaxial_rectangular_strands` | 5.20 % | 10.4 % |
| `overhead_bare_wires` | 0.610 % | 0.817 % |
| `mixed_bare_wires` | 0.609 % | 0.617 % |

The file `two_bare_wires/xue2021.csv` gives `analytical` and `fem` rows of `Z` and `Y` at 1 kHz
for the wires of `two_bare_wires` on a layered earth. The upper layer is 100 Ω·m and 5 m thick
with εr 10. Below it lie 10 m of 500 Ω·m with εr 20 and a half-space of 50 Ω·m with εr 5. Both
rows use the equivalent earth `:xue2021`, of 62.05 Ω·m and relative permittivity −48 962 at that
frequency.

The FEM row deviates by e_Z 0.503 % and e_Y 0.735 %, and the pair by 0.568 % and 0.798 % on its
homogeneous earth. So the field solve computes with the negative permittivity of the equivalent
layer within the FEM accuracy level. Both computations ran once, at the commit that adds the
file. The test suite repeats the analytical one alone.

For `two_bare_wires`, `two_insulated_wires`, `overhead_bare_wires` and
`mixed_bare_wires`, the meshes coincide with the source meshes.
Their rows agree with the frozen `fem` rows within 1e-6 scaled.
The frozen wire-pair rows use a tolerance of 1e-6 of the largest diagonal entry on the reference
build, x86_64 Linux with Gmsh 4.15.2, and 1e-2, the FEM acceptance metric, on other builds.
For the coaxial cases, the maxima against the `fem` rows frozen on 2026-10-25 are:

| Case | Scaled difference in Z | Scaled difference in Y |
|---|---|---|
| `coaxial_circular_strands` | 5.84e-5 | 3.53e-3 |
| `coaxial_rectangular_strands` | 1.08e-5 | 1.36e-3 |

The coaxial meshes are two discretizations of identical definitions.
Many points of the core and armor wires lie on common circles.
Triangulation breaks ties on those circles by round-off.
The port coordinates and native size values differ from source values by one or two units in the last place.
The port composes the data model definitions.
The source evaluated each point at its construction site.

All scaled differences and the spreads `D_Z` and `D_Y` below are dimensionless.
Diagnostic A runs the port solve on the source meshes.
Diagnostic B uses the source numbers in the port inputs and keeps the port names.
It substitutes 367 lines of `physical.geo` for the circular case and 442 for the rectangular case.
It substitutes 11 lines of `model_data.pro` for each case.
The meshes are identical, with 57 816 nodes for the circular case and 45 353 for the rectangular case at the first frequency.
The maximum scaled differences are:

| Case | A, Z | A, Y | B, Z | B, Y |
|---|---|---|---|---|
| `coaxial_circular_strands` | 6.5e-11 | 6.5e-12 | 6.5e-11 | 6.1e-12 |
| `coaxial_rectangular_strands` | 4.4e-10 | 3.2e-10 | 4.4e-10 | 8.2e-11 |

Diagnostic C compares the default mesh with `mesh_size_factor = 0.5` at the first and last frequency.
Across the two cases, D_Z reaches 1.8e-4 at the first frequency and 2.4e-3 at the last.
D_Y is 3.1e-3 for the circular case and 1.1e-3 to 1.7e-3 for the rectangular case.
The difference between the port and source is 0.24 to 1.26 times D_Y.

The source gives the following results against the analytical rows:

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

## Cable systems of the case scripts

The scripts in `test/manual/cases/` declare five more cable systems. Each writes `<case>/matrices.csv` in the columns
above. Four of them hold `analytical` rows of `CoaxialFormulation(bundle = false, kron = false, transposition = false)`
at the 101 frequencies and `fem` rows of `FEM` with the same options at 11 of them. Each conductor of a cable is a
terminal of its own.

The fifth, `cable_30kv_630mm2_trefoil`, computes the cables of the example `examples/cable_30kv_630mm2_trefoil.jl` at
the 50 multiples of 50 Hz up to 2500 Hz. Both formulations remove its grounded screens by Kron reduction and transpose its phases.
Every system lies in earth of 100 Ω·m (εr 10). Each script ran once, with gmsh_jll 4.15.2, at the commit that adds its
file. The run is the folder of its results under `.linecablemodels/cases/<case>/`. Run each script from the repository
root as `julia --project=docs test/manual/cases/<case>.jl`.

| Case | System | Run |
|---|---|---|
| `cable_30kv_630mm2_trefoil` | three 18/30 kV cables with a solid aluminum core of 630 mm² and a copper screen 0.3 mm thick, in trefoil with the cables in contact, the lower two at 1 m depth | `20261104-004949` |
| `cable_30kv_na2xs2y_630mm2_trefoil` | three 18/30 kV cables with a stranded aluminum core of 630 mm² and a screen of 60 copper wires under a copper tape, in trefoil with the cables in contact, about 1 m deep | `20261104-005450` |
| `cable_132kv_630mm2_flathor` | three 132 kV cables with a stranded copper core of 630 mm², a screen of 19 copper wires under a copper tape and an aluminum tape under the jacket, 0.5 m apart in a row at 1 m depth | `20261104-010142` |
| `cable_18kv_1000mm2_trefoil` | three 18 kV cables with a stranded aluminum core of 1000 mm², a screen of 49 copper wires under a copper tape and an aluminum tape under the jacket, in trefoil with 1.1 diameters between centers, about 1 m deep | `20261104-011055` |
| `cable_380kv_2000mm2_flatver` | three 380 kV cables with a stranded copper core of 2000 mm², a screen of 14 copper wires under a copper tape and an aluminum tape under the jacket, 0.5 m apart in a column from 1 m to 2 m depth | `20261104-012742` |

Each stranded core has a central wire and courses of 6, 12, 18 and 24 wires. The courses of 12 and 24 wires start half
a pitch from the others.

FEM against the analytical rows gives these maxima of e_Z and e_Y over the FEM frequencies, with the frequency and the
entry of each maximum and the RMS over the frequencies:

| Case | e_Z | At | RMS | e_Y | At | RMS |
|---|---|---|---|---|---|---|
| `cable_30kv_630mm2_trefoil` | 5.21 % | 350 Hz, (1,1) | 3.34 % | 0.0209 % | all, (1,1) | 0.0209 % |
| `cable_30kv_na2xs2y_630mm2_trefoil` | 10.5 % | 10 MHz, (1,1) | 5.86 % | 16.7 % | 10 MHz, (4,4) | 15.9 % |
| `cable_132kv_630mm2_flathor` | 5.48 % | 10 MHz, (4,4) | 3.76 % | 91.4 % | 10 MHz, (5,5) | 90.6 % |
| `cable_18kv_1000mm2_trefoil` | 9.48 % | 10 MHz, (1,1) | 5.53 % | 86.9 % | 10 MHz, (5,5) | 86.0 % |
| `cable_380kv_2000mm2_flatver` | 5.60 % | 10 MHz, (4,4) | 3.74 % | 94.9 % | 10 MHz, (8,9) | 94.0 % |

The coaxial formulation takes a screen of wires with its copper tape as one tube. The tube shields the field between
the wires. The capacitance of the wires surrounded by the tape stays out of the analytical `Y`, while the FEM
keeps it. Below 10 kHz, this limitation of the coaxial formulation gives e_Y of 15.5 % in
`cable_30kv_na2xs2y_630mm2_trefoil`. In the three cases with an aluminum tape over the screen, it gives 90.4 %, 85.9 %
and 93.8 %.

The analytical rows correct the impedance of the helical core and screen wires for their lays, and the two-dimensional
FEM model has no lay. With every lay straight in the analytical rows (lay ratio 10⁶), the maximum e_Z of the two flat
formations is 1.81 % and 2.35 %. It lies in core entries between 25 Hz and 6.3 kHz. The maximum e_Z of the two
trefoils is 7.51 % for `cable_18kv_1000mm2_trefoil` and 9.00 % for `cable_30kv_na2xs2y_630mm2_trefoil`, at 10 MHz. A FEM run of
`cable_132kv_630mm2_flathor` with straight lays gives its `fem` rows again within 1.3e-8 of the largest diagonal entry,
at 0.1 Hz and 10 MHz.

In an analytical run of `cable_132kv_630mm2_flathor` with straight lays, the 61 core wires become one solid copper core
of their area, inside the same fill. Its e_Z exceeds the stranded core's at every FEM frequency from 0.63 Hz to
39.8 kHz, at most 2.19 % at 158 Hz.

Both trefoils keep their e_Z at 10 MHz from the proximity of the cables. The coaxial formulation distributes the
current of each conductor uniformly around its axis, and the FEM resolves the currents that the neighbors induce. With
the cables five diameters apart and straight lays, e_Z is at most 1.08 % in both, at 6.3 kHz.

For `cable_30kv_630mm2_trefoil`, the script also compares the positive-sequence parameters `z1 = zs - zm` and
`y1 = ys - ym`. FEM against the analytical rows, the resistance `R1` deviates by +5.57 % at 50 Hz, -12.7 % at 650 Hz and
-5.63 % at 2500 Hz. The reactance `X1` deviates by -7.93 % at 250 Hz to +1.42 % at 1950 Hz, and the susceptance `B1` by
-0.0209 % at every frequency. These deviations come from the proximity of the touching cables. With the cables five
diameters apart, e_Z is at most 0.155 %, and `R1` and `X1` lie within 0.24 % and 0.28 % of the analytical values.

## A sector core against its equivalent tube

One copper sector core is alone in its cable, 10 m above earth of 100 Ω·m (εr 10):
`Sector(0, 10.0e-3, 2π / 3, 1.0e-3 / sin(π / 3), 1.5e-3, Placement(0, 0))`. Its sides meet
1.15 mm from its center, and fillets of 1.5 mm round its three corners. Its area is 83.59 mm² and
its perimeter 35.51 mm.

Ametani's equivalent tube has the same area and perimeter. Its inner radius is 2.311 mm and its
outer radius 5.652 mm, with an insulator in its hole. The table compares the FEM of the core with
the coaxial formulation of the tube. The FEM of the tube against its coaxial formulation is the
control.

| Frequency | Core: e_Z | Core: e_Y | Control: e_Z | Control: e_Y |
|---|---|---|---|---|
| 0.1 Hz | 0.0625 % | 0.389 % | 0.0715 % | 0.765 % |
| 1 Hz | 0.0620 % | 0.387 % | 0.0876 % | 0.712 % |
| 10 Hz | 0.112 % | 0.387 % | 0.373 % | 0.701 % |
| 100 Hz | 0.243 % | 0.396 % | 0.596 % | 0.746 % |
| 1 kHz | 0.347 % | 0.386 % | 0.590 % | 0.736 % |
| 10 kHz | 0.355 % | 0.371 % | 0.580 % | 0.706 % |
| 100 kHz | 0.317 % | 0.210 % | 0.554 % | 0.653 % |
| 1 MHz | 0.492 % | 0.479 % | 0.846 % | 0.940 % |

The core and the tube agree within the FEM accuracy level at every frequency. The mesh of the core
has the node counts of its curves and the interior target of a conductor, without the graded
layer. At 0.1 Hz, the resistance of the core is 0.0624 % above `ρ/A` plus the earth return
resistance `ωμ₀/8`. Each computation ran once, at the commit that adds this section.

## Internal approximation

The formula `:wedepohl1973` approximates the surface impedances of `:schelkunoff1934` by the closed
forms of its paper. The copper sheath of the coaxial cases has inner radius 20.5 mm and outer radius
21.5 mm. Its relative deviation over the 101 frequencies of the analytical rows is:

| Surface | Maximum | At | RMS |
| --- | --- | --- | --- |
| `inner` | 1.55e-4 | 11 kHz | 5.95e-5 |
| `outer` | 1.49e-4 | 11 kHz | 5.68e-5 |
| `transfer` | 2.75e-4 | 10 MHz | 1.39e-4 |

The derivation of the paper takes the curvature of the sheath as constant across its wall.

The formula enters `Z` alone, and `Y` of each reference is unchanged to the bit. With
`internal_impedance = :wedepohl1973`, the largest scaled deviation of `Z` from the reference is:

| Case | Scaled difference in Z | At |
| --- | --- | --- |
| `two_bare_wires` | 1.43e-3 | 14 Hz |
| `two_insulated_wires` | 1.43e-3 | 14 Hz |
| `coaxial_circular_strands` | 1.10e-3 | 275 Hz |
| `coaxial_rectangular_strands` | 1.04e-3 | 229 Hz |
| `overhead_bare_wires` | 1.42e-3 | 14 Hz |
| `mixed_bare_wires` | 1.42e-3 | 14 Hz |

Each maximum lies near `|m r| = 3.5`, in the solid conductor. The paper gives that value for the largest
error of the reactive part of its form for a solid conductor.

## Earth-return formulas

Each earth-return formula of this section has one file in each case folder whose pairs it treats,
named after the formula. The file gives the matrix `Zg` of the formula, and `Pg` where the formula
gives it, in the form of `stages.csv`, at the same 101 frequencies. The reference problem of
`problems.jl` computed it with `trace = true`, at the commit that adds the file.

| Files | Formulation |
| --- | --- |
| `gary1976.csv` of `overhead_bare_wires` | `earth_impedance = :gary1976` |
| `wedepohl1973.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :wedepohl1973` |
| `saad1996.csv` of those four cases | `earth_impedance = :saad1996` |
| `carson1926.csv` of `overhead_bare_wires` | `earth_impedance = :carson1926` |
| `wise1934.csv` of `overhead_bare_wires` | `earth_impedance = :wise1934` |
| `wise1948.csv` of `overhead_bare_wires` | `earth_admittance = :wise1948` |
| `pollaczek1926.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :pollaczek1926` |
| `pollaczek1926.csv` of `mixed_bare_wires` | `earth_impedance = (air = :carson1926, earth = :pollaczek1926, mixed = :pollaczek1926)` |
| `papadopoulos2010.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :papadopoulos2010, earth_admittance = :papadopoulos2010` |
| `xue2018.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :xue2018, earth_admittance = :xue2018` |
| `xue2018-gamma.csv` of `two_insulated_wires` | those two slots, with the `Γ_k` of `stages-gamma.csv` |
| `lucca1994.csv` of `mixed_bare_wires` | `earth_impedance = (air = :gary1976, earth = :saad1996, mixed = :lucca1994)` |
| `ametani2009.csv` of `mixed_bare_wires` | `earth_impedance = (air = :gary1976, earth = :saad1996, mixed = :ametani2009)` |

These formulas give `Zg` alone, so each `Pg` and `Y` is that of the reference. The table gives the
deviation e_Z of `Z`. Its maximum and RMS over the 101 frequencies compare with the `analytical`
rows of `matrices.csv`. The last column compares with the `fem` rows over the 11 FEM frequencies:

| Formula | Case | Maximum | RMS | Maximum against FEM |
| --- | --- | --- | --- | --- |
| `:gary1976` | `overhead_bare_wires` | 0.0285 | 0.00748 | 0.0257 |
| `:wedepohl1973` | `two_bare_wires` | 62.2 | 12 | 62.3 |
| `:wedepohl1973` | `two_insulated_wires` | 61.5 | 11.9 | 61.5 |
| `:wedepohl1973` | `coaxial_circular_strands` | 0.105 | 0.0202 | 0.108 |
| `:wedepohl1973` | `coaxial_rectangular_strands` | 0.105 | 0.0203 | 0.109 |
| `:saad1996` | `two_bare_wires` | 0.568 | 0.126 | 0.567 |
| `:saad1996` | `two_insulated_wires` | 0.565 | 0.126 | 0.564 |
| `:saad1996` | `coaxial_circular_strands` | 0.058 | 0.0113 | 0.102 |
| `:saad1996` | `coaxial_rectangular_strands` | 0.0581 | 0.0113 | 0.105 |
| `:carson1926` | `overhead_bare_wires` | 0.0285 | 0.00422 | 0.0257 |
| `:wise1934` | `overhead_bare_wires` | 0.0285 | 0.00422 | 0.0257 |
| `:pollaczek1926` | `two_bare_wires` | 0.568 | 0.126 | 0.567 |
| `:pollaczek1926` | `two_insulated_wires` | 0.565 | 0.126 | 0.564 |
| `:pollaczek1926` | `coaxial_circular_strands` | 0.0599 | 0.011 | 0.101 |
| `:pollaczek1926` | `coaxial_rectangular_strands` | 0.06 | 0.011 | 0.104 |

Each deviation follows from the approximation of the formula:

- `:gary1976` replaces the interface integral by an image at a complex depth.
- `:wedepohl1973` keeps the low orders of an expansion in `m`, within `|m D| < 0.25` and
  `|m r| < 0.25` by its source.
- `:saad1996` approximates the interface integral in closed form, and it takes each cable as a
  filament of the exterior radius.
- `:carson1926` leaves out the displacement currents of the air and the earth. Its deviation grows
  with frequency.
- `:wise1934` is the integral of `:carson1926` at `eps_r = 1`, so the two share their cause.
- `:pollaczek1926` takes each cable as a filament and leaves out both displacement currents.

The implementation at commit `88f5de54` stops the integral at `λ = 1` 1/m and keeps the displacement
current of the air in its root. It differs from `:carson1926` by 0.0105 of the largest diagonal entry
of `Zg` at 13.2 kHz, where the stop of the integral acts. At 10 MHz, where the air acts, it differs by
0.0281.

A mixed formula gives the entries of a pair across the interface, so the air and earth places give
the other entries of `Zg`. The table gives e_Z of `mixed_bare_wires` for each mixed formula, under
two selections of those places. With `:martins2026` in them, e_Z is the deviation of the mixed
entries alone:

| Mixed place | Air and earth places | Maximum | RMS | Maximum against FEM |
| --- | --- | --- | --- | --- |
| `:lucca1994` | `:gary1976`, `:saad1996` | 0.0921 | 0.0502 | 0.0924 |
| `:lucca1994` | `:martins2026` | 0.0907 | 0.0427 | 0.0879 |
| `:ametani2009` | `:gary1976`, `:saad1996` | 0.0921 | 0.0305 | 0.0924 |
| `:ametani2009` | `:martins2026` | 0.0258 | 0.0122 | 0.0239 |
| `:pollaczek1926` | `:carson1926`, `:pollaczek1926` | 0.0921 | 0.0277 | 0.0924 |
| `:pollaczek1926` | `:martins2026` | 0.000341 | 0.000139 | 0.00613 |

At 631 kHz the mixed entry of `Zg` is `-4.61e-4 + 6.80e-4j` Ω/m under `:martins2026`,
`5.98e-2 + 6.35e-2j` under `:lucca1994` and `-3.03e-4 + 5.10e-4j` under `:ametani2009`. `:lucca1994`
and `:ametani2009` differ by at most 0.101 of the largest diagonal entry of `Zg`, at 13.2 kHz.

Each deviation of a mixed formula follows from its approximation:

- `:lucca1994` takes an image at a complex depth and one correction term. Its logarithm does not
  decay with depth.
- `:ametani2009` multiplies that logarithm by the decay `exp(-h_g / p)` over the depth `h_g`.

`:pollaczek1926` in the mixed place is the integral that the two others approximate. Against it,
`:lucca1994` differs by at most 0.0911 of the largest diagonal entry of `Zg`, at 22.9 kHz. `:ametani2009`
differs by at most 0.0258, at 2.51 kHz. At commit `88f5de54` the buried integral stops at `λ = 1` 1/m. On
`two_bare_wires` it differs from `:pollaczek1926` by at most 0.0125 of the largest diagonal entry of `Zg`,
at 9.12 kHz.

`:papadopoulos2010` gives `Zg` and `Pg`, so `Z` and `Y` both follow from it. The table gives e_Z and e_Y.
Their maximum and RMS over the 101 frequencies compare with the `analytical` rows of `matrices.csv`. The
last two columns compare with the `fem` rows over the 11 FEM frequencies:

| Case | e_Z maximum | e_Z RMS | e_Y maximum | e_Y RMS | e_Z against FEM | e_Y against FEM |
| --- | --- | --- | --- | --- | --- | --- |
| `two_bare_wires` | 0.568 | 0.126 | 1.95 | 0.374 | 0.567 | 0.788 |
| `two_insulated_wires` | 0.565 | 0.126 | 0.00202 | 0.000336 | 0.564 | 0.00214 |
| `coaxial_circular_strands` | 0.0778 | 0.0162 | 0.0668 | 0.0249 | 0.118 | 0.117 |
| `coaxial_rectangular_strands` | 0.0779 | 0.0162 | 0.0668 | 0.0249 | 0.121 | 0.118 |

`e_Z` follows from the filament and from the longitudinal constant. `e_Y` follows from the term `G` of
equation (6c) of the source. On `two_bare_wires`, that term of a diagonal entry is 8 to 9 times the rest
of the entry below 1 kHz. On `coaxial_circular_strands` it is 2.4 times the rest at 0.1 Hz. The source
writes the direct and image term as a spectral integral. Its cosine converges only conditionally at
equal depths, so the formula takes that term in closed form.

The logarithm in `G` stops at `|κ₂ σ̃₁/σ̃₂|` or at `k₀ √(ε_r - 1)`, the larger of the two, so `Pg` is
sensitive to `eps_r` near 1 at low frequency. At 0.1 Hz on `two_bare_wires` a diagonal entry moves by
11 % from `eps_r = 1` to `1 + 1e-9`, and by 39 % to `eps_r = 2`.

At commit `88f5de54` the integral of `Zg` stops at `λ = 1` 1/m, and the direct and image term has
`√(s μ₂ σ̃₂)` in place of `√(s μ₂ σ₂)`. The table gives the deviation of its `Zg` and `Pg` from those of
`:papadopoulos2010`, as a fraction of the largest diagonal entry of the matrix:

| Case | `Zg` maximum | At | `Pg` maximum | At |
| --- | --- | --- | --- | --- |
| `two_bare_wires` | 0.0125 | 9.12 kHz | 4.33e-5 | 10 MHz |
| `two_insulated_wires` | 0.0125 | 9.12 kHz | 4.41e-5 | 10 MHz |
| `coaxial_circular_strands` | 0.0539 | 10 MHz | 0.0589 | 10 MHz |
| `coaxial_rectangular_strands` | 0.0539 | 10 MHz | 0.0589 | 10 MHz |

`:wise1948` gives `Pg` alone, so `Z` is that of the reference. The table gives the deviation e_Y of `Y`.
Its maximum and RMS over the 101 frequencies compare with the `analytical` rows of `matrices.csv`. The
last column compares with the `fem` rows over the 11 FEM frequencies:

| Case | Maximum | RMS | Maximum against FEM |
| --- | --- | --- | --- |
| `overhead_bare_wires` | 0.0202 | 0.00297 | 0.0174 |

`:wise1948` is another quantity than `:martins2026`. It is a scalar potential coefficient under a
longitudinal constant fixed at that of the air. `:martins2026` gives the mean of the integrated vertical
field referred to the interface at `Γ = 0`.

Each corrects the ideal images. Over 100 Ω·m with `eps_r = 10`, the ratio of the correction of
`:martins2026` to that of `:wise1948` is -0.0235 - 0.00190j at 0.1 Hz. At 10 MHz it is -0.297 - 0.839j. At low
frequency they have opposite signs, and that of `:wise1948` is the larger.

`:xue2018` gives `Zg` and `Pg`. The table compares them with `stages.csv`, and on its last row with
`stages-gamma.csv`, as a fraction of the largest diagonal entry of the matrix at each frequency. It gives
the maximum and RMS over the 101 frequencies:

| Case | `Zg` maximum | `Zg` RMS | `Pg` maximum | `Pg` RMS |
| --- | --- | --- | --- | --- |
| `two_bare_wires` | 0.568 | 0.126 | 0.568 | 0.126 |
| `two_insulated_wires` | 0.582 | 0.129 | 0.582 | 0.129 |
| `coaxial_circular_strands` | 0.00201 | 0.000364 | 0.00198 | 0.000361 |
| `coaxial_rectangular_strands` | 0.00201 | 0.000364 | 0.00198 | 0.000361 |
| `two_insulated_wires`, `Γ_k` | 0.582 | 0.129 | 0.582 | 0.129 |

The deviation follows from the finite radius. `:xue2018` takes each cable as a filament, and
`:martins2026` closes the current on its circumference. At low frequency the deviation falls with the
square of the radius times the square of the filament's logarithm, in `Zg` and in `Pg`. So `Pg` refers
to the same infinite depth as `:martins2026`.
