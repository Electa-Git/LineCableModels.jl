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
diagonal entry. Each `physical.geo` is byte-identical to the retained file of its case's run of
2026-10-25. The geometry tests compare against these files.

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
| `coaxial_rectangular_strands` | 5.22 % | 10.4 % |
| `overhead_bare_wires` | 0.610 % | 0.817 % |
| `mixed_bare_wires` | 0.609 % | 0.617 % |

For `two_bare_wires`, `two_insulated_wires`, `overhead_bare_wires` and
`mixed_bare_wires`, the meshes coincide with the source meshes.
Their rows agree with the frozen `fem` rows within 1e-6 scaled.
For the coaxial cases, the maxima against the frozen `fem` rows are:

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
| `pollaczek1926.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :pollaczek1926` |
| `pollaczek1926.csv` of `mixed_bare_wires` | `earth_impedance = (air = :carson1926, earth = :pollaczek1926, mixed = :pollaczek1926)` |
| `papadopoulos2010.csv` of `two_bare_wires`, `two_insulated_wires` and the two coaxial cases | `earth_impedance = :papadopoulos2010, earth_admittance = :papadopoulos2010` |
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

At commit `88f5de54` the integral of `Zg` stops at `λ = 1` 1/m, and the direct and image term has
`√(s μ₂ σ̃₂)` in place of `√(s μ₂ σ₂)`. The table gives the deviation of its `Zg` and `Pg` from those of
`:papadopoulos2010`, as a fraction of the largest diagonal entry of the matrix:

| Case | `Zg` maximum | At | `Pg` maximum | At |
| --- | --- | --- | --- | --- |
| `two_bare_wires` | 0.0125 | 9.12 kHz | 4.33e-5 | 10 MHz |
| `two_insulated_wires` | 0.0125 | 9.12 kHz | 4.41e-5 | 10 MHz |
| `coaxial_circular_strands` | 0.0539 | 10 MHz | 0.0589 | 10 MHz |
| `coaxial_rectangular_strands` | 0.0539 | 10 MHz | 0.0589 | 10 MHz |
