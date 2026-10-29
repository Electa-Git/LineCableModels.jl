# Line parameters

```@meta
CurrentModule = LineCableModels
```

In a computation, the problem states what to compute and the formulation selects the formulas.
A backend executes the computation, and `compute(problem, formulation; kw...)` is the one
call.

## The problem

[`LineParametersProblem`](@ref)`(system, earth; temperature, frequencies)` describes the
computation by the cables of a [`CableSystem`](@ref) and the [`EarthModel`](@ref) below the
air, with the temperature of every conductor and dielectric in °C, default 20, and the
frequencies of the scan in Hz. The frequencies must be positive, finite and in strict order from low to high.
The problem promotes every number of the system, the earth, the temperature and the
frequencies with `Float64`, so its number type `T` is at least `Float64`: a `Float32` input
computes in `Float64`.

## The coaxial formulation

[`CoaxialFormulation`](@ref) selects one formula for each physical contribution. Each slot
takes a formula of its family, and its keyword accepts `:id`, `:id => (key = value, ...)` with
the formula's options, or a formula. The defaults:

| Keyword | Family | Default | Contribution |
| --- | --- | --- | --- |
| `internal_impedance` | `InternalImpedance` | `:schelkunoff1934` | surface impedances of the conductors |
| `insulation_impedance` | `InsulationImpedance` | `:ametani1980` | magnetic series term of a dielectric interval |
| `insulation_admittance` | `InsulationAdmittance` | `:lossless` | admittivity of insulator layers |
| `semicon_admittance` | `SemiconAdmittance` | `:lossless` | admittivity of semicon layers |
| `earth_impedance` | `EarthReturn` | `:unified` | earth-return impedance |
| `earth_admittance` | `EarthReturn` | `:unified` | earth-return potential coefficients |
| `shunt` | `Shunt` | `:series` | potential coefficients inside the cables |
| `resistivity` | `Resistivity` | `:linear` | resistivity of the conductors at the temperature |
| `earth_properties` | `EarthProperties` | `:constant` | properties of the air and the earth |

When both earth slots take one identifier, their formulas must be equal. The unified formula
gives the impedance and the potential coefficients from the same circumference currents. The
other keywords are the options of the formulation:

| Option | Default | Meaning |
| --- | --- | --- |
| `Γ` | `0` | longitudinal wavenumber in 1/m, one value or one per frequency |
| `bundle` | `true` | merge the conductors that share a positive connection index |
| `kron` | `true` | remove the conductors of connection index zero by [`eliminate!`](@ref) |
| `transposition` | `false` | apply the ideal cyclic transposition |

The formulas are on the theory pages [Matrix formulation](theory/matrix.md),
[Conductors](theory/conductors.md), [Dielectrics](theory/dielectrics.md) and
[Earth return](theory/earth.md).

## The computation

`compute(problem, formulation; kw...)` runs the default backend of the formulation type, which
the trait [`backend`](@ref) gives: [`LCMBackend`](@ref), the computation of this package, for
[`CoaxialFormulation`](@ref), and [`GmshBackend`](@ref) for [`FEM`](@ref). A backend passed
first, `compute(backend, problem, formulation; kw...)`, chooses another. The keywords are the
computation options of the backend. [`LCMBackend`](@ref) has one, `trace`, default `false`.

At each frequency `f`, with `s = j2πf`, the computation of [`LCMBackend`](@ref):

1. It evaluates the material laws: the resistivity of each conductor at the temperature, the
   admittivity of each dielectric layer and the properties of the air and the earth.
2. It computes the surface impedances of the conductors and the magnetic terms of the
   dielectric intervals, and the loop recursion assembles them into the impedance matrix
   inside the cables, `Zin`.
3. It builds the potential coefficients inside the cables, `Pin`, by the series chain of
   their layers.
4. It computes the earth-return impedance `Zg` and potential coefficients `Pg` of the cables
   from their outer surfaces.
5. It adds the entry of `Zg` and `Pg` for each pair of cables, a cable with itself included,
   to every pair of their conductors.
6. It applies the bundle and Kron reductions to the impedance and potential matrices.
7. It gives the admittance `Y = s P⁻¹` of the reduced potential matrix `P`.
8. It applies the ideal transposition to `Z` and `Y` when `transposition` is `true`.

The docstring of `solve!` in the [Matrix formulation](theory/matrix.md) page states the same
steps.

## The result

[`LineParameters`](@ref) holds the result:

| Field | Meaning | Unit |
| --- | --- | --- |
| `Z` | series impedance by row, column and frequency | Ω/m |
| `Y` | shunt admittance by row, column and frequency | S/m |
| `frequencies` | frequencies of the scan | Hz |
| `terminals` | label of each row | |
| `line_length` | length of the line | m |
| `formulation` | formulation that computed the result | |
| `trace` | matrices from which the backend builds `Z` and `Y` | |

A row label is `cable:<i>:<terminal>`, the index of the cable in the system and the name of the
terminal. A bundle's row joins the labels of its terminals with `+`. With `trace = true`, the
trace of [`LCMBackend`](@ref) holds `Zin` in Ω/m and `Pin` in m/F by conductor, and `Zg` in
Ω/m and `Pg` in m/F by cable, each with one slice per frequency. Without it, the four arrays
are empty.

## Options

Each owner of options declares them with [`defaults`](@ref): a formula or a formulation type
for [`FormulationOptions`](@ref), a backend type for [`ComputationOptions`](@ref). The record
is built once from the defaults and the keywords given. An unknown key throws an `ArgumentError`,
each value takes the type of its default, and [`validate`](@ref) checks the record.

```@example computing
using LineCableModels
LineCableModels.defaults(CoaxialFormulation)
```

## The references

The case scripts of the branch `release/v0.2.0-fem-investigation` computed the reference
systems in `test/references/`:

- `two_bare_wires`: two copper wires of radius 42.5 mm, 1 m apart, buried 1 m deep in earth of
  0.1 Ω·m.
- `two_insulated_wires`: the same wires with 1 mm of insulation.
- `coaxial_circular_strands`: one coaxial cable with a core of circular strands, buried 1 m
  deep in earth of 100 Ω·m.
- `coaxial_rectangular_strands`: the same cable with a core of rectangular strands.
- `overhead_bare_wires`: the bare wires 1 m above earth of 0.1 Ω·m.
- `mixed_bare_wires`: one wire overhead and one buried.

For each case the tests compare `Z` and `Y` with the frozen analytical rows within 1e-9 of the
largest diagonal entry at each frequency, `Zin` and `Pin` within 1e-12 and `Zg` and `Pg`
within 1e-9.
