# Line parameters

```@meta
CurrentModule = LineCableModels
```

A computation of line parameters is composed of a problem, which defines what is computed,
and a formulation, which selects the formulas. The computation is executed by a backend, and
all computations are invoked by the same call, `compute(problem, formulation; kw...)`.

## Problem definition

The problem [`LineParametersProblem`](@ref)`(system, earth; temperature, frequencies)` is
defined by the cables of a [`CableSystem`](@ref), the [`EarthModel`](@ref) below the air, the
temperature of all conductors and dielectrics in °C, equal to 20 by default, and the
frequencies of interest in Hz. The frequencies must be positive, finite and sorted in strictly
increasing order. All numbers of the system, earth, temperature and frequencies are promoted
together with `Float64`, so that the number type `T` of the problem is at least `Float64`.

## Coaxial formulation

The [`CoaxialFormulation`](@ref) selects one formula for each physical contribution. Each
slot takes a formula of its family, and the corresponding keyword accepts `:id`,
`:id => (key = value, ...)` with the options of the formula, or a formula object. The
default formulas are summarized below.

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

When both earth slots take the same identifier, their formulas must be equal, since the
unified formula obtains the impedance and the potential coefficients from the same
circumferential currents. The remaining keywords are the options of the formulation:

| Option | Default | Meaning |
| --- | --- | --- |
| `Γ` | `0` | longitudinal wavenumber in 1/m, one value or one per frequency |
| `bundle` | `true` | merge the conductors that share a positive connection index |
| `kron` | `true` | remove the conductors of connection index zero by [`eliminate!`](@ref) |
| `transposition` | `false` | apply the ideal cyclic transposition |

The formulas are described in the theory pages [Matrix formulation](theory/matrix.md),
[Conductors](theory/conductors.md), [Dielectrics](theory/dielectrics.md) and
[Earth return](theory/earth.md).

## Computation

The call `compute(problem, formulation; kw...)` runs the default backend of the formulation
type, given by the trait [`backend`](@ref): [`LCMBackend`](@ref), the analytical
computation of this package, for [`CoaxialFormulation`](@ref), and [`GmshBackend`](@ref) for
[`FEM`](@ref). Another backend is selected by passing it as the first argument,
`compute(backend, problem, formulation; kw...)`. The keywords are the computation options of
the backend, and [`LCMBackend`](@ref) has one option, `trace`, equal to `false` by default.

At each frequency `f`, with `s = j2πf`, the computation of [`LCMBackend`](@ref) consists of
the following steps:

1. Evaluation of the material laws, which give the resistivity of each conductor at the
   operating temperature, the admittivity of each dielectric layer and the properties of the
   air and the earth.
2. Computation of the surface impedances of the conductors and of the magnetic terms of the
   dielectric intervals, which are assembled by the loop recursion into the impedance matrix
   inside the cables, `Zin`.
3. Computation of the potential coefficients inside the cables, `Pin`, by the series
   combination of their layers.
4. Computation of the earth-return impedance `Zg` and potential coefficients `Pg` of the
   cables, from their outer surfaces.
5. Addition of the entries of `Zg` and `Pg` for each pair of cables, including each cable
   with itself, to every pair of their conductors.
6. Bundle and Kron reductions of the impedance and potential matrices.
7. Computation of the admittance `Y = s P⁻¹` from the reduced potential matrix `P`.
8. Ideal transposition of `Z` and `Y`, when `transposition` is `true`.

The same steps are given in the docstring of `solve!`, in the
[Matrix formulation](theory/matrix.md) page.

## Results

The result of the computation is a [`LineParameters`](@ref) object, with the following
fields:

| Field | Meaning | Unit |
| --- | --- | --- |
| `Z` | series impedance by row, column and frequency | Ω/m |
| `Y` | shunt admittance by row, column and frequency | S/m |
| `frequencies` | frequencies of interest | Hz |
| `terminals` | label of each row | |
| `line_length` | length of the line | m |
| `formulation` | formulation used in the computation | |
| `trace` | matrices from which the backend assembles `Z` and `Y` | |

Each row is labeled as `cable:<i>:<terminal>`, where `i` is the index of the cable in the
system and `terminal` is the name of the terminal. The label of a bundle joins the labels of
its terminals with `+`. With `trace = true`, the trace of [`LCMBackend`](@ref) contains `Zin`
in Ω/m and `Pin` in m/F per conductor, and `Zg` in Ω/m and `Pg` in m/F per cable, each with
one slice per frequency. Otherwise, the four arrays are empty.

## Options

The options of each owner are declared by [`defaults`](@ref), where the owner is a formula or
a formulation type for [`FormulationOptions`](@ref), and a backend type for
[`ComputationOptions`](@ref). The record of options is built once, from the defaults and the
keywords given. An unknown key throws an `ArgumentError`, each value is converted to the type
of its default, and the record is checked by [`validate`](@ref). The defaults of the coaxial
formulation are:

```@example computing
using LineCableModels
LineCableModels.defaults(CoaxialFormulation)
```

## Reference systems

The test suite includes six reference systems, whose line parameters computed by the
analytical and finite-element formulations are stored in `test/references/`:

- `two_bare_wires`: two copper wires with a radius of 42.5 mm, buried at 1 m depth with a
  horizontal spacing of 1 m, in a homogeneous earth of 0.1 Ω·m.
- `two_insulated_wires`: the same wires with an insulation 1 mm thick.
- `coaxial_circular_strands`: one coaxial cable with a core of circular strands, buried at 1 m
  depth in a homogeneous earth of 100 Ω·m.
- `coaxial_rectangular_strands`: the same cable with a core of rectangular strands.
- `overhead_bare_wires`: the two bare wires at a height of 1 m above an earth of 0.1 Ω·m.
- `mixed_bare_wires`: the two bare wires with one placed above the earth and the other
  buried.

The validation consists of comparing the analytical and finite-element formulations on these
systems. The analytical results of the package agree with the reference values of `Z` and
`Y` within 1e-9 of the largest diagonal entry at each frequency, with `Zin` and `Pin` within
1e-12 and `Zg` and `Pg` within 1e-9. The finite-element results are compared with the
analytical ones by means of the relative errors defined in [Finite element solver](fem.md#Accuracy),
at an accuracy level of 1 %.
