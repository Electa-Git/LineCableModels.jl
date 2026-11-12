# Architecture

This page is the architectural specification of LineCableModels. It states the structure
that the code implements and that every change preserves. The [rules and
patterns](rules.md) say how to extend the package within that structure, and the
[conventions](conventions.md) say how the code and its documentation are written.

The package computes the series impedance and the shunt admittance of cable systems per unit
length. A problem describes the cables, the earth, the temperature and the frequencies. The
formulation selects the formula of each step, and the backend executes the computation.

## Three axes

The problem, the formulation and the backend answer three different questions. They are
orthogonal.

- **The problem says what is computed.** `LineParametersProblem` combines a cable system,
  an earth model, the operating temperature and the frequencies of the scan.
- **The formulation says how.** It is a methodology: one formula for each step, with the
  controls of the methodology. `CoaxialFormulation` is the methodology of concentric tubes.
  `FEM` is the methodology of a field solve with the material laws.
- **The backend says who computes.** `LCMBackend` is this package's own implementation, in
  process. `GmshBackend` runs the Gmsh and GetDP processes of the extension.

To solve a problem, a backend implements a formulation, a methodology with its own formula
for each step. A formulation is the same object on every backend. The result records that
same object. A backend implements a formulation by a `compute` method. It implements each
formula of the formulation by the methods of its own route. The trait `backend` gives the
default backend of a formulation type and nothing more.

```julia
result = compute(problem, CoaxialFormulation())               # on LCMBackend, the default
result = compute(problem, FEM(); workers = 4)                  # on GmshBackend, after `using Gmsh`
result = compute(LCMBackend(), problem, CoaxialFormulation(); trace = true)
```

A backend refuses before it runs what it cannot honor: a formula it has no route for, or an
option it cannot apply. Results never claim a computation that did not take place.

Controls live at two levels. A formula's numerical controls belong to the formula. A
methodology's controls belong to the formulation. A backend's process controls belong to the
backend. A control belongs to whoever makes the approximation it controls.

## Modules and ownership

The package is six modules in one load order, and four extensions.

| Module | Owns |
| --- | --- |
| `Commons` | the grammar: `validate`, `resolve`, `build`, `compute`, `defaults`, `select`, the abstract types `AbstractFormulation` and `AbstractBackend`, the option records |
| `Materials` | `Material{T, K}`, the properties of a material with its class `K` as a type parameter |
| `Earth` | `EarthLayer`, `EarthModel`, `homogeneous`: the air and the earth as a stack of layers. `EquivalentEarthLayer` is the medium that a model computes for the layers at a frequency. |
| `BaseGeometry` | the plane geometry of the placement, the segments, arcs and elliptic arcs, the contour and every method on them, and the verb `draw`, which walks a contour to its curves on a target |
| `DataModel` | the cross-sections of a cable, each with its contour, the declared parts, the cable design, the cable system, the blueprint and the builders of the surface syntax |
| `Engine` | the problem, the formulations, the backends, the formula families, the numerical input, the buffers, the loop, the reductions, the spectral tools and the result |
| `LineCableModelsGmshExt` | the finite-element backend, with its model, geometry, meshes, GetDP scan and result |
| `LineCableModelsMeasurementsExt` | the propagation of input uncertainty: methods only |
| `LineCableModelsForwardDiffExt` | the derivatives with respect to the inputs: methods only |
| `LineCableModelsIntervalArithmeticExt` | the enclosures of the Bessel functions on an interval argument: methods only |

The formula families are the eight child modules of `Engine`: `InternalImpedance`,
`InsulationImpedance`, `InsulationAdmittance`, `SemiconAdmittance`, `EarthReturn`, `Shunt`,
`Resistivity` and `EarthProperties`.

The structure follows four rules.

- The root module `LineCableModels` defines nothing. It includes the modules and imports
  their names. It exports or declares public the names users see.
- A module references only the modules before it, its ancestors and its descendants.
- Each concept has one owner, one name and one implementation. Each name has one meaning.
- A module imports each foreign name explicitly from its owner. Its owner exports or
  declares public a name that other modules use.

A method that extends a function of another module mentions a type of its own module.
Constructors belong to the owner of their type. Every file of a module lies under the
directory of its module file. No core method is defined under `ext/`.

## Grammar

`Commons` defines the verbs that the other modules extend with methods.

| Function | Meaning |
| --- | --- |
| `validate(subject, context...)` | reject an invalid input before its consumer uses it, and return `subject` unchanged |
| `resolve(context, declaration)` | give a declared part or layer its geometry in its context, or place a cross-section in a frame |
| `build(Target, declarations...)` | construct one completed object of `Target` from its declarations |
| `compute(problem, formulation; kw...)` | run one computation of a problem by a formulation on the backend that runs it |
| `defaults(owner)` | the declared option keys and default values of a formula, a formulation or a backend |
| `select(default, value)` | convert an option to its declared type, or lift a formula selection into its slot |

The grammar belongs to `Commons` however many modules extend it. Any other definition in
`Commons` is user API, or it is used by methods of two owners outside `Commons`. Each
definition of `Commons` is exported or public and has a docstring, and its name appears in a
test under `test/commons/`.

The options are stored in records of two types. `FormulationOptions{O, NT}` is the record of
a formula or a formulation `O`. `ComputationOptions{B, NT}` is the record of a backend `B`.
Each record is built once, at construction, from `defaults(O)` and the keywords given. An
unknown key throws. Each value takes its declared type. `validate` checks the values. The
type of a record does not depend on which keys the caller passed, so one computation
compiles one specialization.

## Data model

The data model has three layers. It encodes a cable as a declaration and describes the cable
as resolved regions with tables. From the description it draws the equivalences that the
formulations read. It never rejects a physically valid design for a missing formula. A
formula that cannot treat an input throws when it is called.

### Cross-sections

A cross-section is one type with the methods of one unary interface: `area`, `centroid`,
`support`, `distance`, `perimeter`, `boundary`, `resolve` and its exact outline, `Contour`.
`r_in` and `r_ex` exist on circular cross-sections only. No method takes two cross-sections.
The primitives are `Circle`, `Annulus`, `Rectangle`, `Sector`, `Ellipse` and `Polygon`.
`ArcPolygon` is the region inside a closed path of segments and of circular or elliptic
arcs. `Offset(t)` is an operator: a layer declared by
its thickness, resolved on the outline of the part inside it. The layer is an `Annulus`
around a circle and a `PerforatedShape` around another outline. `PerforatedShape` is the
cross-section of a fill, a container without the outlines of its contents, or of such a
layer. A new kind of
cross-section is one type and its methods. The structure above it does not change.

### Declaration

A cable is a tree of parts. Each node has one job.

| Node | Job |
| --- | --- |
| `Region(tag, primitive, material)` | one material in one cross-section, named by its tag |
| `Stack(parts...)` | outward layering: each part resolves against the outline of the part inside it |
| `Terminal(name, part)` | every conductive region of `part` belongs to the terminal `name` |
| `Repeat(pattern, part)` | copies of `part` at the placements of a pattern, such as `CircularPattern` |
| `Lay(helix, repeat)` | each copy follows a `Helix` of a `LayRatio` about the center of the repeat |
| `Enclosure(fill, contents)` | a fill region, whose primitive is the container, around the contents |

The material's class tells whether a region conducts. A conductive region lies in a
terminal. A region of an insulator or a semicon belongs to no terminal. A terminal name
is an identification and nothing else. The circuit role of a terminal belongs to the
connections of the system.

### Description

`build(CableDesign, id, parts...)` resolves the tree into a `CableDesign`. The design stores
the placed regions from the inside out, the outline, and the tables of terminals,
sub-cables, repeats, lays and fills. A table has a row only for a fact that exists. The design's
number type `T` is the promotion of the number types of its cross-sections and its
materials. An uncertain material widens `T`. The design validates itself at construction.
The regions cover the outline without gaps or overlaps. Every region of an enclosure lies
inside its container. Every conductive region belongs to a terminal, and every region
belongs to one sub-cable.

### Equivalence

`flatten(design)` builds a `CableBlueprint` for the coaxial formulations. Each terminal
becomes an equivalent tube. A course of laid wires becomes a tube whose material reproduces
the direct-current resistance and the geometric mean radius of the wires. The overlength of
the helix enters the resistivity. The solenoid effect of the helix enters the relative
permeability of the dielectric interval, `mu_r_ins`. The dielectric layers between two
consecutive conductors form the interval of the inner one. The regions of an assembly
around its members form the interior of that assembly. The lay correction is made in the
equivalence. The earth-return formulas never see it.

### System

`build(CableSystem, designs, placements; connections, id, line_length)` places cables in one
cross-section. `connections` is one named tuple per cable, keyed by terminal name. An index
`k ≥ 1` is a live conductor, a repeated index bundles conductors, and `0` eliminates one.
Every terminal name must appear. The system rejects overlapping cables and keeps the
declared placements. A cable may touch the interface between air and earth and may never
cross it. `:martins2026` and the finite-element model refuse a cable that touches it.

### Surface syntax

`@design` and `@terminal` and the builders `core`, `insulation`, `screen`, `sheath`,
`bedding`, `jacket`, `stranded`, `wires` and `tape` lower to the nodes above, and `trefoil`
and `flat` return the placements of the cables of a system. They do not add a capability.

```julia
copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
polyethylene = Material(:insulator, 1.0e14, 2.3)
wire = @design "insulated_wire" begin
    @terminal :core begin
        core(copper; r = 4.0e-3)
    end
    insulation(polyethylene; t = 1.0e-3)
end
```

## Engine

### Formula families

A formula is one object per family: `Family.Formula{ID, O}`. The identifier `ID` is a
`Symbol` in the type. The options `O` are a `FormulationOptions{Formula{ID}}` record. The
family module is the context: `EarthReturn.Formula{:martins2026}` and
`InternalImpedance.Formula{:schelkunoff1934}` never meet. Each family exports `Formula` and
`FORMULAS`, the tuple of its identifiers. The unknown-identifier error prints that tuple.

A new formula is three parts and nothing else. It adds a `defaults` method for its options,
a line of the registry in `names.toml` mirrored by `FORMULAS`, and its physics methods.
`FormulaDefinition`, `Expression` and `Functor` do not exist. `:default` is never an
identifier: it is the one router symbol, which `select` resolves at construction to the slot's
default formula.

The physics of a formula is its callable. Each mathematical expression is a method of the
formula's own type. A literal tag selects the expression when the formula has more than one. A
component of one expression is a field of the record it returns.

A formula's programmatic methods take the formula object first: `defaults`, `domain`,
`validate`, `initialize_buffers`, `to_buffer!` and `from_buffer`. `to_buffer!` moves the formula's
per-frequency terms into its block of the buffers and does nothing by default. `from_buffer`
returns the block. The formulation's methods of `impedance!` and `potential!` assemble the
matrices of a computation. A formula does not write a matrix of the computation.

`domain` gives, by key, the closed interval of each quantity inside the formula's range of
validity. The fallback `(;)` in `Engine` stands for a formula without bounds.

| Family | Formula | Callable |
| --- | --- | --- |
| `InternalImpedance` | `:schelkunoff1934`, `:wedepohl1973` | `formula(r_in, r_ex, ρ, μ_r, s)`, the inner, outer and transfer surface impedances |
| `InsulationImpedance` | `:ametani1980` | `formula(r_in, r_ex, μ_r, s)`, the magnetic term of an annular dielectric, and `formula(n, r_in, r_ex, μ_r, s)`, its surfaces at an order |
| `InsulationAdmittance` | `:lossless`, `:lossy` | `formula(material, f)`, the admittivity of an insulator |
| `SemiconAdmittance` | `:lossless`, `:lossy` | `formula(material, f)`, the admittivity of a semicon |
| `EarthReturn` | `:martins2026` | the kernels `formula(Val(:J), λ, medium, pair)` and the entries `formula(Val(:Z), Val(:self), i, ...)` |
| `EarthReturn` | `:ideal` | the entries `formula(Val(:P), Val(:self), i, ...)` of an overhead pair |
| `EarthReturn` | `:gary1976`, `:wedepohl1973`, `:saad1996` | the entries `formula(Val(:Z), Val(:self), i, ...)` of the domains each treats |
| `EarthReturn` | `:lucca1994`, `:ametani2009` | the entry `formula(Val(:Z), Val(:mutual), i, j, ...)` of a mixed pair |
| `EarthReturn` | `:carson1926` | the kernel `formula(Val(:J), λ, γ², H, y)` and the entries `formula(Val(:Z), Val(:self), i, ...)` |
| `EarthReturn` | `:pollaczek1926` | the kernels and the entries of its two expressions, `formula(Val(:Z), Val(:earth), i, j, ...)` and `formula(Val(:Z), Val(:mixed), i, j, ...)` |
| `EarthReturn` | `:papadopoulos2010` | the kernels `formula(Val(:J), λ, ...)` and `formula(Val(:P̃), λ, ...)`, and the entries `formula(Val(:Z), Val(:self), i, ...)` and `formula(Val(:P), Val(:self), i, ...)` |
| `EarthReturn` | `:wise1934` | the kernel `formula(Val(:J), λ, κ₁², κ₂², H, y)` and the entries `formula(Val(:Z), Val(:self), i, ...)` |
| `EarthReturn` | `:wise1948` | the kernel `formula(Val(:P̃), λ, κ₁², κ₂², σ̃₁, σ̃₂, H, y)` and the entries `formula(Val(:P), Val(:self), i, ...)` |
| `EarthReturn` | `:xue2018` | the kernels `formula(Val(:J), λ, ...)`, `formula(Val(:M), λ, ...)` and `formula(Val(:P̃), λ, ...)`, and the entries `formula(Val(:Z), Val(:self), i, ...)` and `formula(Val(:P), Val(:self), i, ...)` |
| `Shunt` | `:series` | `formula(r_in, r_ex, κ, s)`, the coefficient of one layer, and `formula(n, r_in, r_ex, κ, s)`, its surfaces at an order, with the chain `potential!(formula, Pin, buffers, input, kappa, s)` |
| `Resistivity` | `:linear` | `formula(material, θ)` |
| `EarthProperties` | `:constant` and the nine soil laws | `formula(layer, f)`, a law: the layer at a frequency |
| `EarthProperties` | `:bottommost`, `:martins2020` and `:xue2021` | `formula(layers, f)`, a model: the equivalent earth of the layers |

### Formulation

`CoaxialFormulation` is a struct of ten typed slots, the tuple of its distinct formulas and
one options record. A slot accepts `:id`, `:id => (key = value, ...)`, a `Formula` of its family
or `:default`. The keyword constructor lifts each selection once, through the slot-tagged
`select`, and one method of `select` per slot names its default formula. The options are `Γ`,
the longitudinal wavenumber the computation assumes, and the reductions `bundle`, `kron` and
`transposition`.

```julia
CoaxialFormulation(; internal_impedance = :default, insulation_impedance = :default,
    insulation_admittance = :default, semicon_admittance = :default,
    earth_impedance = :default, earth_admittance = :default, shunt = :default,
    resistivity = :default, earth_properties = :default, equivalent_earth = :default,
    Γ = complex(0.0), bundle = true, kron = true, transposition = false)
```

Both earth slots are typed by one family, `EarthReturn`. Each records one formula per domain of
a pair of conductors:

- `air`: both conductors lie above the interface.
- `earth`: both conductors lie below the interface.
- `mixed`: one conductor lies above the interface and the other below it.

One selection fills the three places. A selection per domain,
`earth_impedance = (air = :martins2026,)`, fills the places it names, and the others take the
slot's default. A formula declares the pairs it treats. The computation refuses a formula in a
place whose pairs it does not treat, before it runs.

`@formulation` writes the keywords of a constructor as the lines of a block, and the tags `@air`,
`@earth` and `@mixed` write a selection per domain. They lower to the constructor and do not add
a capability. The constructor refuses a line for a slot or option that the formulation lacks.

A formula of the family gives the impedance, the potential coefficients, or both, by the
methods it defines. A formula in a slot whose quantity it does not give is refused at
construction. Formulas of one identifier among the six places must be equal, because a formula
that gives both quantities is one computation. The constructor stores each distinct formula
once in the field `formulas`, so `to_buffer!` runs once per formula and frequency. At a nonzero
`Γ`, a computation warns for each formula that does not read it: its `domain` declares
`Γ = 0.0 .. 0.0`.

The slots `earth_properties` and `equivalent_earth` hold formulas of one family, of two kinds.
The law describes each layer of the earth model at a frequency. Then the model reduces the
described layers below the air to the one medium that the earth-return formulas read, an
`EquivalentEarthLayer`. The computation builds that medium at every frequency before the loop,
and the air passes as declared. A law in the model's slot and a model in the
law's are refused at construction.

`FEM` is the formulation of the field solve. It keeps the five slots that the field solve
reads, `insulation_admittance`, `semicon_admittance`, `resistivity`, `earth_properties` and
`equivalent_earth`, the tuple of those formulas and the same options. Before it runs, the
computation validates the problem against each formula, as on the coaxial path. Every control
of the field solve belongs to the Gmsh backend.

### Backends and dispatch

The backends are singletons of `AbstractBackend`. `compute` has three methods in `Engine`.

```julia
compute(problem, formulation; kw...)                       # the default backend of the formulation type
compute(backend, problem, formulation; kw...)              # a chosen backend, the options built once
compute(backend, problem, formulation, options)            # the fallback: "does not implement"
```

`LCMBackend` computes a `CoaxialFormulation` in process. Without Gmsh, the options of
`GmshBackend` refuse the computation: `defaults(::Type{<:GmshBackend})` throws in `Engine`.
The Gmsh extension adds `defaults(::Type{GmshBackend})`, more specific than the refusal, and
the `compute` method of `GmshBackend` for `FEM`. Loading Gmsh replaces an error by a computation
without overwriting a method. A later backend adds a singleton, a `defaults` method, a
`compute` method and its refusals. The formulation types stay as they are.

### One computation

The scalar `compute` of `LCMBackend` depends on nothing above it. A parametric layer, when
it arrives, calls it per point, never the reverse. Its steps are fixed.

1. `validate(problem, formulation)`: the options against the frequencies and
   `validate(problem, formula)` for each formula.
2. `CoaxialInput(problem, options)` builds the blueprint of every cable, the placements and
   the reduction plan from `reductions(system, options)`.
3. `validate(input, formulation)`: the formula of each earth place against the pairs of layers
   of the input.
4. `validate(presented, formulation)`: one warning for each quantity of the problem and the
   input outside the `domain` of a formula.
5. `initialize_buffers(formulation, problem, input)`: every array once, before the loop.
6. The loop `solve!`, its own method, over the frequencies.
7. `LineParameters` from the filled arrays, validated once.

At each frequency the loop runs these steps in order.

1. The admittivity of each layer, the described layers of the earth model, their equivalent
   earth, and `to_buffer!` of the formulation.
2. The series matrix inside each cable, by the loop recursion of the surface impedances.
3. The potential matrix inside each cable, by the series chain of the layers.
4. `to_buffer!` of each distinct formula, then `impedance!` and `potential!` of the
   formulation, which fill `Zg` and `Pg` pair by pair.
5. The earth matrices added to every pair of conductors.
6. The bundle and Kron reductions.
7. `Y = s P⁻¹` from the reduced potential matrix.
8. The ideal transposition, when selected.

After the first call, the loop of a computation in `Float64` does not allocate. Each LU
decomposition writes its row interchanges into the pivots of the buffers, and the test of
the loop measures zero bytes.

### Buffers and tool levels

`initialize_buffers` is the one action that builds buffers. It builds the base keys and
folds over `formulation.formulas`, so that a formula adds its own arrays. A base key is a
value computed once and read again, or a value one step hands to the next. A value made and
used inside one step is inline. One formula's scratch is never a base key.

A name is defined at the level whose every member uses it, and nowhere higher. The engine
defines two spectral tools. `SpectralIntegral` is an integral of a kernel in λ, and its
constructor chooses the complex path and the points where the quadrature subdivides.
`integrate` evaluates it by a method selected with a `Val` tag, and maps the path and applies
its Jacobian. Every spectral earth-return formula uses them.

A family defines what every formula of the family declares. `EarthReturn` gives each formula
`Λ`, the mean of the direct and image fields, and a formula can override it. Everything else
is the formula's, in the bodies of its physics methods. `:martins2026` writes its kernels as
tagged variants of the formula and its pair terms in `to_buffer!`. The family has one private
operation, `i0m1`.

### Result

`LineParameters{T, F, N}` stores `Z` and `Y` as arrays indexed by row, column and frequency.
It also stores the frequencies, the label of each retained conductor, the line length, the
formulation and the trace. The trace `N` is typed per backend. For `LCMBackend` it is the
stage matrices `Zin`, `Pin`, `Zg` and `Pg` per frequency, or four empty arrays when the
computation did not keep them. A result without a trace stores an empty record, never
`nothing`.

### Numbers

A `Measurement` input keeps its type through the data model and the engine. `BigFloat` is
unsupported: SpecialFunctions has no Bessel methods for it, and the honest error is the
`MethodError` at the first call.

## Extensions

### Finite-element backend

`using Gmsh` loads `LineCableModelsGmshExt`. The extension reads the engine's problem, the
laws of `FEM` and the reduction plan of `reductions`. From the data model it reads the
placed regions, their cross-sections and the curves of their outlines. It defines
`defaults(::Type{GmshBackend})` and the `compute` method of `GmshBackend` for `FEM`.

Its files follow the computation. `options.jl` validates the record and defines
the table of native controls. `model.jl` builds the field model of a
problem: the material groups, the terminals and the medium of each cable. `geometry.jl`
builds the Gmsh geometry through `draw` and writes `physical.geo`. `mesh.jl` writes the
model data and meshes every frequency, one Gmsh process at a time per worker. `getdp.jl`
scans the frequencies with GetDP. `result.jl` reduces the
terminal matrices and builds the result. The twelve native sources under `getdp/` are the
GetDP problem definitions, kept byte for byte.

The backend refuses four things before it runs.

- A Γ vector whose length differs from the number of frequencies.
- A conductor of zero conductivity at the operating temperature.
- A cable tangent to the ground surface.
- A native path argument of 256 bytes or more.

Everything it writes goes under `.linecablemodels/fem/` of the working directory,
removed after the run unless `keep_run_directory` is set. The result is a
`LineParameters{Float64, FEM, N}`. For any number type of the problem, the backend computes
the nominal problem in `Float64` and returns no uncertainty.

The extension is the owner of the native data. The string comparisons that decode Gmsh's and
GetDP's outputs live in it and nowhere else.

### Uncertainty propagation

`using Measurements` loads `LineCableModelsMeasurementsExt`. The extension is methods only.
It defines `nominal` on a `Measurement` and four methods of the SpecialFunctions Bessel
functions on an uncertain complex argument, with their derivative identities. It defines no
name of its own. Everything else is the parametric `T` of the package. Γ and every option
refuse a `Measurement`, because a control is not a physical input. The test tree checks the
nominal parts against the deterministic engine and the sensitivities against finite
differences.

### Derivatives

`using ForwardDiff` loads `LineCableModelsForwardDiffExt`. The extension is methods only.
It defines `nominal` on a `Dual` and four methods of the SpecialFunctions Bessel functions on
a complex dual argument, with the derivative identities of the Measurements extension. It
defines no name of its own. A real dual argument takes ForwardDiff's own rules. The test tree
checks the partial derivatives of the Bessel functions and of the line parameters against
central differences.

### Rigorous enclosures

`using IntervalArithmetic, Arblib` loads `LineCableModelsIntervalArithmeticExt`. The extension
is methods only. It defines `nominal` on an `Interval`, its midpoint, and five methods of the
SpecialFunctions Bessel functions on an interval argument. Each method goes through the
conversions of IntervalArithmetic to and from Arb. It defines no name of its own. The test
tree checks each enclosure against the definition in Arb at 256 bits.

## Register and lexicon

`test/architecture/names.toml` is the one record of every name the package defines. It has
two halves.

The lexicon lists the sanctioned terms, globally and per module, each with one role: verb,
noun, qualifier, symbol or preposition. Each name splits into its terms in exactly one way.
A function is a verb, a verb and its object, a quantity, a qualifier, or a predicate, in at
most two terms. A field, keyword, constant or module is a noun phrase of at most two terms.
A type is a noun phrase of at most three terms after `Abstract`. Each rejected word has a
stated replacement: `check` is `validate`, `get` is the name of the quantity, `data` is the
name of what it stores.

The register lists, per module and per extension, every function, type and constant the
module defines, as public or private, each with its meaning. The families list their
identifiers with a one-line reference each. The lists `lifts` mark the functions that lift a
value to a type with a `Val`.

```toml
[register."Engine.EarthReturn".formulas]
martins2026 = "the circumferentially averaged earth return of air, earth and mixed systems ..."
```

A name enters the register with the maintainers' confirmation, in the change that introduces
it, and leaves with its definition. Everything else is a method of a registered function or
type, or of a function of Base, a standard library or a dependency.

## Guards

`test/architecture/` enforces the structure on every test run. Each guard blocks at zero.
A baseline, a ratchet or a count does not exist.

| Guard | What it checks |
| --- | --- |
| ownership | a method that extends a foreign function mentions a type of its own module |
| placement | a module's files lie under its home, and no core method lies under `ext/` |
| names | functions and types of different modules have different names, unless they descend from one abstract type of `Commons` |
| shadowing | a name the root exports means what Base or a standard library exports under it |
| explicit imports | each foreign name is imported from its owner, and no import is stale |
| direction | a module references only earlier modules, its ancestors and its descendants |
| validate returns its subject | every `validate` method returns its first argument or throws |
| reserved verbs | no function name starts with `check_`, `verify_`, `ensure_`, `require_` or `assert_` |
| storage vocabulary | no `work`, `scratch`, `storage` or `cache`, and no type ending in `Buffer`, `Storage`, `Scratch` or `Cache` |
| repeated access names | an access chain's base and its first field differ |
| symbol switches and probes | a comparison against a Symbol literal, `applicable` and `@eval` are refused, and `isa` appears only in its admitted forms |
| decisions by absence | `isnothing`, `ismissing`, `something`, a comparison with `nothing` or `missing`, and a `Nothing` argument or parameter are refused |
| file hygiene | no `tempdir`, `mktempdir`, `tempname` or `homedir` under `src/` or `ext/` |
| Commons admission | a `Commons` definition is grammar, user API, or shared by two owners |
| reserved vocabulary | the physical constants and the matrix reductions are used by their names, never redefined |
| literal fingerprints | the digits of the physical constants appear only in their home |
| clones of Commons bodies | no function body repeats a public `Commons` body |
| clones of every body | no function body of at least 90 tokens repeats three quarters of another body's shingles. The rule admits a pair of constructors of package types and a pair of methods of formulas or formulations |
| the method table over the curves | a generic with a method on one member of the closed `Curve` union has a method for every other member, or one on the union |
| narrow floats | no number type narrower than `Float64`, and no literal of one, under `src/`, `ext/`, `test/` or `docs/src/` |
| offset reads | only the methods that apply, resolve, validate or map an offset read its field `t` |
| distance across on a concrete primitive | no `dimension` method on a concrete primitive type: the distance across is one computation over the contour |
| distance across outside `dimension` and `section` | under `ext/`, no hand computation of the distance across, such as `ro - ri`, `min(w, h)`, `2r` or a bounding box |
| register | the package defines exactly the names of `names.toml` |
| local functions | no method is defined inside a function body, and no local is bound to an anonymous function |
| forwarders | no function is only a second name for one call |
| root module | the root module only includes, imports, exports and declares names public |
| boxed captures | no package method's lowered code contains a `Core.Box` |
| typed values | every field has a concrete type or a bounded parameter, and every constructor types its positional values |
| `Val` tags | a `Val` appears in three places only: a lift, a literal tag after a formula and a type parameter |
| formula registry | the identifiers with a `defaults` method, `FORMULAS` and the register agree, and each has a physics method |
| naming rule | every name is built from sanctioned terms, in the sanctioned form |

A set of planted probes, the negative controls, proves that each guard reports a violation.
Beside the guards, Aqua checks the project, JET and `@inferred` check the inference of
`compute`, and the loop's allocation is measured exactly. Under `test/references/`, ten reference cases hold
the physics to their rows of record, analytical and finite-element.
