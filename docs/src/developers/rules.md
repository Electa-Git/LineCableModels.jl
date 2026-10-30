# Rules and patterns

This page states the principles behind LineCableModels, the rules that every change
follows, with the reason for each, and the patterns that extend the package within them. The guards of
`test/architecture/` enforce the rules mechanically. A guard failure blocks a change. There
is no baseline, and an exemption exists only in the rule's own text.

Every rule serves one intent. Dispatch comes first, and the vocabulary is narrow, with short
action verbs. Dispatch gives a function its subject, and the module namespace and the
dispatch type give its context. Types select behavior, and names are few with one meaning
each. Structure replaces machinery.

## Principles

These principles underlie the rules below. The vocabulary of the package shows them in a few
verbs and many types, with modules made of methods of verbs that already existed.

### Protocols are the interface

A type joins the package by implementing protocols that already exist, those of Julia, such as
`convert`, `promote`, `==`, `map`, `eltype`, `write`, `read`, `run`, `in` and `show`, and those of the
grammar, `validate`, `resolve`, `build`, `compute`, `defaults` and `select`. The set of
verbs is closed. The set of types is open. A new capability is a new type with methods of
existing verbs, placed with the owner of the type. The package has no interface of its own
beyond its grammar, and a reader who knows the grammar knows every module.

### A verb exists only for an act with no native name

Before a verb is proposed, its call is written with a function of Base or of the grammar. If
that call reads true, the verb does not exist. Starting a process is `run`, writing a file is
`write`, decoding one is `read`, converting a record is `map`, comparing is `==`. The package
defines verbs only for acts Julia has no word for. Its verbs flatten a design, draw a
contour, integrate a spectral integrand, mesh a cross-section and scan frequencies with a
solver. An act
made of two acts is two lines, never one name.

### Predictability is the measure of a module

A module has the right form when a reader can write the signature of each definition from
its purpose before reading it. The definitions follow the computation in order, each a
method of a known verb on an owned type. Variety of names is a cost, never a feature. A
module that reads as a tedious sequence of expected overloads is a module that is right.

### Declarations are passive, and few verbs act on them

The nouns of the package are types with typed fields and no behavior: materials, layers,
parts, designs, systems, problems, formulations, records and results. The verbs act on them.
No declaration computes, caches or decides. A choice that selects behavior is a type
parameter set at construction, because dispatch reads it. A name that only identifies, such
as a terminal name or a region tag, stays a value, because only a table reads it.

### The structure is designed for the domain and implemented for the subset released

The data model is designed for every cable construction the field knows, and the engine for
every formula of its families, before either is implemented for the cases a release covers.
A later feature adds a type and its methods inside that structure. It never adds a second
structure, a flag or a mode, and never re-implements a concept for one case.

### Capability is refused where it is missing, at the call

The structure above a formula never rejects a physically valid input for a formula's sake.
The data model encodes, describes and draws equivalences. A formula that cannot treat an
input throws when it is called. A backend refuses before it runs what it cannot honor. The
refusal states the thing that is missing.

### Each component treats its own limitations

A shared layer never bends to the limitation of one consumer. A mesh that cannot resolve
contact between cables meshes the contact region itself. A formula that needs a cable in one
medium checks that itself. A check in a shared layer exists for a physical or structural fact
and for nothing else.

### Operators are named for their action, elements for what they are

A physical element takes a physical noun: `Circle`, `Annulus`, `CircularPattern`,
`Helix`. An operator on an element takes the name of its action: `Offset`, `Repeat`, `Lay`,
`Enclosure`. A selector symbol stands for a definite act, law or source: `:discover`,
`:download`, `:unified`, `:schelkunoff1934`. No selector is `:auto` or a vague quality.
`:default` selects no formula: it routes to the slot's default formula at construction.

### Equality is equality of the computation

Selections are equal when they produce one computation. A formula is its identifier and its
options, records compare by value, and equal formulas have one type. The same formula
object means the same computation on every backend. A result records the formulation that
computed it, so that a result never claims a computation that did not take place.

### Selection happens once, evaluation reads types

Each choice is made at construction. A Symbol becomes a type by its lift, a keyword fills a
typed record, and a formula fills a typed slot. `:default` resolves to the slot's default
formula, which then fills the slot. The computation that follows reads types and never
decides. What a reader sees inside a loop is physics and arithmetic.

### Reference results are frozen with their origin, and deviations are explained

A reference is a file with the commit, the script and the run that produced it, and a
computation is held to it at a stated level. Each deviation is recorded with its cause as the
known status of that case. Tolerances are not chosen from the difference that was achieved,
and a failing case stays in the set.

### A measurement precedes a design decision

A choice between structures rests on a number, such as an allocation count, a timing, an
inference report or a difference against a reference. A passing check, an existing test or a page of
documentation is evidence of what it encodes, never the reason a structure stays. Physics and
structure are the arguments.

### Generality is bounded by the mathematics

A generalization with a physical meaning is implemented when its mathematics is within
reach, never deferred because no present case needs it, whereas a convenience keeps the
scope it was given and grows into nothing.

### The notation of the physics is the vocabulary

The symbols of the equations are terms of the lexicon: `rho`, `eps_r`, `mu_r`, `r_in`,
`r_ex`, `Γ`, `s`. The code writes the equation with the symbols the docstring defines, once,
inline, so that the physics is read where it is computed.

## Rules

### One owner and one name per concept, with one meaning per name

Each concept has one owner module, one name and one implementation. Each name has one
meaning. Functions owned by different modules have different names. So do types, unless each
descends from one abstract type that `Commons` owns, as the eight `Formula` types do. A name
the root exports means the same object as the name that Base exports under it.

Reason: readers lose the program where a concept has a second name or a name stands for two
concepts. The naming guards check this on every run.

### The register lists each name

The package defines only the names that `test/architecture/names.toml` lists. Everything
else is a method of a registered function or type, or of a function of Base, a standard
library or a dependency. A name enters the register with the maintainers' confirmation, in
the change that introduces it, with its meaning. It leaves with its definition.

Reason: a private helper exists because naming a fragment costs its writer less than finding
the fragment's owner. Every threshold on the form of a helper is a door, and the cheapest
form that passes the threshold is the one that gets written. So the gate judges no name. It
brings every name to the maintainers for confirmation, and no other path can define behavior
under a new name.

A helper takes these forms.

- A forwarder or an alias.
- A reader or a repacker of an owned value.
- A wrapper around Base or a dependency.
- One function per step of one fixed algorithm.
- A normalizer between two owned representations.
- A shared function extracted from look-alike code.
- A function barrier around lost inference.
- A private dispatcher beside the dispatch structure.
- A sentinel for an absent value.

Its disguises are a missing underscore, a move into `Commons`, nesting inside the caller, a
field-less type used as a namespace, a new label such as stage, kernel or adapter, and a
test that pins it.

These are not helpers.

- A method of a Base or dependency function on an owned type.
- A method of a grammar verb, or a formula method.
- A domain operation with its own name.
- The loop of a computation as its own method.

Before proposing a name, look for the effect in three places. The first is a method of an
existing generic. The second is a Base function overloaded on the owned type. The third is
inline code at the one place that needs it. A proposal states the name, its meaning, and why
none of the three applies.

### No function is a second name for one call

A function whose every method passes its own arguments, unchanged, to one call does not
exist in this package. No method is defined inside a function body, a `do` block, an
anonymous function, a `let` or a macro. A local is never bound to an anonymous function.
Anonymous functions appear only as arguments.

### Types select behavior

The code does not use `applicable` or `@eval`, and it compares nothing against a `Symbol`
literal. A type test appears only in a `catch` clause, or as the condition of a selection
from a collection. Such a selection is a comprehension filter, or the predicate of `filter`,
`findall`, `findfirst` or `findlast`. A choice known at construction is lifted once to a
type or a type parameter and dispatched after that.

```julia
Material(:conductor, rho, eps_r, mu_r, T0, alpha)   # lifts :conductor to Material{T, :conductor}
select(EarthReturn.Formula, :unified)                # lifts :unified to Formula{:unified}
```

Reason: a Symbol comparison inside a computation is a hidden switch that dispatch cannot
see, and each one hides a default. The one admitted home of a Symbol or string
comparison is the module that owns external data. The Gmsh extension decodes the native
outputs of its own solver runs, and no package module decodes external data.

### No value stands for absence

The code calls none of `isnothing`, `ismissing` and `something`. It compares nothing with
`nothing` or `missing`, and no argument or type parameter is `Nothing` or `Missing`. The one
exemption is the result of a Julia function that returns `nothing` by its own protocol, such
as `findfirst`, `match`, `tryparse` or `Sys.which`, tested directly.

A fact that does not hold is a row that does not exist. A result without a trace stores an
empty record. An optional external resource takes a positive default: a public type that
names where the resource comes from, or a path.

```julia
getdp_executable = GetdpBinary{:discover}()   # the default: PATH, else the download
getdp_executable = :download                  # the pinned artifact
getdp_executable = "/opt/getdp/bin/getdp"     # a path, used as is
```

Reason: a decision by absence is a hidden branch with a hidden default, and the type system
cannot see either.

### Typed fields and typed constructors

Each field of a struct has a concrete type, such as `Vector{T}`, or a type parameter with a
bound other than `Any`. A `Union` of concrete types may be a parameter of a field's type,
never a field's own type. Each constructor types each positional value with a type or a
bound other than `Any`. A struct without type parameters declares an inner constructor,
because Julia's converting default constructor would accept anything.

### Input checks use `validate`

`validate(subject, context...)` is the one checking verb. No function name starts with
`validate_`, `check_`, `require_`, `assert_`, `verify_` or `ensure_`. A `validate` method
checks user input, or an input derived from it, before a consumer uses it. It throws on
invalid input with a native exception that states the field, the value and the condition. It
may warn on input that is computable but questionable. It returns its subject unchanged.

A check runs where its fact is built: the inner constructor of a record calls `validate` on
the new object. The context is the consumer, such as the problem, the formulation, the
physical model or the result being built, and never a plumbing function or a `Val` tag.

```julia
function EarthLayer{T}(rho::T, eps_r::T, mu_r::T, thickness::T) where {T <: Real}
    return validate(new{T}(rho, eps_r, mu_r, thickness))
end
```

### Options

Options follow two levels with one grammar, `compute(problem, formulation; kw...)`.
Formulation options are the controls a formulation determines or needs, in
`FormulationOptions{O}`, owned by a formula or a formulation type. Computation options
control the backend, whatever the formulation, in `ComputationOptions{B}`, owned by a
backend. Each record is built once, from `defaults(O)` and the keywords given. An unknown
key throws. Each value takes its declared type. `validate` checks the values. A record is
never rebuilt in a loop.

Keys are flat per owner. A selectable sub-method is one section keyed by its selector, such
as `integral = (method = :quad, rtol = 1e-8)`. A Symbol choice is lifted at construction.
Neither `Nothing` nor a function is stored in a record.

### A `Val` has three homes

A `Val` has no behavior of its own. It appears in a lift, in a function the register marks
as a lift, such as `Formula(::Val{ID}; kw...)` or `Material(::Val{K}, ...)`. It appears as a
literal tag after a formula object, such as `integrate(integral, Val(:quad), ...)`. It
appears as a type parameter set at construction. A method whose first positional argument is
a `Val` is a constructor or a registered lift. A `Val(x)` whose `x` is not a literal appears
only in a constructor, a lift or an options constructor.

Reason: a `Val` with a literal is free. A `Val` with a runtime value is a dispatch at run
time, harmless once per computation and harmful per frequency. With these homes, every lift
happens at construction.

### Buffers, inference and the loop

`initialize_buffers` is the one action that builds buffers. It builds every buffer once,
before the loop, by folding over the formulation's formulas. A formula that needs arrays of
its own adds its method. Storage types and storage-building functions for one formula do not
exist. A named `*Buffers` type exists only in `Commons`, for a component that several owners
share. No function name ends in `_workspace` or `_storage`. No field, key, parameter or
local is named `work`, `scratch`, `storage` or `cache`.

When inference is lost, find where the type is lost, such as a boxed captured variable, a
state initialized as `nothing` or a type discovered inside a loop. Fix the structure: make
the decision first and allocate before the loop, which is its own method. Never wrap the
symptom in a helper. No package method's lowered code contains a `Core.Box`.

### Levels of tools

Each name is defined at the level whose every member uses it, and nowhere higher. An engine
tool is one every formula of its concern uses, delivered through the engine's buffers and a
docstring that states its purpose. Family tools are those that every formula of the family
declares, a no-op default included. Everything else is the formula's, in the bodies of its
physics methods. Inline mathematics is preferred to named fragments, so that the physics is
not concealed by the code that writes it. A term goes to a buffer when it is computed once
and reused.

The tools grow from the engine itself. A name at the level of the engine or of a family is
declared by every formula of that level.

### `compute` computes

`compute` runs as its controls say and reports a requested target it did not meet. It leaves
error budgets, refinement loops, audits and resource guards out. Scientific validation
belongs to the tests and to the user's comparisons with finite-element results.

### Numbers

Parametric types where they matter, type stability, and no over-constrained `<:Real`. A
`Float32` or mixed-type input must not crash. Precision retention is not a goal: the problem
computes in `Float64`.

### A native name or none

An act that Julia already names is written with Julia's name, inline. Starting an external
process is Base's `run` on a pipeline. Walking a record to convert its numbers is a method
of Base's `map`. Writing a file is a method of Base's `write`. No package verb bundles the
preparation of inputs with the start of a process, because that verb would mean something
different at every use.

### Files

Anything the package writes goes to the path the user gives, else under the untracked
`.linecablemodels/<kind>` folder of the current directory. The package never writes to the
home folder, to `tempdir()` or to the package source tree.

### Names

Each name is built from terms sanctioned in `names.toml`, globally or for its module. Each
term has one role. Compounds such as `line_length` or `earth_return` are one term each. A
function is a verb, a verb and its object, a quantity, a qualifier, or `is` or `has` and a
term, in at most two terms. A field, keyword, constant or module is a noun phrase of at most
two terms. A type is a noun phrase of at most three terms after `Abstract`. A preposition is
a term only as a whole field or keyword name, such as `at`.

Dispatch gives the subject: a function's name and a package type its methods take do not
share a term. Namespaces give the context: a function reached through its module does not
share a word with the module's path. A global verb is defined only in `Commons`, a module's
verb only in that module. No function of the package takes a name that Base or Core exports.
Methods of Base and other packages are free.

Reason: a vocabulary that repeats the subject or the context inside names grows without
bound, and one concept soon has several names. A narrow vocabulary with dispatch as the
subject keeps one name per concept.

### Tests

Every test calls feature code and checks a feature result. A test asserts only what the code
under test can produce. It leaves out what a constructor or `validate` already blocks, and
defensive branches for unrealistic failures stay out of it. Never loosen a tolerance, delete
or weaken an assertion, skip a test or loosen a guard. A `Float32` input is tested to build
and to compute with inferred types, never for precision. Coverage is reached by tests that
check feature results, not by tests written to touch lines.

### Documentation

Docstrings follow the [conventions](conventions.md). Every physical argument, return value,
field and constant states its SI unit. A docstring describes implemented behavior and states
each fact once. Vale and cspell pass on every maintained file. A finding is fixed by
rewriting the sentence, never by a directive, and a word joins the dictionary only by a
decision of the maintainers.

## Patterns

Each pattern below is the complete list of what a change adds. Anything more is a change to
the architecture and is discussed with the maintainers first.

### Add a formula to a family

Three parts. The example is an earth-return formula that gives the impedance alone.

```julia
# in src/engine/earthreturn/EarthReturn.jl
const FORMULAS = (:unified, :carson1926)

# in src/engine/earthreturn/carson1926.jl
defaults(::Type{Formula{:carson1926}}) = (;)

function impedance!(formula::Formula{:carson1926}, Zg::AbstractMatrix, buffers, input,
        media, s, k)
    # the closed form, written here, from the per-frequency arguments
    return Zg
end
```

The registry line in `names.toml` gives the identifier with a one-line reference, under
`[register."Engine.EarthReturn".formulas]`. Without potential coefficients, this formula is
refused by the formulation in the `earth_admittance` slot, by dispatch. A formula with its
own controls declares them in `defaults` and checks them in `validate(options)`. A formula
with shared state between its quantities defines `prepare!` and reads its state in
`impedance!` and `potential!`. A formula that needs arrays adds
`initialize_buffers(formula, problem, input, buffers)` and returns the buffers merged with
its own, keyed by its identifier. Its equations go to the "Formulas" section of the family
module's docstring, as the [conventions](conventions.md) say.

### Add a slot to the formulation

A new slot is a field with a default in `CoaxialFormulation`, typed by its family, and a
keyword with its default identifier in the keyword constructor. The family is a child module
of `Engine` with `Formula`, `FORMULAS` and the lift. The function the family's methods
extend is declared in `Engine`, with its docstring, and the loop calls it at the step where
its quantity enters.

### Add a backend

A backend is a singleton of `AbstractBackend` in `Engine`. It adds four things.

- `defaults(::Type{B})`, its computation options.
- `compute(::B, problem, formulation; kw...)`, which builds the options once.
- `compute(::B, problem, formulation, options::ComputationOptions{B})`, which computes.
- `validate(formulation, ::B)` per formulation it implements, with a throwing fallback.

A backend implemented in an extension defines its methods there and leaves its singleton in
`Engine`, as `GmshBackend` does. A formulation's default backend is one `backend` method.

### Add a cross-section

A cross-section is one type under `AbstractPrimitive` with the methods of the unary
interface. The methods are `area`, `centroid`, `support`, `distance` and `boundary`,
`resolve` against an outline and in a frame, and `r_in` and `r_ex` when the cross-section is
circular. The `CAD` module adds its `Contour`. The Gmsh extension adds its native section.
Nothing in the parts, the design or the engine changes.

### Add a part or a builder

A part is a node type under `AbstractCablePart` with one job and one `resolve` method.
Builders are functions that lower to existing nodes and do not add a capability. Both enter
the register with their meanings.

### Add an extension

An extension declares its weak dependency in `Project.toml` and imports each name it uses
from its owner. It defines methods of registered functions and of Base. A name of its own
enters a `[register.<Extension>]` table of `names.toml`, with its terms sanctioned in the
extension's own section. The extension never defines a core method, and the rules apply to
it as to a package module.

### Add a name

A new name is proposed with its meaning and its terms. The proposal says why the effect is
not inline code, a method of a Base function, or a method of a registered name. The
maintainers confirm it. The register and the lexicon are updated in the same change as the
definition.
