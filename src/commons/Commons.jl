"""
    LineCableModels.Commons

Define the grammar that the other modules extend with methods, and the option records of
formulations and backends.
"""
module Commons

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

export build, resolve, validate, compute, defaults, select
export AbstractFormulation, AbstractBackend, FormulationOptions, ComputationOptions

"""
    validate(subject, context...)

Reject an invalid input before its consumer uses it. Return `subject`, the first positional
argument, unchanged.

`subject` is user input or an input derived from it. The optional context says what the
subject is checked for: its consumer, such as the problem, formulation, physical model or
result being built, or, for the result of a law or record, the function that produced it.
Data that the check needs follows the context. A method can warn on a value that it can
compute but that is questionable, such as one outside a formula's range of validity.

# Errors

- Throws a native exception that identifies the invalid field, value, and required
  condition.
"""
function validate end

"""
    resolve(context, declaration)

Return the geometry that `declaration` takes in `context`. The outline of a part gives its
geometry to the layer or part declared around it, and a placement locates a cross-section in
its frame.
"""
function resolve end

"""
    build(Target, declarations...)

Construct one completed `Target` from its declarations, such as the design of a cable from
its declared parts.
"""
function build end

"""
    compute(problem, formulation; kw...)

Run one computation of `problem` by `formulation`. The backend that runs the formulation
defines the computation as a method of this function, and the keywords are the options of
that backend.
"""
function compute end

"""
    defaults(owner)

Return the option keys of `owner` with their default values, as a `NamedTuple`. The type of
each default is the type of its option. The owner is a formula, formulation or backend type.
"""
function defaults end

"""
    select(default, value)

Fit the selection `value` to the slot or option key whose default is `default`. A scalar
option takes the type of its default, and a slot lifts a formula identifier to the formula
it names.
"""
function select end

"""
    AbstractFormulation

Supertype of the formulations and of the formula families.
"""
abstract type AbstractFormulation end

"""
    AbstractBackend

Supertype of the singletons that run a formulation.
"""
abstract type AbstractBackend end

include("options.jl")

end
