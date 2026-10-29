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
    compute(backend, problem, formulation; kw...)

Run one computation of `problem` with `formulation`. A method of `backend` selects the
default backend. A backend passed first chooses another.

# Arguments

- `problem`: problem to compute.
- `formulation`: formulas and their controls.

# Keywords

- `kw...`: computation options of the selected backend.

# Returns

- The result of the selected backend.
"""
function compute end

"""
    defaults(owner)

Return the declared option keys and default values of an owner.

# Arguments

- `owner`: formula, formulation, or backend type.

# Returns

- A `NamedTuple` whose value types define the option types.
"""
function defaults end

"""
    select(default, value)

Convert an option to its declared type or lift a formula selection into its slot.

# Arguments

- `default`: default value or formula family that determines the selection type.
- `value`: supplied option value or formula selection. A selection is an identifier, an
  identifier paired with options, or a formula.

# Returns

- The converted option or selected formula.
"""
function select end

"""
$(TYPEDEF)

Represent the common supertype of formulations and formula families.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Represent the common supertype of backend singletons that execute a formulation.
"""
abstract type AbstractBackend end

include("options.jl")

end
