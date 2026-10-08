"""
$(TYPEDEF)

Supertype for complete LineCableModels computation inputs.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Supertype for scientific and higher-order computation selections.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Supertype for completed LineCableModels computation results.
"""
abstract type AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for direct results owned by LineCableModels computations.
"""
abstract type AbstractCoreResult <: AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for completed finite collections whose element type is `T`.
"""
abstract type AbstractResultSpace{T} <: AbstractProblemResult end

"""
$(TYPEDEF)

Supertype for deterministic result spaces whose element type is `T`.
"""
abstract type AbstractParametricResult{T} <: AbstractResultSpace{T} end

"""
$(TYPEDEF)

Supertype for uncertainty result spaces whose element type is `T`.
"""
abstract type AbstractUncertaintyResult{T} <: AbstractResultSpace{T} end

"""
$(TYPEDEF)

Retain formulation-owned inputs. The selected formulation defines defaults and
validation. Construction preserves the supplied named tuple, including the
identity of mutable values within it.

$(TYPEDFIELDS)
"""
struct FormulationOptions{NT <: NamedTuple}
    "Supplied formulation options."
    data::NT
    FormulationOptions(data::NamedTuple) = new{typeof(data)}(data)
end

"""
$(TYPEDSIGNATURES)

Construct formulation inputs from keywords, or an empty record with no keywords.
Owner-specific interpretation and validation occur when the formulation resolves
the options, not when this record is constructed.
"""
FormulationOptions(; kwargs...) = FormulationOptions((; kwargs...))

"""
$(TYPEDEF)

Retain computation-owned inputs. The receiving computation or backend owns
defaults and validation.

$(TYPEDFIELDS)
"""
struct ComputationOptions{NT <: NamedTuple}
    "Supplied computation options."
    data::NT
    ComputationOptions(data::NamedTuple) = new{typeof(data)}(data)
end

"""
$(TYPEDSIGNATURES)

Construct computation inputs from keywords, or an empty record with no keywords.
The receiving computation validates the supported keys and values.
"""
ComputationOptions(; kwargs...) = ComputationOptions((; kwargs...))

"""
$(TYPEDEF)

Retain computation-owned supplemental output. The producing computation defines its contents. Immutability is shallow: arrays and other mutable payload values
are neither copied nor frozen. Nested diagnostic type bounds are preserved.

$(TYPEDFIELDS)
"""
struct ComputationDetails{NT <: NamedTuple}
    "Named-tuple supplemental output, accessed explicitly through `.data`."
    data::NT
    ComputationDetails(data::NamedTuple) = new{typeof(data)}(data)
end

"""
$(TYPEDSIGNATURES)

Construct supplemental output from keywords. With no keywords, return an empty
record indicating that the computation supplied no supplemental output.
"""
ComputationDetails(; kwargs...) = ComputationDetails((; kwargs...))

"""
$(TYPEDEF)

Store one passive formula selection until its defining formulation resolves the
identifier, model parameters, and formulation-owned physical and numerical options.

$(TYPEDFIELDS)
"""
struct FormulaDefinition{ID, Order, P <: NamedTuple, O <: FormulationOptions, E}
    "Explicit model parameters without evaluated physical state."
    parameters::P
    "Explicit physical choices and numerical controls owned by the selected equation."
    options::O
    "Optional equivalent homogeneous-earth selection."
    equivalent_earth::E
end

"""
$(TYPEDEF)

Hold the formula object, the operation and the `Val` selectors of one expression of a
formula. Calling the expression passes the formula, the selectors and the runtime arguments
to the operation, in that order.

$(TYPEDFIELDS)
"""
struct Expression{S, F, A <: Tuple}
    "Selected formula whose concrete type selects the operation's method."
    selection::S
    "Native domain method accepting the selected formulation first."
    method::F
    "Semantic Val selectors inserted before runtime arguments."
    arguments::A
end

"""
$(TYPEDSIGNATURES)

Build the expression of `selection` for the operation `method` and optional `Val` selectors.
Throw `ArgumentError` if a selector is not a `Val` instance.
"""
function Expression(selection::S, method::F, arguments...) where {S, F}
    all(argument -> argument isa Val, arguments) || throw(ArgumentError(
        "Expression semantic selectors must be Val instances"))
    return Expression{S, F, typeof(arguments)}(selection, method, arguments)
end

@inline function (expression::Expression)(arguments...)
    return expression.method(expression.selection, expression.arguments..., arguments...)
end

"""
$(TYPEDEF)

Hold a formula resolved at one evaluation point: the formula, the input of that point and the
values that the formula's expressions share there. `expression(functor, workspace)` evaluates
an expression of the formula at that point. The state holds plain values. Arrays come from
the buffers of the workspace.

$(TYPEDFIELDS)
"""
struct Functor{F, I <: NamedTuple, S <: NamedTuple}
    "Formula whose expressions this Functor evaluates."
    formula::F
    "Inputs of the evaluation point, with the options of the expression evaluated there."
    input::I
    "Plain values that the expressions of the formula share at that point."
    state::S
end

"""
$(TYPEDSIGNATURES)

Build the Functor of `formula` at the point that `input` describes. A formula without a method
of its own does not share values, and its state is empty. A family or a formula that checks its
input, shares values or reads arrays from the buffers of `workspace` adds a method on its own
type.
"""
Functor(formula, input::NamedTuple; workspace = nothing) = Functor(formula, input, (;))

"""
$(TYPEDSIGNATURES)

Return the Functor of the same formula and state at a point whose input extends
`functor.input` with `extension`, such as one conductor pair of an earth calculation.
"""
function Functor(functor::Functor, extension::NamedTuple)
    return Functor(functor.formula, merge(functor.input, extension), functor.state)
end
