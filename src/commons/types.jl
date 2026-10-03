"""
$(TYPEDEF)

Supertype for complete LineCableModels calculation inputs.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Supertype for scientific and higher-order calculation selections.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Supertype for completed LineCableModels calculation results.
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

Retain computation-owned inputs. The receiving calculation or backend owns
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

Bind a concrete formulation and semantic selectors to its native domain method.
Calling the binding passes the selection, selectors, and runtime arguments in
that order.

$(TYPEDFIELDS)
"""
struct FormulaMethod{S, F, A <: Tuple}
    "Selected formulation whose concrete type owns equation dispatch."
    selection::S
    "Native domain method accepting the selected formulation first."
    method::F
    "Semantic Val selectors inserted before runtime arguments."
    arguments::A
end

"""
$(TYPEDSIGNATURES)

Bind a selected formulation to a native method and optional `Val` selectors.
Throw `ArgumentError` if a selector is not a `Val` instance.
"""
function FormulaMethod(selection::S, method::F, arguments...) where {S, F}
    all(argument -> argument isa Val, arguments) || throw(ArgumentError(
        "FormulaMethod semantic selectors must be Val instances"))
    return FormulaMethod{S, F, typeof(arguments)}(selection, method, arguments)
end

@inline function (bound::FormulaMethod)(arguments...)
    return bound.method(bound.selection, bound.arguments..., arguments...)
end
