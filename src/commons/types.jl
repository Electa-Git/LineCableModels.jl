"""
$(TYPEDEF)

Supertype for complete LineCableModels computation inputs.
"""
abstract type AbstractProblemDefinition end

"""
$(TYPEDEF)

Supertype of the core formulations: the scientific and numerical methods selected to compute
one problem.
"""
abstract type AbstractFormulation end

"""
$(TYPEDEF)

Supertype of the studies that compute every point of a problem space with an `inner`
formulation: [`Combinatorial`](@ref), `LinearError` and `MonteCarlo`. The inner formulation
is a core formulation, an `AbstractFormulation`, and a study is a formulation of a problem
space, not of a problem. Each subtype has the field `inner`. A study whose formulation
determines or needs controls of its own, such as the trials of `MonteCarlo`, stores them in
`options`. They are scalars, and `formulation_options` of its type normalizes them. The
computation options of the run, `verbosity`, `on_result` and `retain_details`, are the
`options` of `compute(problem::ParametricProblem, study; options)`, which
`computation_options(T, options)` normalizes once, one method for every study.
"""
abstract type AbstractTraversal end

"""
$(TYPEDEF)

Supertype of the spaces of core problems that a study computes:
[`ParametricProblem`](@ref). A problem space stores a core problem, an
`AbstractProblemDefinition`, or a space of them, and the computation options of the backend
for every point.
"""
abstract type AbstractProblemSpace end

"""
$(TYPEDEF)

Supertype of the builders of the user-facing language, such as `insulation`, `wires` and
`formula`. Each builder is a singleton that returns a `Target`.

Calling a builder lifts it over all its arguments and keywords, including those inside the
tuples, vectors, named tuples and pairs it receives. With no `Grid` or `Gridspace` among
them, the call returns `build(builder, arguments...; keywords...)`. With one, it returns
the `Gridspace{Target}` whose points do, composed by `combine = :product` or `:zip`. The
points of a product vary the first source fastest, in the order the call writes them.
The owner of a builder, or a user for an own type, adds `build` methods, which receive
scalars only. A call method of a builder's own type would take its calls before the
grammar, so a `Grid` keyword would reach it unlifted: extend `build` instead.
"""
abstract type Builder{Target} end

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
validation. Construction copies the supplied named tuple once, with `deepcopy`, which makes
the record a value that an edit of an array the caller passed does not reach. Uncertain
values keep their correlations.

$(TYPEDFIELDS)
"""
struct FormulationOptions{NT <: NamedTuple}
    "Supplied formulation options."
    data::NT
    # The parameter is the type of the data, so no record holds an abstract parameter.
    function FormulationOptions{NT}(data::NamedTuple{K, V}) where {K, V, NT <: NamedTuple{K, V}}
        return new{NT}(data)
    end
end

"""
$(TYPEDSIGNATURES)

Construct formulation inputs from a named tuple. A `Grid` or `Gridspace` at any depth of
`data` returns a `Gridspace{FormulationOptions}`. Owner-specific interpretation and
validation occur when the formulation resolves the options, not when this record is
constructed.
"""
function FormulationOptions(data::NamedTuple)
    caller = data -> (copied = deepcopy(data); FormulationOptions{typeof(copied)}(copied))
    return parameterize(FormulationOptions, caller, (data,))
end

"""
$(TYPEDSIGNATURES)

Construct formulation inputs from keywords, as from the named tuple of the keywords, or an
empty record with no keywords.
"""
FormulationOptions(; kwargs...) = FormulationOptions((; kwargs...))
# A record passes through its own constructor. A record of another kind has no method.
FormulationOptions(options::FormulationOptions) = options

"""
$(TYPEDEF)

Retain computation-owned inputs. The receiving computation or backend owns
defaults and validation.

$(TYPEDFIELDS)
"""
struct ComputationOptions{NT <: NamedTuple}
    "Supplied computation options."
    data::NT
    ComputationOptions(data::NamedTuple) = validate(new{typeof(data)}(data))
end

# Computation options control one computation. They do not vary, so no Grid or Gridspace
# enters them.
function validate(options::ComputationOptions)
    has_grid(options.data) && throw(ArgumentError(
        "computation options do not vary; vary the problem or the formulation instead"
    ))
    return options
end

"""
$(TYPEDSIGNATURES)

Construct computation inputs from keywords, or an empty record with no keywords.
The receiving computation validates the supported keys and values.
"""
ComputationOptions(; kwargs...) = ComputationOptions((; kwargs...))
# A record passes through its own constructor. A record of another kind has no method.
ComputationOptions(options::ComputationOptions) = options

"""
$(TYPEDEF)

Retain computation-owned supplemental output. The producing computation defines its contents. Immutability is shallow: arrays and other mutable payload values
are neither copied nor frozen. Every record has this one type: its contents change the keys
of `data`, never the type of the record or of the result that holds it.

$(TYPEDFIELDS)
"""
struct ComputationDetails
    "Named-tuple supplemental output, accessed explicitly through `.data`."
    data::NamedTuple
    ComputationDetails(data::NamedTuple) = new(data)
end

"""
$(TYPEDSIGNATURES)

Construct supplemental output from keywords. With no keywords, return an empty
record indicating that the computation supplied no supplemental output.
"""
ComputationDetails(; kwargs...) = ComputationDetails((; kwargs...))
# A record passes through its own constructor. A record of another kind has no method.
ComputationDetails(details::ComputationDetails) = details

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
struct Expression{S, F, A <: Tuple{Vararg{Val}}}
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
"""
function Expression(selection::S, method::F, arguments::Val...) where {S, F}
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
