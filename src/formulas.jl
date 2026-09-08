"Return a short scientific description of a registered formulation."
function description end

"""
$(TYPEDEF)

Store one declarative formula selection until its owning formulation resolves
the identifier and overrides into a concrete formula type.

`FormulaDefinition` is produced by [`formula`](@ref). It does not participate in a
numerical loop.

$(TYPEDFIELDS)
"""
struct FormulaDefinition{ID, Order, P <: NamedTuple, H <: NamedTuple, O <: NamedTuple, E}
    "Explicit formula parameters, without evaluated physical state."
    parameters::P
    "Explicit callable overrides, without numerical workspaces."
    hooks::H
    "Explicit numerical sections owned by the consuming equation."
    options::O
    "Optional equivalent homogeneous-earth selection owned by this formula."
    equivalent_earth::E
end

"""
$(TYPEDEF)

Bind one formula identity and optional semantic selectors to a domain method.

Calling the bound method inserts `Val(ID)` before the stored selectors and
runtime arguments. Formula catalogues use this invariant to retain owner-local
dispatch while carrying the selected formula identity as concrete type
information.

$(TYPEDFIELDS)
"""
struct FormulaMethod{ID, F, A <: Tuple}
    "Owner-local domain method selected by the formula."
    method::F
    "Semantic `Val` selectors inserted before runtime arguments."
    arguments::A
end

"""
$(TYPEDSIGNATURES)

Bind a formula identity and optional semantic selectors to a domain method.

# Arguments

- `identifier`: Formula identity carried as `Val{:ID}`.
- `method`: Owner-local domain method whose first argument accepts that
  identity.
- `arguments`: Optional semantic selectors inserted before runtime arguments.

# Returns

- A callable [`FormulaMethod`](@ref).

# Errors

- Throws `ArgumentError` when a stored semantic selector is not a `Val`.
"""
function FormulaMethod(::Val{ID}, method::F, arguments...) where {ID, F}
    all(argument->argument isa Val, arguments) || throw(ArgumentError(
        "FormulaMethod semantic selectors must be Val instances"
    ))
    return FormulaMethod{ID, F, typeof(arguments)}(method, arguments)
end

@inline function (bound::FormulaMethod{ID})(arguments...) where {ID}
    return bound.method(Val(ID), bound.arguments..., arguments...)
end

"""
$(TYPEDSIGNATURES)

Select a registered formula without exposing its owner module or concrete
wrapper type. The receiving formulation determines the formula family from the
keyword slot in which the selection appears.

# Arguments

- `identifier`: Stable formula identifier.
  `:default` requests the applicable choice from the resolved problem, geometry,
  earth characteristics and backend. It is not a fallback after a failed formula.
  Cable-insulation and semicon-admittance defaults explicitly select lossless
  dielectric relations. Unsupported contexts fail before frequency evaluation.

# Keywords

- `order`: Position of an equivalent homogeneous-earth reduction relative to
  material frequency dependence. `:before` applies EquivalentHomogeneous before FrequencyDependent, `:after`
  applies EquivalentHomogeneous after FrequencyDependent, and `:default` selects the receiving formulation's
  default. Non-EquivalentHomogeneous formula slots accept only `:default`.
- `parameters=(;)`: Explicit model parameters accepted by the owning formula.
- `hooks=(;)`: Callable overrides at the owning formula's documented variation points.
- `options=(;)`: Numerical operation sections, such as `integration=(method=:quad, options=(;))`.
- `equivalent_earth=nothing`: Explicit reduction for a compatible external formula.

# Returns

- A concrete declarative selection resolved before computation.

# Examples

```julia
earth = formula(:Carson1926)
soil = formula(:default)
equivalent = formula(:default; order=:before)
```
"""
function formula(identifier::Symbol; order::Symbol = :default,
        parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;), equivalent_earth = nothing)
    order in (:default, :before, :after) || throw(ArgumentError(
        "formula order must be :default, :before, or :after"
    ))
    return FormulaDefinition{identifier, order, typeof(parameters), typeof(hooks),
        typeof(options), typeof(equivalent_earth)}(
        parameters, hooks, options, equivalent_earth)
end

"Return the stable formula identifier of a formula value."
function formula_id end

formula_id(::FormulaDefinition{ID}) where {ID} = ID
