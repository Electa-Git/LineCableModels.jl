"""
Return owned scientific text for a registered formula or formulation.
`compact=true` selects its short display name; `formula_id` remains its identity.
"""
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
runtime arguments. Formula catalogs use this invariant to retain owner-local
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
  Cable-insulation and semicon-admittance `:default` selections route to the
  explicit `:lossless` dielectric relation. Unsupported contexts fail before
  frequency evaluation.

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
earth = formula(:carson1926)
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

"""
Return the stable formula identifier of a formula value.
"""
function formula_id end

"Return the stable formula identifier carried by a bound formula hook."
formula_id(::FormulaMethod{ID}) where {ID} = ID

formula_id(::FormulaDefinition{ID}) where {ID} = ID

"""Expose a requested formula identifier and its explicit parameters and overrides."""
function Base.NamedTuple(value::FormulaDefinition{ID,Order}) where {ID,Order}
    return (identifier=ID, order=Order, parameters=value.parameters, hooks=value.hooks,
        options=value.options, equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end

"""
$(TYPEDSIGNATURES)

Describe a retained selection through its owning type. This does not construct
a formula, evaluate overrides, or substitute current defaults.
"""
description(source::Pair{<:Type,<:NamedTuple}; compact::Bool=false) = description(first(source); compact)
formula_id(source::Pair{<:Type,<:NamedTuple}) = formula_id(first(source))

description(::Missing; compact::Bool=false) = "method unavailable"
description(::Nothing; compact::Bool=false) = "none"
formula_id(::Missing) = missing
formula_id(::Nothing) = :none

"""
$(TYPEDSIGNATURES)

Label ordered formulation selections without interpreting their representation.

# Arguments

- `sources`: Native formulations or owner-bound retained declarations.

# Keywords

- `roles`: One `:reference`, `:candidate`, or `:none` role per source.
- `quantity`: Physical quantity selected through the observation grammar, or
  `nothing` for all equation choices.
- `compact=true`: Use the short owned names for both root and child selections.

# Returns

- One text label per source, in input order, without candidate numbering.
  Consumers select unique quantity-relevant formulations before presentation.

# Notes

Owners expose native children through `pairs(source; quantity)` and retained
children through `pairs(owner, record; quantity)`, scoped
identity through `formula_id`, controls through `formulation_options`, and
text through `description`. This formatter only decides common-field omission
and reference prefixes. Unsupported owner methods are not caught.
"""
function description(sources::AbstractVector;
        roles=fill(:candidate,length(sources)),
        quantity=nothing, compact::Bool=true)
    length(roles)==length(sources) ||
        throw(DimensionMismatch("one role is required per formulation"))
    all(in((:reference,:candidate,:none)),roles) ||
        throw(ArgumentError("description roles must be reference, candidate, or none"))
    isempty(sources) && return String[]
    selections=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...; quantity)) for source in sources]
    complete=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...)) for source in sources]
    candidates=findall(!=(:reference),roles)
    # Scientific comparisons use owner-scoped identifiers, never display text.
    scopes=unique([scope for index in candidates for (scope,_) in complete[index]
        if !isempty(last(scope))])
    varying=filter(scopes) do scope
        values=[[(formula_id(value),formulation_options(value)) for (key,value) in complete[index]
            if key==scope] for index in candidates]
        !all(value -> isequal(value,first(values)),values)
    end
    identities=[[scope for (scope,_) in entries if isempty(last(scope))] for entries in complete]
    candidate_owners=unique(first(ids) for ids in identities[candidates] if !isempty(ids))
    inner_owners=unique(last(ids) for ids in identities if length(ids)>1)
    return map(eachindex(sources)) do index
        prefix=roles[index]===:reference ? "Reference" : ""
        parts=String[]
        for (scope,value) in selections[index]
            if isempty(last(scope))
                position=findfirst(==(scope),identities[index])
                show_identity=position==1 ?
                    (roles[index]!==:candidate || length(identities[index])>1 || length(candidate_owners)>1) :
                    length(inner_owners)>1
                show_identity && push!(parts,description(value;compact))
                controls=formulation_options(value)
                peer=[other for entries in selections for (key,other) in entries if key==scope]
                if !isempty(controls) && any(other -> !isequal(formulation_options(other),controls),peer)
                    push!(parts,description(scope,value;compact,settings=true))
                end
            else
                peer=[other for entries in selections for (key,other) in entries if key==scope]
                changed=any(other -> !isequal((formula_id(other),formulation_options(other)),
                    (formula_id(value),formulation_options(value))),peer)
                explicit=!ismissing(formula_id(value)) && formula_id(value) ∉ (:default,:none)
                # Explicit branch structure and controls are meaningful even
                # when every candidate shares them or the leaf ID is default.
                (length(last(scope))>1 || !isempty(formulation_options(value)) ||
                    scope in varying || changed || explicit) &&
                    push!(parts,description(scope,value;compact))
            end
        end
        if isempty(parts)
            unavailable=ismissing(sources[index]) ||
                any(entry -> ismissing(formula_id(last(entry))),selections[index])
            push!(parts,unavailable ? description(missing;compact) : "default")
        end
        text=join(parts,"; ")
        isempty(prefix) ? text : isempty(text) ? prefix : prefix*" · "*text
    end
end
