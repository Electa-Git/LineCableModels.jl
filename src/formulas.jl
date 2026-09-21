"""
$(TYPEDSIGNATURES)

Select a registered formula without exposing its owner module or concrete
wrapper type. The receiving formulation determines the formula family from the
keyword slot in which the selection appears.

# Arguments

- `identifier`: Stable formula identifier.
  `:default` routes to the owning family's explicit default implementation for
  the chosen backend. The resulting selection is checked against the problem;
  it is not a fallback after a failed or inapplicable formula.
  Cable-insulation and semicon-admittance `:default` selections route to the
  explicit `:lossless` dielectric relation. Unsupported contexts fail before
  frequency evaluation.

# Keywords

- `order`: Position of an equivalent homogeneous-earth reduction relative to
  material frequency dependence. `:before` applies EquivalentHomogeneous before FrequencyDependent, `:after`
  applies EquivalentHomogeneous after FrequencyDependent, and `:default` selects the receiving formulation's
  default. Non-EquivalentHomogeneous formula slots accept only `:default`.
- `parameters=(;)`: Explicit model parameters accepted by the owning formula.
- `options=(;)`: Formulation-owned physical choices and numerical controls, such
  as Unified's `Γ` \\[1/m\\] and `integration=(method=:quad, options=(;))`.
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
        parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions(), equivalent_earth = nothing)
    options = options isa NamedTuple ? FormulationOptions(options) : options
    order in (:default, :before, :after) || throw(ArgumentError(
        "formula order must be :default, :before, or :after"
    ))
    return FormulaDefinition{identifier, order, typeof(parameters),
        typeof(options), typeof(equivalent_earth)}(
        parameters, options, equivalent_earth)
end

"Return the selected formulation's identifier without evaluating its equation."
formula_id(bound::FormulaMethod) = formula_id(bound.selection)

formula_id(::FormulaDefinition{ID}) where {ID} = ID
formula_id(::Type{<:FormulaDefinition{ID}}) where {ID} = ID

"""Describe an opaque retained identity without claiming or reconstructing an implementation."""
description(::Type{<:FormulaDefinition{ID}}; compact::Bool=false) where {ID} = string(ID)

"""Expose a requested formula identifier and its explicit model and numerical controls."""
function Base.NamedTuple(value::FormulaDefinition{ID,Order}) where {ID,Order}
    return (identifier=ID, order=Order, parameters=value.parameters,
        options=value.options.data, equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end

"""
$(TYPEDSIGNATURES)

Describe a retained selection through its owning type. This does not construct
a formula, evaluate equations, or substitute current defaults.
"""
description(source::Pair{<:Type,<:NamedTuple}; compact::Bool=false) = description(first(source); compact)
formula_id(source::Pair{<:Type,<:NamedTuple}) = formula_id(first(source))

description(::Missing; compact::Bool=false) = "method unavailable"
description(::Nothing; compact::Bool=false) = "none"

"""Describe a scalar, explicitly selected control at its owning formula."""
description(::Type, ::Val{key}, value::Union{Number,Symbol,AbstractString};
    compact::Bool=false) where {key} = string(key,"=",value)
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
identity through `formula_id`, explicit controls in those selection pairs, and
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
    retained(value) = value isa Pair ? last(value) : (;)
    selections=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...; quantity)) for source in sources]
    complete=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...)) for source in sources]
    candidates=findall(!=(:reference),roles)
    # Scientific comparisons use owner-scoped identifiers, never display text.
    scopes=unique([scope for index in candidates for (scope,_) in complete[index]
        if !isempty(last(scope))])
    varying=filter(scopes) do scope
        values=[[(formula_id(value),retained(value)) for (key,value) in complete[index]
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
                controls=retained(value)
                peer=[other for entries in selections for (key,other) in entries if key==scope]
                if !isempty(controls) && any(other -> !isequal(retained(other),controls),peer)
                    push!(parts,description(scope,value;compact,settings=true))
                end
            else
                peer=[other for entries in selections for (key,other) in entries if key==scope]
                changed=any(other -> !isequal((formula_id(other),retained(other)),
                    (formula_id(value),retained(value))),peer)
                standalone=length(sources)==1 && !ismissing(formula_id(value)) &&
                    formula_id(value)!==:none
                # Common scalar choices are omitted independently of their IDs.
                # Branch structure and controls remain visible in comparisons;
                # a standalone description shows its concrete selections.
                (length(last(scope))>1 || !isempty(retained(value)) ||
                    scope in varying || changed || standalone) &&
                    push!(parts,description(scope,value;compact))
            end
        end
        if isempty(parts)
            unavailable=ismissing(sources[index]) ||
                any(entry -> ismissing(formula_id(last(entry))),selections[index])
            push!(parts,unavailable ? description(missing;compact) : description(sources[index];compact))
        end
        text=join(parts,"; ")
        isempty(prefix) ? text : isempty(text) ? prefix : prefix*" · "*text
    end
end
