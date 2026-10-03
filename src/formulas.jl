"""
$(TYPEDSIGNATURES)

Select a registered formula. The receiving formulation determines its family from the
keyword slot in which the selection appears.

# Arguments

- `identifier`: stable formula identifier.
  `:default` routes to the defining family's explicit default implementation for
  the chosen backend. The resulting selection is checked against the problem
  before computation.
  Cable-insulation and semicon-admittance `:default` selections route to the
  explicit `:lossless` dielectric relation. Unsupported contexts fail before
  frequency evaluation.

# Keywords

- `order`: position of an equivalent homogeneous-earth reduction relative to
  material frequency dependence. `:before` applies EquivalentHomogeneous before FrequencyDependent, `:after`
  applies EquivalentHomogeneous after FrequencyDependent, and `:default` selects the receiving formulation's
  default. Non-EquivalentHomogeneous formula slots accept only `:default`.
- `parameters=(;)`: explicit model parameters accepted by the defining formula.
- `options=(;)`: formulation-owned physical choices and numerical controls, such
  as Unified's `Γ` \\[1/m\\] and `integration=(method=:quad, options=(;))`.
- `equivalent_earth=nothing`: explicit reduction for a compatible external formula.

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

"Return the selected formulation's identifier."
formula_id(bound::FormulaMethod) = formula_id(bound.selection)

formula_id(::FormulaDefinition{ID}) where {ID} = ID
formula_id(::Type{<:FormulaDefinition{ID}}) where {ID} = ID

"""Describe a retained formula identifier."""
description(::Type{<:FormulaDefinition{ID}}; compact::Bool=false) where {ID} = string(ID)

"""Expose a requested formula identifier and its explicit model and numerical controls."""
function Base.NamedTuple(value::FormulaDefinition{ID,Order}) where {ID,Order}
    return (identifier=ID, order=Order, parameters=value.parameters,
        options=value.options.data, equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end

"""
$(TYPEDSIGNATURES)

Describe a retained selection through its defining type.
"""
description(source::Pair{<:Type,<:NamedTuple}; compact::Bool=false) = description(first(source); compact)
formula_id(source::Pair{<:Type,<:NamedTuple}) = formula_id(first(source))

description(::Missing; compact::Bool=false) = "method unavailable"
description(::Nothing; compact::Bool=false) = "none"

"""Describe a scalar, explicitly selected control at its defining formula."""
description(::Type, ::Val{key}, value::Union{Number,Symbol,AbstractString};
    compact::Bool=false) where {key} = string(key,"=",value)
formula_id(::Missing) = missing
formula_id(::Nothing) = :none

"""
$(TYPEDSIGNATURES)

Label ordered formulation selections using their owner-supplied descriptions.

# Arguments

- `sources`: native formulations or owner-bound retained declarations.

# Keywords

- `roles`: one `:reference`, `:result`, or `:none` role per source.
- `quantity`: physical quantity selected through the observation grammar, or
  `nothing` for all equation choices.
- `compact=true`: use the short owned names for both root and child selections.

# Returns

- One text label per source, in input order, without result numbering.
  Consumers select unique quantity-relevant formulations before presentation.

# Notes

Owners expose native children through `pairs(source; quantity)` and retained
children through `pairs(owner, record; quantity)`, scoped
identity through `formula_id`, explicit controls in those selection pairs, and
text through `description`. This formatter only decides common-field omission
and reference prefixes. Unsupported owner methods are not caught.
"""
function description(sources::AbstractVector;
        roles=fill(:result,length(sources)),
        quantity=nothing, compact::Bool=true)
    length(roles)==length(sources) ||
        throw(DimensionMismatch("one role is required per formulation"))
    all(in((:reference,:result,:none)),roles) ||
        throw(ArgumentError("description roles must be reference, result, or none"))
    isempty(sources) && return String[]
    retained(value) = value isa Pair ? last(value) : (;)
    selections=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...; quantity)) for source in sources]
    complete=[ismissing(source) ? Pair[] : collect(pairs((source isa Pair ? Tuple(source) : (source,))...)) for source in sources]
    formulation_indices=findall(!=(:reference),roles)
    # Scientific comparisons use owner-scoped identifiers, never display text.
    scopes=unique([scope for index in formulation_indices for (scope,_) in complete[index]
        if !isempty(last(scope))])
    varying=filter(scopes) do scope
        values=[[(formula_id(value),retained(value)) for (key,value) in complete[index]
            if key==scope] for index in formulation_indices]
        !all(value -> isequal(value,first(values)),values)
    end
    identities=[[scope for (scope,_) in entries if isempty(last(scope))] for entries in complete]
    formulation_owners=unique(first(ids) for ids in identities[formulation_indices] if !isempty(ids))
    inner_owners=unique(last(ids) for ids in identities if length(ids)>1)
    return map(eachindex(sources)) do index
        prefix=roles[index]===:reference ? "Reference" : ""
        parts=String[]
        for (scope,value) in selections[index]
            if isempty(last(scope))
                position=findfirst(==(scope),identities[index])
                show_identity=position==1 ?
                    (roles[index]!==:result || length(identities[index])>1 || length(formulation_owners)>1) :
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
                # Branch structure and controls remain visible in comparisons.
                # A standalone description shows its concrete selections.
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
