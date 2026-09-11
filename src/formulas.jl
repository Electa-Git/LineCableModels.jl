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

"""Expose a requested formula identifier and its explicit parameters and overrides."""
function Base.NamedTuple(value::FormulaDefinition{ID,Order}) where {ID,Order}
    return (identifier=ID, order=Order, parameters=value.parameters, hooks=value.hooks,
        options=value.options, equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end

"""
$(TYPEDSIGNATURES)

Describe complete formulation records with stable indices. Common fields are
omitted from curve labels and retained in the supplied records. Every differing
nested field, including air/earth/mixed choices and numerical controls, remains
identifiable. This operation never selects an equation or merges equal curves.
"""
function description(records::AbstractVector; indices=collect(eachindex(records)), prefix="F")
    length(indices) == length(records) || throw(DimensionMismatch("one index is required per formulation"))
    flatten = function (record)
        fields=Pair{String,String}[]
        function visit(value,path)
            if value isa NamedTuple
                for (key,item) in pairs(value)
                    visit(item,isempty(path) ? string(key) : path*"."*string(key))
                end
            elseif value isa AbstractDict
                for key in sort!(collect(keys(value));by=string)
                    visit(value[key],isempty(path) ? string(key) : path*"."*string(key))
                end
            elseif value isa Union{Nothing,Missing,Number,AbstractString,Symbol,Bool}
                push!(fields,path=>sprint(show,value;context=:compact=>true))
            elseif value isa Type
                visit(sprint(show,value;context=(:module=>nothing,:compact=>false)),path)
            elseif value isa Union{Tuple,AbstractArray}
                for (index,item) in enumerate(value)
                    visit(item,path*"[$index]")
                end
            else
                visit(sprint(show,typeof(value);context=(:module=>nothing,:compact=>false)),path*".type")
                for key in fieldnames(typeof(value))
                    visit(getfield(value,key),path*".fields."*string(key))
                end
            end
        end
        if record isa NamedTuple && haskey(record,:requested)
            visible=NamedTuple{Tuple(key for key in keys(record) if key in (:backend,:requested,:options,:execution))}(
                Tuple(value for (key,value) in pairs(record) if key in (:backend,:requested,:options,:execution)))
            visit(visible,"")
        else
            visit(record,"")
        end
        return Dict(fields)
    end
    fields=map(flatten,records)
    paths=sort!(unique([key for record in fields for key in keys(record)]))
    differences=filter(paths) do path
        if length(records)==1
            value=get(first(fields),path,nothing)
            return path == "backend" ||
                ((startswith(path,"requested.") || endswith(path,".identifier")) && value ∉ (":default","nothing"))
        end
        !all(record -> get(record,path,nothing) == get(first(fields),path,nothing),fields)
    end
    names=map(differences) do path
        replace(path,r"^requested\."=>"",r"\.identifier$"=>"",
            "internal_impedance"=>"internal Z","insulation_impedance"=>"insulation Z",
            "insulation_admittance"=>"insulation Y","semicon_admittance"=>"semicon Y",
            "earth_impedance"=>"earth Z","earth_admittance"=>"earth Y",
            "pipe_impedance"=>"pipe Z","earth_properties"=>"soil law",
            "temperature_dependence"=>"temperature law","equivalent_earth"=>"equivalent earth",
            ".options."=>".",".integration.method"=>".integration")
    end
    return map(enumerate(fields)) do (index,record)
        selected=map(zip(differences,names)) do (path,name)
            value=get(record,path,"unspecified")
            startswith(value,":") && (value=chop(value;head=1,tail=0))
            "$name=$value"
        end
        "$prefix$(indices[index])" * (isempty(selected) ? "" : " · " * join(selected,"; "))
    end
end
