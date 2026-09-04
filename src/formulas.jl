"Return a short scientific description of a registered formulation."
function description end

"""
$(TYPEDEF)

Store one declarative formula selection until its owning formulation resolves
the identifier and overrides into a concrete formula type.

`FormulaSpec` is produced by [`formula`](@ref). It does not participate in a
numerical loop.

$(TYPEDFIELDS)
"""
struct FormulaSpec{ID, Order, O <: NamedTuple}
    "Formula-specific route or assumption overrides."
    overrides::O
end

"""
$(TYPEDEF)

Bind one formula identity and optional semantic selectors to a domain method.

Calling the bound method inserts `Val(ID)` before the stored selectors and
runtime arguments. Formula catalogues use this invariant to retain owner-local
dispatch while carrying the selected literature identity as concrete type
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

# Keywords

- `order`: Position of an equivalent homogeneous-earth reduction relative to
  material frequency dependence. `:before` applies EHEM before FD, `:after`
  applies EHEM after FD, and `:default` selects the receiving formulation's
  default. Non-EHEM formula slots accept only `:default`.
- `kwargs`: Formula-specific route or assumption overrides.

# Returns

- A concrete declarative selection resolved before computation.

# Examples

```julia
earth = formula(:Papadopoulos2010)
soil = formula(:CIGRE2019; epsilon_infinity=10.0)
equivalent = formula(:Xue2021; order=:before)
```
"""
function formula(identifier::Symbol; order::Symbol = :default, kwargs...)
    order in (:default, :before, :after) || throw(ArgumentError(
        "formula order must be :default, :before, or :after"
    ))
    overrides = (; kwargs...)
    return FormulaSpec{identifier, order, typeof(overrides)}(overrides)
end

"Return the stable literature identifier of a formula value."
function formula_id end

formula_id(::FormulaSpec{ID}) where {ID} = ID
