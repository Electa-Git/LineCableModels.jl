"""
$(TYPEDSIGNATURES)

Return the stored unit for a physical quantity.

The generic method returns a dimensionless unit. Calculation modules define
methods for the quantities they publish.
"""
default_unit(::QuantityTag{Q}) where {Q} = Units()

"""
$(TYPEDSIGNATURES)

Return the default display unit for a physical quantity.

The generic method returns [`default_unit`](@ref). Calculation modules may
select a different metric prefix for display.
"""
display_unit(q::QuantityTag{Q}) where {Q} = default_unit(q)

function _with_basis(unit::Units, basis::Symbol; length_prefix::Symbol = :base)
    basis in (:per_length, :total) || throw(
        ArgumentError("basis must be :per_length or :total"),
    )
    has_length_denominator = any(item -> item.name == :meter, unit.per)
    retained = Unit[item for item in unit.per if item.name != :meter]
    if basis === :per_length && has_length_denominator
        push!(retained, Unit(name = :meter, prefix = length_prefix))
    end
    return Units(base = copy(unit.base), per = retained)
end

function default_unit(quantity::QuantityTag, basis::Symbol)
    _with_basis(default_unit(quantity), basis; length_prefix = :base)
end
function display_unit(
        quantity::QuantityTag,
        basis::Symbol;
        length_prefix::Symbol = :kilo
)
    _with_basis(display_unit(quantity), basis; length_prefix)
end

function scale_factor(quantity::QuantityTag, basis::Symbol, target::Units)
    return scale_factor(default_unit(quantity, basis), target)
end

"""
$(TYPEDSIGNATURES)

Return the quantity label without units.

The generic method returns the quantity-tag symbol as text.
"""
get_label(::QuantityTag{Q}) where {Q} = String(Q)

"""
$(TYPEDSIGNATURES)

Return the mathematical symbol used to label a quantity.

The generic method returns the quantity-tag symbol as text.
"""
function get_symbol(::QuantityTag{Q}) where {Q}
    Q isa Tuple ? get_symbol(QuantityTag{first(Q)}()) : String(Q)
end

# --------------------------------------------------------------------------
# Unit scaling between arbitrary composite units
# --------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return the metric-prefix scale from `from_unit` to `to_unit`:

```math
v_{\\mathrm{to}} =
v_{\\mathrm{from}} 10^{e_{\\mathrm{from}}-e_{\\mathrm{to}}}.
```

The unit names must already describe compatible physical dimensions. The
method converts prefixes only.
"""
function scale_factor(from_unit::Units, to_unit::Units)
    ex_from = get_exp(from_unit)
    ex_to = get_exp(to_unit)
    return 10.0^(ex_from - ex_to)
end

"""
$(TYPEDSIGNATURES)

Return the metric-prefix scale from the stored unit of `q` to `to_unit`.
"""
function scale_factor(q::QuantityTag, to_unit::Units)
    return scale_factor(default_unit(q), to_unit)
end

# --------------------------------------------------------------------------
# Prefix exponent aggregation for composite units
# --------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return the net base-10 exponent of a composite unit:

```math
e = \\sum_{u \\in \\mathrm{base}} e_u -
    \\sum_{u \\in \\mathrm{per}} e_u.
```
"""
function get_exp(u::Units)::Int
    num = 0
    for b in u.base
        num += _prefix_exp(b.prefix)
    end
    den = 0
    for p in u.per
        den += _prefix_exp(p.prefix)
    end
    return num - den
end

"""
$(TYPEDSIGNATURES)

Return the scale from unprefixed units to `u`.
"""
scale_factor(u::Units) = 10.0^(-get_exp(u))

# Accept `Val{:frequency}` and the other scientific quantity keys.
default_unit(::Val{Q}) where {Q} = default_unit(QuantityTag{Q}())

display_unit(::Val{Q}) where {Q} = display_unit(QuantityTag{Q}())

get_label(::Val{Q}) where {Q} = get_label(QuantityTag{Q}())

# Symbol methods dispatch through Val.
default_unit(q::Symbol) = default_unit(Val(q))

display_unit(q::Symbol) = display_unit(Val(q))

get_label(q::Symbol) = get_label(Val(q))
