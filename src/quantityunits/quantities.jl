# --------------------------------------------------------------------------
# Quantity semantics
# --------------------------------------------------------------------------

"""
$(TYPEDEF)

Identify a physical quantity without storing a value.

The type parameter `Q` is a symbol such as `:frequency`, `:resistance`,
`:inductance`, or `:capacitance`. Calculation modules define the mapping from
their accessors to these tags.
"""
struct QuantityTag{Q} end

QuantityTag(::Val{Q}) where {Q} = QuantityTag{Q}()
QuantityTag(::Type{QuantityTag{Q}}) where {Q} = QuantityTag{Q}()
_quantity_name(::QuantityTag{Q}) where {Q} = Q

"""
$(TYPEDSIGNATURES)

Return the physical quantity selected by an owner-defined accessor and
component selector.

# Arguments

- `accessor`: Calculation accessor with a `quantity` method.
- `selector`: One of `:re`, `:im`, `:abs`, or `:angle`.

# Returns

- The corresponding [`QuantityTag`](@ref).

# Errors

- Throws `ArgumentError` when `selector` is unsupported or the accessor does
  not define that selection.
"""

function quantity(accessor, selector::Symbol)
    selector in (:re, :im, :abs, :angle) || throw(
        ArgumentError("selector must be :re, :im, :abs, or :angle"),
    )
    applicable(quantity, accessor, Val(selector)) || throw(
        ArgumentError("selector :$selector is not defined for accessor $accessor"),
    )
    return quantity(accessor, Val(selector))
end

function quantity(accessor)
    throw(ArgumentError("unsupported physical-quantity accessor $accessor"))
end

quantity(::Val{:series_impedance}) = QuantityTag{:series_impedance}()
quantity(::Val{:resistance}) = QuantityTag{:resistance}()
quantity(::Val{:reactance}) = QuantityTag{:reactance}()
quantity(::Val{:inductance}) = QuantityTag{:inductance}()
quantity(::Val{:shunt_admittance}) = QuantityTag{:shunt_admittance}()
quantity(::Val{:conductance}) = QuantityTag{:conductance}()
quantity(::Val{:susceptance}) = QuantityTag{:susceptance}()
quantity(::Val{:capacitance}) = QuantityTag{:capacitance}()

quantity(::Val{:R}) = quantity(Val(:resistance))
quantity(::Val{:X}) = quantity(Val(:reactance))
quantity(::Val{:L}) = quantity(Val(:inductance))
quantity(::Val{:G}) = quantity(Val(:conductance))
quantity(::Val{:B}) = quantity(Val(:susceptance))
quantity(::Val{:C}) = quantity(Val(:capacitance))
quantity(::Val{:Z_re}) = quantity(Val(:resistance))
quantity(::Val{:Z_im}) = quantity(Val(:reactance))
quantity(::Val{:Z_abs}) = quantity(Val(:series_impedance))
quantity(::Val{:Y_re}) = quantity(Val(:conductance))
quantity(::Val{:Y_im}) = quantity(Val(:susceptance))
quantity(::Val{:Y_abs}) = quantity(Val(:shunt_admittance))

quantity(::Val{:frequency}) = QuantityTag{:frequency}()
quantity(::Val{:angle}) = QuantityTag{:angle}()
function quantity(::Val{:series_impedance_absolute_error})
    QuantityTag{:series_impedance_absolute_error}()
end
function quantity(::Val{:series_impedance_relative_error})
    QuantityTag{:series_impedance_relative_error}()
end
function quantity(::Val{:shunt_admittance_absolute_error})
    QuantityTag{:shunt_admittance_absolute_error}()
end
function quantity(::Val{:shunt_admittance_relative_error})
    QuantityTag{:shunt_admittance_relative_error}()
end

quantity(::Val{:Z_angle}) = quantity(Val(:angle))
quantity(::Val{:Y_angle}) = quantity(Val(:angle))

function quantity(key::Symbol)
    applicable(quantity, Val(key)) || throw(
        ArgumentError("unsupported scientific quantity :$key"),
    )
    return quantity(Val(key))
end

"""
$(TYPEDSIGNATURES)

Return the quantity tag, numerator unit name, and default display prefix for a
line-parameter component key.

# Arguments

- `key`: Published line-parameter component, such as `:R` or `:C`.

# Returns

- A named tuple containing `semantic`, `tag`, `unit_name`, and `prefix`.

# Errors

- Throws `ArgumentError` when the key is unknown or its stored/display unit
  has more than one numerator unit.
"""
function line_component_quantity(key::Symbol)
    tag = quantity(key)
    native = default_unit(tag)
    displayed = display_unit(tag)
    length(native.base) == 1 || throw(ArgumentError(
        "scientific quantity :$key must have one numerator unit",
    ))
    length(displayed.base) == 1 || throw(ArgumentError(
        "scientific quantity :$key must have one display numerator unit",
    ))
    return (;
        semantic = _quantity_name(tag),
        tag,
        unit_name = only(native.base).name,
        prefix = only(displayed.base).prefix
    )
end

function _line_quantity_prefix(quantity_units, component, semantic, fallback)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    if quantity_units isa NamedTuple || quantity_units isa AbstractDict
        haskey(quantity_units, component) && return quantity_units[component]
        semantic isa Symbol && haskey(quantity_units, semantic) &&
            return quantity_units[semantic]
        return fallback
    end
    throw(ArgumentError("quantity_units must be a Symbol, NamedTuple, dictionary, or nothing"))
end

"""
$(TYPEDSIGNATURES)

Resolve the display unit and numeric scale factor for one line-parameter
component.

# Arguments

- `key`: Published line-parameter component.
- `basis`: `:per_length` or `:total` storage basis.

# Keywords

- `length_unit`: Metric prefix for the denominator of per-length quantities.
  Default: `:kilo`.
- `quantity_units`: Numerator-prefix override as one symbol or a component
  mapping. Default: `nothing`.

# Returns

- A named tuple containing the quantity tag, display units, and multiplicative
  scale.
"""
function line_component_unit(
        key::Symbol,
        basis::Symbol;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    quantity = line_component_quantity(key)
    prefix = _line_quantity_prefix(
        quantity_units,
        key,
        quantity.semantic,
        quantity.prefix
    )
    prefix isa Symbol || throw(
        ArgumentError("quantity unit prefixes must be symbols"),
    )
    target = if quantity.unit_name in (:degree, :hertz, :dimensionless)
        units(prefix, quantity.unit_name)
    elseif basis === :per_length
        units(prefix, quantity.unit_name; per = (length_unit, :meter))
    else
        units(prefix, quantity.unit_name)
    end
    return (;
        quantity = quantity.tag,
        units = target,
        scale = scale_factor(quantity.tag, basis, target)
    )
end
