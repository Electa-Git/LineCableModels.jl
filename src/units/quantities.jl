"""
$(TYPEDEF)

Identify one scientific quantity without storing its value.
"""
struct QuantityTag{Q} end

QuantityTag(::Val{Q}) where {Q} = QuantityTag{Q}()
QuantityTag(::Type{QuantityTag{Q}}) where {Q} = QuantityTag{Q}()
_quantity_name(::QuantityTag{Q}) where {Q} = Q

"Return the physical quantity selected by a scientific function object."
function quantity end

function quantity(selector, transform::Symbol)
    normalized = transform === :abs ? abs : transform === :angle ? angle : nothing
    normalized === nothing && throw(ArgumentError("transform must be :abs or :angle"))
    return quantity(selector, normalized)
end

# Temporary translation for presentation boundaries migrated in later commits.
quantity(::Val{:frequency}) = QuantityTag{:frequency}()
quantity(::Val{:series_impedance}) = QuantityTag{:series_impedance}()
quantity(::Val{:resistance}) = QuantityTag{:series_resistance}()
quantity(::Val{:reactance}) = QuantityTag{:series_reactance}()
quantity(::Val{:inductance}) = QuantityTag{:series_inductance}()
quantity(::Val{:shunt_admittance}) = QuantityTag{:shunt_admittance}()
quantity(::Val{:conductance}) = QuantityTag{:shunt_conductance}()
quantity(::Val{:susceptance}) = QuantityTag{:shunt_susceptance}()
quantity(::Val{:capacitance}) = QuantityTag{:shunt_capacitance}()
quantity(::Val{:angle}) = QuantityTag{:phase_angle}()
quantity(::Val{:R}) = quantity(Val(:resistance))
quantity(::Val{:X}) = quantity(Val(:reactance))
quantity(::Val{:L}) = quantity(Val(:inductance))
quantity(::Val{:G}) = quantity(Val(:conductance))
quantity(::Val{:B}) = quantity(Val(:susceptance))
quantity(::Val{:C}) = quantity(Val(:capacitance))
quantity(::Val{:Z_re}) = quantity(Val(:resistance))
quantity(::Val{:Z_im}) = quantity(Val(:reactance))
quantity(::Val{:Z_abs}) = QuantityTag{(:series_impedance, :magnitude)}()
quantity(::Val{:Z_angle}) = QuantityTag{(:series_impedance, :phase_angle)}()
quantity(::Val{:Y_re}) = quantity(Val(:conductance))
quantity(::Val{:Y_im}) = quantity(Val(:susceptance))
quantity(::Val{:Y_abs}) = QuantityTag{(:shunt_admittance, :magnitude)}()
quantity(::Val{:Y_angle}) = QuantityTag{(:shunt_admittance, :phase_angle)}()
quantity(::Val{:series_impedance_absolute_error}) =
    QuantityTag{:series_impedance_absolute_error}()
quantity(::Val{:series_impedance_relative_error}) =
    QuantityTag{:series_impedance_relative_error}()
quantity(::Val{:shunt_admittance_absolute_error}) =
    QuantityTag{:shunt_admittance_absolute_error}()
quantity(::Val{:shunt_admittance_relative_error}) =
    QuantityTag{:shunt_admittance_relative_error}()

function quantity(key::Symbol)
    applicable(quantity, Val(key)) || throw(
        ArgumentError("unsupported scientific quantity :$key"),
    )
    return quantity(Val(key))
end

function line_component_quantity(key::Symbol)
    tag = quantity(key)
    native = native_unit(tag)
    displayed = display_unit(tag)
    length(native.numerator) == 1 || throw(
        ArgumentError("scientific quantity :$key must have one numerator unit"),
    )
    length(displayed.numerator) == 1 || throw(
        ArgumentError("scientific quantity :$key must have one display numerator unit"),
    )
    return (;
        semantic = _quantity_name(tag),
        tag,
        unit_name = only(native.numerator).name,
        prefix = only(displayed.numerator).prefix
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

function line_component_unit(
        key::Symbol,
        basis::Symbol;
        length_unit::Symbol = :kilo,
        quantity_units = nothing
)
    specification = line_component_quantity(key)
    prefix = _line_quantity_prefix(
        quantity_units,
        key,
        specification.semantic,
        specification.prefix
    )
    prefix isa Symbol || throw(ArgumentError("quantity unit prefixes must be symbols"))
    target = if specification.unit_name in (:degree, :radian, :hertz, :dimensionless)
        units(prefix, specification.unit_name)
    elseif basis === :per_length
        units(prefix, specification.unit_name; per = (length_unit, :meter))
    elseif basis === :total
        units(prefix, specification.unit_name)
    else
        throw(ArgumentError("basis must be :per_length or :total"))
    end
    return (;
        quantity = specification.tag,
        units = target,
        scale = scale_factor(specification.tag, basis, target)
    )
end
