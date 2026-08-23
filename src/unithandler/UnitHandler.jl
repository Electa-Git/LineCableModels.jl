module UnitHandler

using Base: @kwdef
import ..LineCableModels: Z, Y, R, X, L, G, B, C

export Unit, Units, units, get_label, get_symbol, get_exp,
       METRIC_PREFIX_EXPONENT, METRIC_PREFIX_SYMBOL, UNIT_SYMBOL,
       QuantityTag, quantity, default_unit, display_unit, scale_factor

# --------------------------------------------------------------------------
# Metric prefixes
# --------------------------------------------------------------------------

const METRIC_PREFIX_EXPONENT = Dict(
    :yocto => -24,
    :zepto => -21,
    :atto => -18,
    :femto => -15,
    :pico => -12,
    :nano => -9,
    :micro => -6,
    :milli => -3,
    :centi => -2,
    :deci => -1,
    :base => 0,
    :deca => 1,
    :hecto => 2,
    :kilo => 3,
    :mega => 6,
    :giga => 9,
    :tera => 12,
    :peta => 15,
    :exa => 18,
    :zetta => 21,
    :yotta => 24
)

const METRIC_PREFIX_SYMBOL = Dict(
    :yocto => "y",
    :zepto => "z",
    :atto => "a",
    :femto => "f",
    :pico => "p",
    :nano => "n",
    :micro => "μ",
    :milli => "m",
    :centi => "c",
    :deci => "d",
    :base => "",
    :deca => "da",
    :hecto => "h",
    :kilo => "k",
    :mega => "M",
    :giga => "G",
    :tera => "T",
    :peta => "P",
    :exa => "E",
    :zetta => "Z",
    :yotta => "Y"
)

const UNIT_SYMBOL = Dict(
    :ohm => "Ω",
    :henry => "H",
    :farad => "F",
    :siemens => "S",
    :meter => "m",
    :hertz => "Hz",
    :degree => "°",
    :dimensionless => ""    # extend this as your sadism requires
)

@inline function _prefix_exp(prefix::Symbol)
    return get(METRIC_PREFIX_EXPONENT, prefix) do
        throw(ArgumentError("unsupported metric prefix :$prefix"))
    end
end

@inline function _prefix_symbol(prefix::Symbol)
    return get(METRIC_PREFIX_SYMBOL, prefix) do
        throw(ArgumentError("unsupported metric prefix :$prefix"))
    end
end

# --------------------------------------------------------------------------
# Core types
# --------------------------------------------------------------------------

"""
A single physical unit with a metric prefix.

- `name`   : symbolic unit name, e.g. :ohm, :henry, :farad, :siemens, :meter, :hertz, :degree, :dimensionless
- `prefix` : one of the metric prefix symbols (:base, :kilo, :milli, ...).
"""
@kwdef struct Unit
    name::Symbol = :dimensionless
    prefix::Symbol = :base
end

"""
A composite unit, e.g., "Ω/km".

- `base` : numerator units.
- `per`  : denominator units.
"""
@kwdef struct Units
    base::Vector{Unit} = [Unit()]  # dimensionless by default
    per::Vector{Unit} = Unit[]
end

# --------------------------------------------------------------------------
# Convenience constructor
# --------------------------------------------------------------------------

"""
Convenience constructor for `Units`.

Examples:

    u1 = units(:base, :ohm)                          # Ω
    u2 = units(:base, :ohm; per = (:kilo, :meter))   # Ω/km
"""
function units(
        prefix::Symbol,
        name::Symbol;
        per::Union{Nothing, Tuple{Symbol, Symbol}} = nothing
)
    b = Unit(name = name, prefix = prefix)
    if per === nothing
        return Units(base = [b], per = Unit[])
    else
        pfx2, name2 = per
        d = Unit(name = name2, prefix = pfx2)
        return Units(base = [b], per = [d])
    end
end

# --------------------------------------------------------------------------
# Formatting
# --------------------------------------------------------------------------

"""
Render a `Unit` as a readable string, e.g. "kΩ", "mH", "μF", "Hz".
"""
function get_label(u::Unit)
    prefix_str = _prefix_symbol(u.prefix)
    unit_str = get(UNIT_SYMBOL, u.name, String(u.name))
    return string(prefix_str, unit_str)
end

"""
Render composite `Units` as a readable string:

- base units concatenated with "."
- per units concatenated with "."
- if more than one per, denominator is shown as "(...)"
- :dimensionless units disappear; all-dimensionless → ""
"""
function get_label(u::Units)
    base_units = [x for x in u.base if x.name != :dimensionless]
    per_units = [x for x in u.per if x.name != :dimensionless]

    base_strs = [get_label(b) for b in base_units if get_label(b) != ""]
    per_strs = [get_label(p) for p in per_units if get_label(p) != ""]

    if isempty(base_strs) && isempty(per_strs)
        return ""
    end

    base_str = isempty(base_strs) ? "" : join(base_strs, ".")

    if isempty(per_strs)
        return base_str
    else
        if isempty(base_str)
            base_str = "1"
        end
        if length(per_strs) == 1
            return string(base_str, "/", per_strs[1])
        else
            return string(base_str, "/(", join(per_strs, "."), ")")
        end
    end
end

# --------------------------------------------------------------------------
# Quantity semantics
# --------------------------------------------------------------------------

"""
QuantityTag{Q} is a semantic tag for a physical quantity.

Examples of Q (to be extended as needed):

- :frequency, :series_impedance, :resistance, :inductance
- :shunt_admittance, :conductance, :capacitance
- :angle, :index, :distance, :magnitude...

These tags are used across plotting, exporting, reporting, etc.
"""
struct QuantityTag{Q} end

QuantityTag(::Val{Q}) where {Q} = QuantityTag{Q}()
QuantityTag(::Type{QuantityTag{Q}}) where {Q} = QuantityTag{Q}()
_quantity_name(::QuantityTag{Q}) where {Q} = Q

"""
    quantity(accessor[, selector])

Return the physical quantity represented by a result accessor. `selector` may
be `:re`, `:im`, `:abs`, or `:angle` for the complex `Z` and `Y` accessors.

# Arguments

- `accessor`: One of `Z`, `Y`, `R`, `X`, `L`, `G`, `B`, or `C`.
- `selector`: Optional Cartesian or polar component of `Z` or `Y`.

# Returns

- A `QuantityTag` used to resolve primary units, display units, labels, and
  scaling.

# Errors

Throws `ArgumentError` for an unsupported accessor or selector.

# Examples

```julia
LineCableModels.UnitHandler.quantity(Z)       # series impedance
LineCableModels.UnitHandler.quantity(Z, :re)  # series resistance
LineCableModels.UnitHandler.quantity(Y, :abs) # shunt-admittance magnitude
```
"""
quantity(::typeof(Z)) = QuantityTag{:series_impedance}()
quantity(::typeof(Y)) = QuantityTag{:shunt_admittance}()
quantity(::typeof(R)) = QuantityTag{:resistance}()
quantity(::typeof(X)) = QuantityTag{:reactance}()
quantity(::typeof(L)) = QuantityTag{:inductance}()
quantity(::typeof(G)) = QuantityTag{:conductance}()
quantity(::typeof(B)) = QuantityTag{:susceptance}()
quantity(::typeof(C)) = QuantityTag{:capacitance}()

quantity(::typeof(Z), ::Val{:re}) = quantity(R)
quantity(::typeof(Z), ::Val{:im}) = quantity(X)
quantity(::typeof(Z), ::Val{:abs}) = QuantityTag{(:series_impedance, :abs)}()
quantity(::typeof(Z), ::Val{:angle}) = QuantityTag{(:series_impedance, :angle)}()
quantity(::typeof(Y), ::Val{:re}) = quantity(G)
quantity(::typeof(Y), ::Val{:im}) = quantity(B)
quantity(::typeof(Y), ::Val{:abs}) = QuantityTag{(:shunt_admittance, :abs)}()
quantity(::typeof(Y), ::Val{:angle}) = QuantityTag{(:shunt_admittance, :angle)}()

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

const _SCIENTIFIC_QUANTITY = Dict(
    :frequency => (QuantityTag{:frequency}(), :hertz, :base),
    :series_impedance => (quantity(Z), :ohm, :base),
    :shunt_admittance => (quantity(Y), :siemens, :base),
    :resistance => (quantity(R), :ohm, :base),
    :reactance => (quantity(X), :ohm, :base),
    :inductance => (quantity(L), :henry, :milli),
    :conductance => (quantity(G), :siemens, :base),
    :susceptance => (quantity(B), :siemens, :base),
    :capacitance => (quantity(C), :farad, :micro),
    :angle => (QuantityTag{:angle}(), :degree, :base),
    :series_impedance_absolute_error =>
        (QuantityTag{:series_impedance_absolute_error}(), :ohm, :base),
    :series_impedance_relative_error =>
        (QuantityTag{:series_impedance_relative_error}(), :dimensionless, :base),
    :shunt_admittance_absolute_error =>
        (QuantityTag{:shunt_admittance_absolute_error}(), :siemens, :base),
    :shunt_admittance_relative_error =>
        (QuantityTag{:shunt_admittance_relative_error}(), :dimensionless, :base)
)

"""
    line_component_quantity(key)

Return the quantity tag, physical unit name, and default display prefix for a
scientific observable key.
"""
function line_component_quantity(key::Symbol)
    tag, unit_name, prefix = get(_SCIENTIFIC_QUANTITY, key) do
        throw(ArgumentError("unsupported scientific quantity :$key"))
    end
    return (;
        semantic = key,
        tag,
        unit_name,
        prefix
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
    line_component_unit(component, basis; length_unit=:kilo, quantity_units=nothing)

Resolve the display unit and numeric scale factor for a line-parameter
component. `length_unit` is the metric prefix of the denominator for
per-length quantities. `quantity_units` may override numerator prefixes with a
single symbol or with component/semantic keys in a named tuple or dictionary.
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

"""
Native (storage/computational) unit for a quantity.

Default: dimensionless, must be overridden in domain code where appropriate.
"""
default_unit(::QuantityTag{Q}) where {Q} = Units()

"""
Default display unit for a quantity.

By default, equal to `default_unit(q)`, but you override this when
you want prettier axes/reports.

Example override (outside this module):

    default_unit(::QuantityTag{:Z}) = units(:base, :ohm; per = (:base, :meter))   # Ω/m
    display_unit(::QuantityTag{:Z}) =
        units(:base, :ohm; per = (:kilo, :meter))                                # Ω/km
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
Human-readable label for the quantity, without units.

Default: `String(Q)`, i.e., :freq → "freq". You override this with
more civilized labels.
"""
get_label(::QuantityTag{Q}) where {Q} = String(Q)

"""
SI symbol label for the quantity.

Default: fall back to the quantity tag name `Q` as a string.
There is no way to guess, so you must override it with meaningful labels.
"""
function get_symbol(::QuantityTag{Q}) where {Q}
    Q isa Tuple ? get_symbol(QuantityTag{first(Q)}()) : String(Q)
end

# --------------------------------------------------------------------------
# Unit scaling between arbitrary composite units
# --------------------------------------------------------------------------

"""
Numeric scale factor to convert values from `from_unit` to `to_unit`.

If `v_raw` is expressed in `from_unit`, then:

    v_to = v_raw * scale_factor(from_unit, to_unit)

because 1[to] = 10^(exp_to - exp_from) [from] ⇒ v_to = v_from * 10^(exp_from - exp_to).
"""
function scale_factor(from_unit::Units, to_unit::Units)
    ex_from = get_exp(from_unit)
    ex_to = get_exp(to_unit)
    return 10.0^(ex_from - ex_to)
end

"""
Scale factor to convert from the native unit of quantity `q` to a
target display unit `to_unit`:

    v_display = v_native * scale_factor(q, to_unit)
"""
function scale_factor(q::QuantityTag, to_unit::Units)
    return scale_factor(default_unit(q), to_unit)
end

# --------------------------------------------------------------------------
# Prefix exponent aggregation for composite units
# --------------------------------------------------------------------------

"""
Return the net base-10 exponent from the prefixes in a composite unit:

    exp = (sum prefix exponents over `base`) - (sum over `per`)
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
Scale factor to convert a raw value expressed in base-prefixed units
(e.g. Ω/m, S/m, Hz) to the given `Units`.

i.e. q_disp = q_raw * scale_factor(u)
"""
scale_factor(u::Units) = 10.0^(-get_exp(u))

# Convenience: call with `Val{:frequency}` and other scientific keys.
default_unit(::Val{Q}) where {Q} = default_unit(QuantityTag{Q}())

display_unit(::Val{Q}) where {Q} = display_unit(QuantityTag{Q}())

get_label(::Val{Q}) where {Q} = get_label(QuantityTag{Q}())

# Convenience: call with a Symbol. This just wraps to Val.
default_unit(q::Symbol) = default_unit(Val(q))

display_unit(q::Symbol) = display_unit(Val(q))

get_label(q::Symbol) = get_label(Val(q))

# --------------------------------------------------------------------------
# Fundamental quantities
# --------------------------------------------------------------------------

# Frequency
default_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
display_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
get_label(::QuantityTag{:frequency}) = "Frequency"
get_symbol(::QuantityTag{:frequency}) = "f"

# Series resistance
default_unit(::QuantityTag{:resistance}) = units(:base, :ohm; per = (:base, :meter))    # Ω/m
display_unit(::QuantityTag{:resistance}) = units(:base, :ohm; per = (:kilo, :meter))    # Ω/km
get_label(::QuantityTag{:resistance}) = "Series resistance"
get_symbol(::QuantityTag{:resistance}) = "R"

# Series inductance
default_unit(::QuantityTag{:inductance}) = units(:base, :henry; per = (:base, :meter))
display_unit(::QuantityTag{:inductance}) = units(:milli, :henry; per = (:kilo, :meter))
get_label(::QuantityTag{:inductance}) = "Series inductance"
get_symbol(::QuantityTag{:inductance}) = "L"

# Shunt capacitance
default_unit(::QuantityTag{:capacitance}) = units(:base, :farad; per = (:base, :meter))
display_unit(::QuantityTag{:capacitance}) = units(:micro, :farad; per = (:kilo, :meter))
get_label(::QuantityTag{:capacitance}) = "Shunt capacitance"
get_symbol(::QuantityTag{:capacitance}) = "C"

# Shunt conductance
default_unit(::QuantityTag{:conductance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:conductance}) = units(:base, :siemens; per = (:kilo, :meter))
get_label(::QuantityTag{:conductance}) = "Shunt conductance"
get_symbol(::QuantityTag{:conductance}) = "G"

# Series impedance
default_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:kilo, :meter))
get_label(::QuantityTag{:series_impedance}) = "Series impedance"
get_symbol(::QuantityTag{:series_impedance}) = "Z"

# Shunt admittance
function default_unit(::QuantityTag{:shunt_admittance})
    units(:base, :siemens; per = (:base, :meter))
end
function display_unit(::QuantityTag{:shunt_admittance})
    units(:base, :siemens; per = (:kilo, :meter))
end
get_label(::QuantityTag{:shunt_admittance}) = "Shunt admittance"
get_symbol(::QuantityTag{:shunt_admittance}) = "Y"

# Inductive reactance
default_unit(::QuantityTag{:reactance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:reactance}) = units(:base, :ohm; per = (:kilo, :meter))
get_label(::QuantityTag{:reactance}) = "Inductive reactance"
get_symbol(::QuantityTag{:reactance}) = "X"

# Capacitive susceptance
default_unit(::QuantityTag{:susceptance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:susceptance}) = units(:base, :siemens; per = (:kilo, :meter))
get_label(::QuantityTag{:susceptance}) = "Capacitive susceptance"
get_symbol(::QuantityTag{:susceptance}) = "B"

# Angle
default_unit(::QuantityTag{:angle}) = units(:base, :degree)
display_unit(::QuantityTag{:angle}) = units(:base, :degree)
get_label(::QuantityTag{:angle}) = "Angle"
get_symbol(::QuantityTag{:angle}) = "∠"

# magnitude uses same unit as base quantity
function default_unit(::QuantityTag{(:series_impedance, :re)})
    default_unit(QuantityTag{:resistance}())
end
function default_unit(::QuantityTag{(:series_impedance, :im)})
    default_unit(QuantityTag{:reactance}())
end
function default_unit(::QuantityTag{(:series_impedance, :abs)})
    default_unit(QuantityTag{:series_impedance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :re)})
    default_unit(QuantityTag{:conductance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :im)})
    default_unit(QuantityTag{:susceptance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :abs)})
    default_unit(QuantityTag{:shunt_admittance}())
end

# angle unit
function default_unit(::QuantityTag{(:series_impedance, :angle)})
    default_unit(QuantityTag{:angle}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :angle)})
    default_unit(QuantityTag{:angle}())
end

# labels
get_label(::QuantityTag{(:series_impedance, :abs)}) = "Series impedance magnitude"
get_label(::QuantityTag{(:series_impedance, :angle)}) = "Series impedance angle"
get_label(::QuantityTag{(:shunt_admittance, :abs)}) = "Shunt admittance magnitude"
get_label(::QuantityTag{(:shunt_admittance, :angle)}) = "Shunt admittance angle"

for (key, source, label, symbol) in (
    (:series_impedance_absolute_error, :series_impedance,
    "Series impedance absolute error", "ΔZ"),
    (:shunt_admittance_absolute_error, :shunt_admittance,
    "Shunt admittance absolute error", "ΔY")
)
    @eval begin
        default_unit(::QuantityTag{$(QuoteNode(key))}) = default_unit(QuantityTag{$(QuoteNode(source))}())
        display_unit(::QuantityTag{$(QuoteNode(key))}) = display_unit(QuantityTag{$(QuoteNode(source))}())
        get_label(::QuantityTag{$(QuoteNode(key))}) = $label
        get_symbol(::QuantityTag{$(QuoteNode(key))}) = $symbol
    end
end

for (key, label, symbol) in (
    (:series_impedance_relative_error, "Series impedance relative error", "εZ"),
    (:shunt_admittance_relative_error, "Shunt admittance relative error", "εY")
)
    @eval begin
        default_unit(::QuantityTag{$(QuoteNode(key))}) = units(:base, :dimensionless)
        display_unit(::QuantityTag{$(QuoteNode(key))}) = units(:base, :dimensionless)
        get_label(::QuantityTag{$(QuoteNode(key))}) = $label
        get_symbol(::QuantityTag{$(QuoteNode(key))}) = $symbol
    end
end

end # module UnitHandler
