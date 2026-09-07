"Return the series-impedance values of a line-parameter result."
function Z end

"Return the shunt-admittance values of a line-parameter result."
function Y end

function X end
function G end
function B end
function series_impedance end
function shunt_admittance end
function reactance end
function conductance end
function susceptance end
function frequencies end
function nconductors end
function nfrequencies end
"Return absolute numerical errors from an owned comparison result."
function absolute_error end
"Return reference-normalised numerical errors from an owned comparison result."
function relative_error end

"Abstract tag for the physical domain represented by line-parameter matrices."
abstract type LineParamsDomain end

"""
$(TYPEDEF)

Store one earth-return matrix interaction and its resolved physical layers.

$(TYPEDFIELDS)
"""
struct EarthPair{T <: Real}
    "Destination row."
    row::Int
    "Destination column."
    column::Int
    "Conductor heights relative to the air-earth interface \\[m\\]."
    heights::Tuple{T, T}
    "Horizontal separation, or cable outer radius for a self interaction \\[m\\]."
    separation::T
    "Physical layer indices of the source and target conductors."
    layers::Tuple{Int, Int}
end

"""
$(TYPEDSIGNATURES)

Check the resolved geometry of an earth-return interaction.

# Arguments

- `pair`: Matrix indices, signed heights, horizontal separation (self radius),
  and physical layer indices. Lengths are in \\[m\\]; layer 1 is air.

# Returns

- The same `pair`, without altering its geometry.

# Errors

- Throws `ArgumentError` for inconsistent indices or layer assignments.
- Throws `DomainError` for nonfinite lengths, nonpositive self radius,
  coincident distinct conductors, or conductors on the air-earth interface.
"""
function validate(pair::EarthPair)
    pair.row > 0 && pair.column > 0 || throw(ArgumentError(
        "earth-pair matrix indices must be positive"))
    all(>(0), pair.layers) || throw(ArgumentError(
        "earth-pair physical layer indices must be positive; layer 1 is air"))
    all(isfinite, pair.heights) && isfinite(pair.separation) || throw(DomainError(
        (pair.heights, pair.separation), "earth-pair lengths must be finite"))
    pair.separation >= zero(pair.separation) || throw(DomainError(
        pair.separation, "earth-pair horizontal separation must be nonnegative"))
    for index in eachindex(pair.heights)
        height = pair.heights[index]
        iszero(height) && throw(DomainError(
            height, "a conductor on the air-earth interface has no physical layer"))
        (height > zero(height)) == (pair.layers[index] == 1) || throw(ArgumentError(
            "earth-pair conductor $index height and physical layer disagree; layer 1 is air"))
    end
    if pair.row == pair.column
        pair.heights[1] == pair.heights[2] && pair.layers[1] == pair.layers[2] ||
            throw(ArgumentError("an earth self interaction must repeat the same conductor"))
        pair.separation > zero(pair.separation) || throw(DomainError(
            pair.separation, "an earth self interaction requires a positive cable outer radius"))
    else
        iszero(pair.separation) && pair.heights[1] == pair.heights[2] &&
            throw(DomainError((pair.heights, pair.separation),
                "distinct earth-return conductors cannot have coincident centres"))
    end
    return pair
end

"""
$(TYPEDSIGNATURES)

Check an earth pair before invoking a user-supplied route. The generic method
checks geometry only: a custom callable owns its additional restrictions and
may specialize this method without executing its numerical calculation.
Registered native routes provide their own methods beside their equations.

# Arguments

- `pair`: Resolved earth-return geometry.
- `route`: Selected callable, including any bound `FormulaMethod` selectors.
- `formula`: Recipe or frequency functor carrying the selected leaf routes.

# Returns

- The same validated `pair`.
"""
validate(pair::EarthPair, route, formula) = validate(pair)

"""
$(TYPEDSIGNATURES)

Check conductor depths against their resolved horizontal earth layers.

# Arguments

- `pair`: Resolved earth-return geometry.
- `thickness`: Layer thicknesses \\[m\\], including semi-infinite air first.

# Returns

- The same `pair`.

# Errors

- Throws `ArgumentError` when a layer index or conductor depth disagrees with
  the supplied horizontal layer inventory.
"""
function validate(pair::EarthPair, thickness::Union{Tuple, AbstractVector})
    validate(pair)
    for position in eachindex(pair.layers)
        layer = pair.layers[position]
        layer <= length(thickness) || throw(ArgumentError(
            "earth-pair conductor $position refers to absent physical layer $layer"))
        layer == 1 && continue
        depth = -pair.heights[position]
        top = sum((thickness[index] for index in 2:(layer - 1)); init=zero(depth))
        local_depth = depth - top
        local_depth >= zero(depth) &&
            (!isfinite(thickness[layer]) || local_depth <= thickness[layer]) ||
            throw(ArgumentError(
                "earth-pair conductor $position depth $depth m is outside its resolved earth layer $layer"))
    end
    return pair
end

"Tag line parameters expressed in the physical phase domain."
struct PhaseDomain <: LineParamsDomain end
"""
Store the coordinate system of a calculated modal transformation.

The operator tensor type parameterizes the domain because inverse transforms
consume it numerically. The owning transform module may use one formula-family
parameter for provenance storage, so different concrete formula identities can
share one concrete result-space element type.
"""
struct ModalDomain{O, F} <: LineParamsDomain
    "Frequency-dependent phase-to-modal voltage and current operators."
    operators::O
    "Formula that resolved mode order, scaling, and phase convention."
    formula::F
end

"Return a domain value restricted to selected frequency samples."
selectdomain(domain::LineParamsDomain, _) = domain

@inline domain(::Type{PhaseDomain}) = PhaseDomain
@inline domain(::Type{<:ModalDomain}) = ModalDomain

"Return the domain tag type of a value, or `nothing` when it has no domain."
@inline domain(::Type) = nothing
@inline domain(value) = domain(typeof(value))
