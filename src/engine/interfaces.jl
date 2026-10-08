"""
Return the series-impedance values of a line-parameter result.
"""
function Z end

"""
Return the shunt-admittance values of a line-parameter result.
"""
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
"""
Return absolute numerical errors from an owned comparison result.
"""
function absolute_error end
"""
Return reference-normalized numerical errors from an owned comparison result.
"""
function relative_error end

"""
Return electrical conductivity \\[S/m\\] from resistivity \\[Ω·m\\], including
the open-circuit limit.
"""
@inline conductivity(rho) = isinf(rho) ? zero(rho) : inv(rho)

"""
Return a real scalar magnitude for numerical error estimates and physical-state
comparisons. Deterministic values use absolute magnitude. Uncertainty extensions
also account for uncertain contributions with zero nominal value. Evaluated
physical quantities retain their correlations.
"""
@inline numerical_magnitude(z) = abs(complex(nominal(real(z)), nominal(imag(z))))

# Equality of physical state includes correlation, not just nominal values.
same_physical_state(a, b) = isequal(a, b)
same_physical_state(a::BigFloat, b::BigFloat) = isequal(a, b)
function same_physical_state(a::Number, b::Number)
    a === b || (isequal(a, b) && iszero(numerical_magnitude(a-b)))
end
function same_physical_state(a::Tuple, b::Tuple)
    # Tuple map preserves each field's type in heterogeneous physical records.
    length(a) == length(b) && all(map(same_physical_state, a, b))
end
function same_physical_state(a::NamedTuple, b::NamedTuple)
    keys(a) == keys(b) && same_physical_state(values(a), values(b))
end
function same_physical_state(a::AbstractArray, b::AbstractArray)
    axes(a) == axes(b) && all(pair -> same_physical_state(pair...), zip(a, b))
end

"""
Resolve a placed conductor's earth-layer index from the problem and coordinates
[m], or read `(source, target)` indices from an initialized `EarthPair`. Air is
layer 1. Subsequent indices retain the physical earth model's layer order.
"""
function layer_index end

"""
Resolve the real scalar representation needed by an active formulation before
allocating numerical storage. Arguments of active formulations may widen `T`.
Frequency-aligned arguments are
validated by their scientific owner against `frequencies`.
"""
computation_type(::Type{T}, ::AbstractFormulation, frequencies) where {T <: Real} = T
computation_type(::Type{T}, ::Nothing, frequencies) where {T <: Real} = T
# Each selection of a record, and of a recipe within it, can widen `T`.
function computation_type(::Type{T}, selections::NamedTuple, frequencies) where {T <: Real}
    return foldl(values(selections); init = T) do scalar, selected
        computation_type(scalar, selected, frequencies)
    end
end

"""Identify the local formula selections that determine blueprint coefficients."""
function blueprint_dependencies end

"""Calculate the selected local shunt response while constructing a blueprint."""
function internal_shunt_response end
"""
$(TYPEDSIGNATURES)

Calculate the selected earth contributions at one frequency from completed material
inputs: for every calculation of the earth plan, or for one calculation. Given a formula and
the Functor of its calculation, return the physical matrices of the calculation after its
parts, by quantity. A formula whose parts write the physical coefficients returns none. A
coupled formula converts its parts' coefficients for both quantities in one calculation.
Matrix rows are receivers and columns are sources. Impedance contributions are \\[Ω/m\\].
Potential-coefficient contributions are \\[m/F\\].
"""
function earth! end
function homogenize! end

"""
$(TYPEDSIGNATURES)

Evaluate the material quantities required by a coaxial computation into its
allocated buffers. Material-law methods receive the original material values.
Field equations consume these completed quantities without reevaluating them.
"""
function materials! end

"""
Abstract tag for the physical domain represented by line-parameter matrices.
"""
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
    "Source and target heights relative to the air-earth interface, in that order \\[m\\]."
    heights::Tuple{T, T}
    "Horizontal distance between conductor centers \\[m\\]."
    separation::T
    "Physical layer indices of the source and target conductors."
    layers::Tuple{Int, Int}
    "Conductor outer radius for a self interaction [m]. Nothing for a mutual. "
    radius::Union{Nothing, T}
end

function EarthPair(row::Integer, column::Integer, heights::Tuple{T, T}, separation::T,
        layers::Tuple{Int, Int}; radius = nothing) where {T <: Real}
    return EarthPair{T}(row, column, heights, separation, layers, radius)
end

"""
$(TYPEDSIGNATURES)

Return `pair` on the earth that its calculation consumes. The layered earth model keeps the
physical layers.
"""
EarthPair(pair::EarthPair, ::EarthModel) = pair

"""
$(TYPEDSIGNATURES)

Return `pair` on an equivalent homogeneous earth: air stays layer 1, and every earth layer
becomes layer 2.
"""
function EarthPair(pair::EarthPair, ::EquivalentHomogeneous.AbstractSequence)
    return EarthPair(pair.row, pair.column, pair.heights, pair.separation,
        map(layer -> layer == 1 ? 1 : 2, pair.layers); radius = pair.radius)
end

"""
$(TYPEDSIGNATURES)

Check the resolved geometry of an earth-return interaction.

# Arguments

- `pair`: matrix indices and physical layer indices, with signed heights and
  horizontal separation. An explicit self radius completes the geometry.
  Lengths are in \\[m\\]. Layer 1 is air.

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
        iszero(pair.separation) ||
            throw(ArgumentError("a self interaction has zero horizontal separation"))
        pair.radius !== nothing && isfinite(pair.radius) && pair.radius > 0 ||
            throw(DomainError(pair.radius, "a self interaction requires an explicit positive conductor radius"))
    else
        pair.radius === nothing ||
            throw(ArgumentError("a mutual interaction has no self radius"))
        iszero(pair.separation) && pair.heights[1] == pair.heights[2] &&
            throw(DomainError((pair.heights, pair.separation),
                "distinct earth-return conductors cannot have coincident centers"))
    end
    return pair
end

"""
$(TYPEDSIGNATURES)

Check conductor depths against their resolved horizontal earth layers.

# Arguments

- `pair`: resolved earth-return geometry.
- `thickness`: layer thicknesses \\[m\\], including semi-infinite air first.

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
        top = sum((thickness[index] for index in 2:(layer - 1)); init = zero(depth))
        local_depth = depth - top
        local_depth >= zero(depth) &&
        (!isfinite(thickness[layer]) || local_depth <= thickness[layer]) ||
            throw(ArgumentError(
                "earth-pair conductor $position depth $depth m is outside its resolved earth layer $layer"))
    end
    return pair
end

"""
Tag line parameters expressed in the physical phase domain.
"""
struct PhaseDomain <: LineParamsDomain end
"""
Supertype of the modal-to-phase voltage and current bases of a modal transformation.
A concrete subtype implements `Base.size`, which returns the
`(modes, modes, frequencies)` dimensions of its operator tensors.
"""
abstract type AbstractModalOperators end

"""
Store the coordinate system of a calculated modal transformation.

The operator tensor type parameterizes the domain because inverse transforms
consume it numerically. The defining transform module may use one formula-family
parameter to record the selected formula, so different concrete formula identities can
share one concrete element type for the result space.
"""
struct ModalDomain{O, G} <: LineParamsDomain
    "Frequency-dependent modal-to-phase voltage and current bases."
    operators::O
    "Aligned propagation roots in the coefficient basis."
    gamma::G

    function ModalDomain(operators::O, gamma::G) where {O, G}
        return validate(new{O, G}(operators, gamma))
    end
end

function validate(domain::ModalDomain)
    domain.gamma isa AbstractMatrix || throw(DimensionMismatch(
        "modal roots must be a mode×frequency matrix"))
    domain.operators isa AbstractModalOperators || throw(ArgumentError(
        "modal domain requires a validated ModalOperators value"))
    dimensions=size(domain.operators)
    size(domain.gamma)==(dimensions[2], dimensions[3]) || throw(DimensionMismatch(
        "modal roots must align with operator modes and frequency samples"))
    return domain
end

"""
Return a domain value restricted to selected frequency samples.
"""
selectdomain(domain::LineParamsDomain, _) = domain
selectdetails(details, ::LineParamsDomain, _) = details

@inline domain(::Type{PhaseDomain}) = PhaseDomain
@inline domain(::Type{<:ModalDomain}) = ModalDomain

"""
Return the domain tag type of a value, or `nothing` when it has no domain.
"""
@inline domain(::Type) = nothing
@inline domain(value) = domain(typeof(value))
