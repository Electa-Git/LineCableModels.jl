"""
$(TYPEDSIGNATURES)

Resolves radius parameters for cable components, converting from various input formats to standardized inner radius, outer radius, and thickness values.

This function serves as a high-level interface to the radius resolution system. It processes inputs through a two-stage pipeline:
1. First normalizes input parameters to consistent forms using [`_parse_radius_operand`](@ref).
2. Then delegates to specialized `_do_resolve_radius` implementations based on the component type.

# Arguments

- `param_in`: Inner boundary parameter (defaults to radius) \\[m\\].
  Can be a number, a [`Diameter`](@ref) , a [`Thickness`](@ref), or an [`AbstractCablePart`](@ref).
- `param_ext`: Outer boundary parameter (defaults to radius) \\[m\\].
  Can be a number, a [`Diameter`](@ref) , a [`Thickness`](@ref), or an [`AbstractCablePart`](@ref).
- `object_type`: Type associated to the constructor of the new [`AbstractCablePart`](@ref).

# Returns

- `r_in`: Normalized inner radius \\[m\\].
- `r_ex`: Normalized outer radius \\[m\\].
- `thickness`: Computed thickness or specialized dimension depending on the method \\[m\\].
  For [`CircStrands`](@ref) components, this value represents the wire radius instead of thickness.

"""
@inline _normalize_radii(::Type{T}, rin, rex) where {T} = _do_normalize_radii(
    _parse_radius_operand(rin, T), _parse_radius_operand(rex, T), T)

"""
$(TYPEDSIGNATURES)

Parses input values into radius representation based on object type and input type.

# Arguments

- `x`: Input value that can be a raw number, a [`Diameter`](@ref), a [`Thickness`](@ref), or other convertible type \\[m\\].
- `object_type`: Type parameter used for dispatch.

# Returns

- Parsed radius value in appropriate units \\[m\\].

# Examples

```julia
radius = $(FUNCTIONNAME)(10.0, ...)   # Direct radius value
radius = $(FUNCTIONNAME)(Diameter(20.0), ...)  # From diameter object
radius = $(FUNCTIONNAME)(Thickness(5.0), ...)  # From thickness object
```

# Methods

$(METHODLIST)

"""
function _parse_radius_operand end

# ------------ Input parsing
@inline _parse_radius_operand(x::Number, ::Type{T}) where {T} = x
@inline _parse_radius_operand(d::Diameter, ::Type{T}) where {T} = d.value / 2
@inline _parse_radius_operand(p::Thickness, ::Type{T}) where {T} = p
# A part proxy denotes the exact physical interface shared by two adjacent
# layers.  Preserve the complete Measurement derivative graph of that boundary;
# stripping it here makes the same radius statistically different on each side
# of the interface and breaks cumulative-geometry covariance.
@inline _parse_radius_operand(p::AbstractCablePart, ::Type{T}) where {T} = getproperty(p, :r_ex)
@inline _parse_radius_operand(x::AbstractString, ::Type{T}) where {T} = throw(
    ArgumentError(
    "[$(nameof(T))] radius parameter must be numeric, not String: $(repr(x))",
),
)
@inline _parse_radius_operand(x, ::Type{T}) where {T} = throw(
    ArgumentError(
    "[$(nameof(T))] unsupported radius parameter $(typeof(x)): $(repr(x))",
),
)

# ------------ Input parsing
@inline function _do_normalize_radii(
        r_in::Number,
        r_ex::Number,
        ::Type{T}
) where {T}
    return r_in, r_ex
end

@inline function _do_normalize_radii(
        r_in::Number,
        thickness::Thickness,
        ::Type{T}
) where {T}
    return r_in, (r_in + thickness.value)
end

@inline function _do_normalize_radii(
        r_in::Number,
        radius_wire::Number,
        ::Type{AbstractStrandsLayer}
)
    return r_in, r_in + (2 * radius_wire)
end

# These intersections make the intended `Thickness` semantics explicit for
# strand layers and keep dispatch unambiguous.
@inline _do_normalize_radii(r_in::Number, thickness::Thickness,
    ::Type{AbstractStrandsLayer}) = (r_in, r_in + thickness.value)

@inline function _do_normalize_radii(
        thickness::Thickness,
        r_ex::Number,
        ::Type{AbstractStrandsLayer}
)
    r_in = r_ex - thickness.value
    r_in >= 0 || throw(
        ArgumentError(
        "[AbstractStrandsLayer] thickness $(thickness.value) exceeds outer radius $(r_ex).",
    ),
    )
    return r_in, r_ex
end

@inline function _do_normalize_radii(
        ::Thickness,
        ::Thickness,
        ::Type{AbstractStrandsLayer}
)
    throw(
        ArgumentError(
        "[AbstractStrandsLayer] cannot specify thickness for both inner and outer radii.",
    ),
    )
end

@inline function _do_normalize_radii(t::Thickness, rex::Number, ::Type{T}) where {T}
    rin = rex - t.value
    rin >= 0 || throw(
        ArgumentError("[$(nameof(T))] thickness $(t.value) exceeds outer radius $(rex)."),
    )
    return rin, rex
end

# NEW: reject thickness on BOTH ends
@inline function _do_normalize_radii(::Thickness, ::Thickness, ::Type{T}) where {T}
    throw(
        ArgumentError(
        "[$(nameof(T))] cannot specify thickness for both inner and outer radii.",
    ),
    )
end
