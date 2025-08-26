"""
$(TYPEDSIGNATURES)

Resolves radius parameters for cable components, converting from various input formats to standardized inner radius, outer radius, and thickness values.

This function serves as a high-level interface to the radius resolution system. It processes inputs through a two-stage pipeline:
1. First normalizes input parameters to consistent forms using [`_parse_inputs_radius`](@ref).
2. Then delegates to specialized implementations via [`_do_resolve_radius`](@ref) based on the component type.

# Arguments

- `param_in`: Inner boundary parameter (defaults to radius) \\[m\\].
  Can be a number, a [`Diameter`](@ref) , a [`Thickness`](@ref), or an [`AbstractCablePart`](@ref).
- `param_ext`: Outer boundary parameter (defaults to radius) \\[m\\].
  Can be a number, a [`Diameter`](@ref) , a [`Thickness`](@ref), or an [`AbstractCablePart`](@ref).
- `object_type`: Type associated to the constructor of the new [`AbstractCablePart`](@ref).

# Returns

- `radius_in`: Normalized inner radius \\[m\\].
- `radius_ext`: Normalized outer radius \\[m\\].
- `thickness`: Computed thickness or specialized dimension depending on the method \\[m\\].
  For [`WireArray`](@ref) components, this value represents the wire radius instead of thickness.

# See also

- [`Diameter`](@ref)
- [`Thickness`](@ref)
- [`AbstractCablePart`](@ref)
"""
@inline _normalize_radii(::Type{T}, rin, rex) where {T} =
  _do_normalize_radii(_parse_inputs_radius(rin, T), _parse_inputs_radius(rex, T), T)

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

$(_CLEANMETHODLIST)

# See also

- [`Diameter`](@ref)
- [`Thickness`](@ref)
- [`strip_uncertainty`](@ref)
"""
function _parse_inputs_radius end

# ------------ Input parsing
@inline _parse_inputs_radius(x::Number, ::Type{T}) where {T} = x
@inline _parse_inputs_radius(d::Diameter, ::Type{T}) where {T} = d.value / 2
@inline _parse_inputs_radius(p::Thickness, ::Type{T}) where {T} = p
@inline function _parse_inputs_radius(p::AbstractCablePart, ::Type{T}) where {T}
  r = getfield(p, :radius_ext)                     # outer radius of prior layer
  return (typeof(p) == T) ? r : strip_uncertainty(r)
end
@inline _parse_inputs_radius(x::AbstractString, ::Type{T}) where {T} =
  throw(ArgumentError("[$(nameof(T))] radius parameter must be numeric, not String: $(repr(x))"))
@inline _parse_inputs_radius(x, ::Type{T}) where {T} =
  throw(ArgumentError("[$(nameof(T))] unsupported radius parameter $(typeof(x)): $(repr(x))"))



# ------------ Input parsing
@inline function _do_normalize_radii(radius_in::Number, radius_ext::Number, ::Type{T}) where {T}
  return radius_in, radius_ext
end

@inline function _do_normalize_radii(radius_in::Number, thickness::Thickness, ::Type{T}) where {T}
  return radius_in, (radius_in + thickness.value)
end

@inline function _do_normalize_radii(radius_in::Number, radius_wire::Number, ::Type{AbstractWireArray})
  return radius_in, radius_in + (2 * radius_wire)
end