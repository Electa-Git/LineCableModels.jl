"""
    _normalize_radii(::Type{Part}, r_in, r_ex)

Validate the strict materialized radial boundary. Both values must be real
numbers. Interpretation of diameter, layer objects, or backward thickness
inference is deliberately unsupported.
"""
@inline function _normalize_radii(::Type{Part}, r_in::Real, r_ex::Real) where {Part}
    return r_in, r_ex
end

@inline function _normalize_radii(::Type{Part}, r_in, r_ex) where {Part}
    throw(ArgumentError(
        "[$(nameof(Part))] radial inputs must be real numeric values; got " *
        "$(typeof(r_in)) and $(typeof(r_ex))",
    ))
end

"""
    _resolve_outer_radius(::Type{Part}, r_in; radius=nothing, thickness=nothing)

Resolve one forward radial declaration. Exactly one of `radius` (an absolute
outer radius) and `thickness` (an outward increment) must be supplied.
"""
function _resolve_outer_radius(
    ::Type{Part},
    r_in::Real;
    radius::Union{Nothing,Real}=nothing,
    thickness::Union{Nothing,Real}=nothing,
) where {Part}
    xor(radius === nothing, thickness === nothing) || throw(ArgumentError(
        "[$(nameof(Part))] provide exactly one of radius or thickness",
    ))
    r_ex = radius === nothing ? r_in + thickness : radius
    r_ex > r_in || throw(ArgumentError(
        "[$(nameof(Part))] outer radius must exceed inner radius; got r_in=$r_in, r_ex=$r_ex",
    ))
    return r_ex
end

@inline function _resolve_group_radial_args(
    ::Type{Part},
    r_in::Real,
    args::Tuple,
    kwargs::NamedTuple,
) where {Part}
    has_radius = haskey(kwargs, :radius)
    has_thickness = haskey(kwargs, :thickness)
    (has_radius || has_thickness) || return args, kwargs
    Validation.has_radii(Part) || throw(ArgumentError(
        "[$(nameof(Part))] radius/thickness keywords are not valid for this part",
    ))
    radius = has_radius ? kwargs.radius : nothing
    thickness = has_thickness ? kwargs.thickness : nothing
    r_ex = _resolve_outer_radius(Part, r_in; radius, thickness)
    cleaned = Base.structdiff(kwargs, (; radius=nothing, thickness=nothing))
    return (r_ex, args...), cleaned
end
