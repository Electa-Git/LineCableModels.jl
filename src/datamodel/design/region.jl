"""
$(TYPEDEF)

Represent one homogeneous physical domain with intrinsic geometry.

$(TYPEDFIELDS)
"""
struct Region{P, M} <: AbstractCablePart
    "Physical identity within its containing cable object."
    tag::Symbol
    "Intrinsic cross-sectional geometry, or `nothing` to take the outline of its context."
    primitive::P
    "Constitutive material data."
    material::M

    function Region(tag::Symbol, primitive::P, material::M) where {
            P <: Union{AbstractPrimitive, Shell, Nothing}, M <: AbstractMaterial
    }
        tag === Symbol("") && throw(ArgumentError("region tag cannot be empty"))
        return new{P, M}(tag, primitive, material)
    end
end

"""
Bind a primitive to a material under a tag, or declare a finite space of such regions.
"""
function Region(
        tag::Union{Symbol, AbstractGrid, Gridspace},
        primitive::Union{AbstractPrimitive, Shell, Nothing, AbstractGrid, Gridspace},
        material::Union{AbstractMaterial, AbstractGrid, Gridspace};
        combine::Symbol = :product
)
    return parameterize(Region, Region, (tag, primitive, material); combine)
end

# A material for the cable part of `role` in a region, restricted to the `allowed` kinds.
function validate(material, ::Type{Region}, role::Symbol, allowed::Tuple)
    throw(ArgumentError(
        "$role material must resolve to AbstractMaterial; got $(nameof(typeof(material)))"
    ))
end
function validate(material::AbstractMaterial, ::Type{Region}, role::Symbol, allowed::Tuple)
    material.kind in allowed || throw(ArgumentError(
        "$role material must have kind $(join(string.(allowed), " or ")); " *
        "got :$(material.kind)"
    ))
    return material
end

function Base.:(==)(left::Region, right::Region)
    left.tag == right.tag &&
        left.primitive == right.primitive &&
        left.material == right.material
end

function Base.isequal(left::Region, right::Region)
    isequal(left.tag, right.tag) &&
        isequal(left.primitive, right.primitive) &&
        isequal(left.material, right.material)
end

function Base.hash(region::Region, h::UInt)
    hash((:Region, region.tag, region.primitive, region.material), h)
end

"""
$(TYPEDEF)

Pair one physical region with its resolved primitive and retained terminal, with resolved declarations for placement and paths.

Each placement retains its physical scope through the existing owner type.
`Assembly` placements preserve independent terminals and are not strand courses.
An enclosing `Group` coalesces those placements into its own terminal scope.
Member indices, placement poses and longitudinal paths remain unchanged.

$(TYPEDFIELDS)
"""
struct PlacedRegion{
    R <: Region,
    S <: AbstractShape,
    P <: NamedTuple{(:patterns,), <:Tuple{Tuple{Vararg{
        NamedTuple{(:owner, :pattern, :member, :pose), <:Tuple{Type, Any, Integer, Pose2}}
    }}}},
    H <: Tuple{Vararg{NamedTuple{(:path, :radius, :axis), <:Tuple{Any, Real, Tuple{Real, Real}}}}}
}
    "Authoritative homogeneous physical region."
    source::R
    "Resolved cross-sectional primitive with its absolute pose."
    primitive::S
    "Retained electrical terminal, or `nothing` for a nonconductive region."
    terminal::Union{Nothing, Symbol}
    "Ordered `(owner, pattern, member, pose)` placements. `Group` coalesces terminals, `Assembly` preserves them."
    placement::P
    "Ordered `(path, radius, axis)` declarations traversed by the region. A helix winds about
    its axis, the center of the part that declared it, at the radius."
    paths::H

    function PlacedRegion(
            source::R,
            primitive::S,
            terminal::Union{Nothing, Symbol},
            placement::P,
            paths::H
    ) where {
            R <: Region,
            S <: AbstractShape,
            P <: NamedTuple{(:patterns,), <:Tuple{Tuple{Vararg{
                NamedTuple{(:owner, :pattern, :member, :pose),
                    <:Tuple{Type, Any, Integer, Pose2}}
            }}}},
            H <: Tuple{Vararg{NamedTuple{(:path, :radius, :axis),
                <:Tuple{Any, Real, Tuple{Real, Real}}}}}
    }
        terminal === nothing || terminal !== Symbol("") ||
            throw(
                ArgumentError("placed-region terminal cannot be empty")
            )
        all(placement.patterns) do entry
            entry.owner <: AbstractCablePart && entry.member > 0
        end || throw(ArgumentError(
            "placed-region patterns must contain a physical owner type, positive member indices and poses"
        ))
        all(paths) do entry
            isfinite(entry.radius) && entry.radius >= 0
        end || throw(ArgumentError(
            "placed-region paths must contain finite nonnegative radii"
        ))
        return new{R, S, P, H}(source, primitive, terminal, placement, paths)
    end
end

function Base.:(==)(left::PlacedRegion, right::PlacedRegion)
    left.source == right.source && left.primitive == right.primitive &&
        left.terminal == right.terminal && left.placement == right.placement &&
        left.paths == right.paths
end

function PlacedRegion(source::Region, primitive::AbstractShape)
    return PlacedRegion(
        source,
        primitive,
        nothing,
        (patterns = (),),
        ()
    )
end

boundary(placed::PlacedRegion) = boundary(placed.primitive)
area(placed::PlacedRegion) = area(placed.primitive)
centroid(placed::PlacedRegion) = centroid(placed.primitive)
support(placed::PlacedRegion, φ::Real) = support(placed.primitive, φ)

# A shell grows outward from the outline of a preceding part. At the cable center there is
# none.
validate(region::Region, ::Union{EmptyBoundary, AbstractShape}) = region
function validate(region::Region{<:Shell}, ::EmptyBoundary)
    throw(ArgumentError(
        "the shell :$(region.tag) has no preceding part to surround; start a tube with " *
        "core(material, Annulus(ri, ro)), or with a filler and then the shell"
    ))
end

# A part declared without its outline takes it from its context. Only the center of a
# Milliken supplies one, the disk that its segments leave.
function validate(region::Region{Nothing}, ::Union{EmptyBoundary, AbstractShape})
    throw(ArgumentError(
        "the region :$(region.tag) has no primitive; a part declared without its outline " *
        "takes it from its context, and only the center of a milliken supplies one"
    ))
end

function resolve(context::Union{EmptyBoundary, AbstractShape}, region::Region)
    primitive = resolve(context, validate(region, context).primitive)
    placed = PlacedRegion(region, primitive)
    return CableGeometry(PlacedRegion[placed], boundary(placed))
end
