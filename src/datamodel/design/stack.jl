"""
$(TYPEDEF)

Store an ordered outward composition of cable parts.

$(TYPEDFIELDS)
"""
struct Stack{V <: AbstractVector} <: AbstractCablePart
    "Physical members in inner-to-outer order."
    items::V

    function Stack(items::V) where {V <: AbstractVector}
        isempty(items) && throw(ArgumentError("a stack requires at least one item"))
        all(item -> item isa AbstractCablePart, items) ||
            throw(ArgumentError("stack items must be cable parts"))
        return new{V}(copy(items))
    end
end

Stack(item::AbstractCablePart) = Stack(AbstractCablePart[item])

"""
$(TYPEDEF)

Store the resolved physical regions and outer boundary of one cable root.

$(TYPEDFIELDS)
"""
struct CableGeometry{V <: AbstractVector, B <: AbstractPrimitive}
    "Resolved regions in physical traversal order."
    regions::V
    "Resolved outer boundary of the complete root."
    outer::B

    function CableGeometry(regions::V, outer::B) where {
            V <: AbstractVector,
            B <: AbstractPrimitive
    }
        isempty(regions) && throw(ArgumentError(
            "cable geometry requires at least one placed region"
        ))
        all(region -> region isa PlacedRegion, regions) || throw(ArgumentError(
            "cable geometry regions must be PlacedRegion values"
        ))
        all(regions) do region
            region_area = area(region)
            centre = centroid(region)
            region_area isa Real && isfinite(region_area) && region_area > 0 &&
                all(isfinite, centre)
        end || throw(ArgumentError(
            "cable geometry requires finite, positive-area placed regions"
        ))
        extent = support(outer)
        extent isa Real && isfinite(extent) && extent > 0 || throw(ArgumentError(
            "cable geometry requires a finite positive outer extent"
        ))
        return new{V, B}(regions, outer)
    end
end

boundary(geometry::CableGeometry) = geometry.outer

function resolve(context::Union{EmptyBoundary, AbstractPrimitive}, stack::Stack)
    regions = PlacedRegion[]
    state = context
    for item in stack.items
        result = resolve(state, item)
        append!(regions, result.regions)
        state = boundary(result)
    end
    return CableGeometry(regions, state)
end
