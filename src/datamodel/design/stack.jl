"""
$(TYPEDEF)

Store an ordered outward composition of cable parts.

$(TYPEDFIELDS)
"""
struct Stack{V <: AbstractVector{<:AbstractCablePart}} <: AbstractCablePart
    "Physical members in inner-to-outer order."
    items::V

    function Stack(items::V) where {V <: AbstractVector{<:AbstractCablePart}}
        isempty(items) && throw(ArgumentError("a stack requires at least one item"))
        return new{V}(copy(items))
    end
end

Stack(item::AbstractCablePart) = Stack(AbstractCablePart[item])
Stack(items::AbstractCablePart...) = Stack(AbstractCablePart[items...])

Base.:(==)(left::Stack, right::Stack) = left.items == right.items

"""
$(TYPEDEF)

Store the resolved physical regions and outer boundary of one cable root.

$(TYPEDFIELDS)
"""
struct CableGeometry{V <: AbstractVector{<:PlacedRegion}, B <: AbstractShape}
    "Resolved regions in physical traversal order."
    regions::V
    "Resolved outer boundary of the complete root."
    outer::B

    function CableGeometry(regions::V, outer::B) where {
            V <: AbstractVector{<:PlacedRegion},
            B <: AbstractShape
    }
        isempty(regions) && throw(ArgumentError(
            "cable geometry requires at least one placed region"
        ))
        all(regions) do placed
            region_area = area(placed)
            center = centroid(placed)
            isfinite(region_area) && region_area > 0 && all(isfinite, center)
        end || throw(ArgumentError(
            "cable geometry requires finite, positive-area placed regions"
        ))
        extent = support(outer)
        isfinite(extent) && extent > 0 || throw(ArgumentError(
            "cable geometry requires a finite positive outer extent"
        ))
        return new{V, B}(regions, outer)
    end
end

boundary(geometry::CableGeometry) = geometry.outer

Base.:(==)(left::CableGeometry, right::CableGeometry) =
    left.regions == right.regions && left.outer == right.outer

function resolve(context::Union{EmptyBoundary, AbstractShape}, stack::Stack)
    regions = PlacedRegion[]
    state = context
    for item in stack.items
        result = resolve(state, item)
        append!(regions, result.regions)
        state = boundary(result)
    end
    return CableGeometry(regions, state)
end
