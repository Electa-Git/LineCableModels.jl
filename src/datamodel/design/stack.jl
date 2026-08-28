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

"Store aligned resolved regions, their outer boundary, and terminal order."
struct ResolvedPart{V <: AbstractVector, B, N <: AbstractVector{Symbol}}
    regions::V
    outer::B
    terminals::N
end

boundary(part::ResolvedPart) = part.outer

function resolve(context::Union{EmptyBoundary, AbstractShape}, stack::Stack)
    regions = ResolvedRegion[]
    terminals = Symbol[]
    state = context
    for item in stack.items
        result = resolve(state, item)
        append!(regions, result.regions)
        for terminal in result.terminals
            terminal in terminals || push!(terminals, terminal)
        end
        state = boundary(result)
    end
    return ResolvedPart(regions, state, terminals)
end
