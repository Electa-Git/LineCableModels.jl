"Report that a cable-part type has no packing-limit implementation."
@noinline function maxfill(::Type{T}, args...) where {T}
    throw(ArgumentError(
        "maxfill is not implemented for $(nameof(T)) with geometry " *
        "$(map(typeof, args))",
    ))
end

"Declare a wire-count field constrained by owner-local packing geometry."
struct PhysicalFillLimit <: Validation.Rule
    count::Symbol
    geometry::Tuple{Vararg{Symbol}}
end

function Validation.apply(rule::PhysicalFillLimit, value, owner::Type)
    count = getproperty(value, rule.count)
    count isa Integer || return nothing
    geometry = map(field -> getproperty(value, field), rule.geometry)
    all(candidate -> candidate isa Real, geometry) || return nothing
    limit = maxfill(owner, geometry...)
    count <= limit || throw(DomainError(
        count,
        "[$(nameof(owner))] $(rule.count) exceeds packing limit $limit"
    ))
    return nothing
end
