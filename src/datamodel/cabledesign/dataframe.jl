import DataFrames: DataFrame

"""
$(TYPEDSIGNATURES)

Return one physical row for every placed region in a completed cable design.

The table reports declarations and resolved geometry. Electrical parameters
belong to computed [`CableConstants`](@ref) or `LineParameters` results.

# Arguments

- `design`: Completed cable design.

# Returns

- A `DataFrame` in physical traversal order.
"""
function DataFrame(design::CableDesign)::DataFrame
    regions = design.geometry.regions
    return DataFrame(
        tag = Symbol[placed.source.tag for placed in regions],
        terminal = Union{Missing, Symbol}[
            placed.terminal === nothing ? missing : placed.terminal
            for placed in regions
        ],
        primitive = Symbol[nameof(typeof(placed.source.primitive)) for placed in regions],
        material_kind = Symbol[placed.source.material.kind for placed in regions],
        area = [area(placed.shape) for placed in regions],
        centroid_x = [first(centroid(placed.shape)) for placed in regions],
        centroid_y = [last(centroid(placed.shape)) for placed in regions],
        overlength = [
            prod(
                entry -> overlength(entry.path, entry.radius),
                placed.paths;
                init = one(eltype(placed.shape))
            )
            for placed in regions
        ]
    )
end
