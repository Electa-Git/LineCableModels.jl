function Region(tag, primitive, material; combine::Symbol = :product)
    values = (tag, primitive, material)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.Region, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Region}(DataModel.Region, grids; combine)
end

function Stack(items::DataModel.AbstractCablePart...; combine::Symbol = :product)
    return DataModel.Stack(DataModel.AbstractCablePart[items...])
end

function Stack(items...; combine::Symbol = :product)
    any(item -> item isa Union{AbstractGrid, Gridspace}, items) ||
        throw(ArgumentError("stack items must resolve to cable parts"))
    grids = map(items) do item
        item isa Union{AbstractGrid, Gridspace} ? item : Grid((item,))
    end
    build = (resolved...) -> DataModel.Stack(DataModel.AbstractCablePart[resolved...])
    return Gridspace{DataModel.Stack}(build, grids; combine)
end

function Group(
        name::Symbol,
        item;
        at = DataModel.Pose2(0, 0, 0),
        pattern = nothing,
        path = nothing,
        compact = nothing,
        combine::Symbol = :product
)
    values = (name, at, item, pattern, path, compact)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Group}(DataModel.Group, grids; combine)
    end
    return DataModel.Group(values...)
end

function Assembly(
        item;
        at = DataModel.Pose2(0, 0, 0),
        pattern = nothing,
        path = nothing,
        compact = nothing,
        names,
        combine::Symbol = :product
)
    values = (at, item, pattern, path, compact, names)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Assembly}(DataModel.Assembly, grids; combine)
    end
    return DataModel.Assembly(values...)
end

function Enclosure(
        tag::Symbol,
        item;
        at = DataModel.Pose2(0, 0, 0),
        primitive,
        fill,
        wall = nothing,
        combine::Symbol = :product
)
    values = (tag, at, primitive, item, fill, wall)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Enclosure}(DataModel.Enclosure, grids; combine)
    end
    return DataModel.Enclosure(values...)
end

function build(
        ::Type{DataModel.CableDesign},
        cable_id,
        root;
        combine::Symbol = :product
)
    values = (cable_id, root)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) || throw(
        MethodError(build, (DataModel.CableDesign, cable_id, root))
    )
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    caller = (identifier, part) -> build(
        DataModel.CableDesign,
        identifier,
        part
    )
    return Gridspace{DataModel.CableDesign}(caller, grids; combine)
end
