function Pose2(;
        x = 0,
        y = 0,
        φ = 0,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    values = (x, y, φ)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{DataModel.Pose2}(DataModel.Pose2, grids; combine)
    end
    return DataModel.Pose2(values...)
end

function Disk(r::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.Disk}(DataModel.Disk, (r,); combine)
end

function Rectangle(w, h; combine::Symbol = :product)
    values = (w, h)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.Rectangle, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Rectangle}(DataModel.Rectangle, grids; combine)
end

function Ellipse(a, b; combine::Symbol = :product)
    values = (a, b)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.Ellipse, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Ellipse}(DataModel.Ellipse, grids; combine)
end

function Sector(ri, ro, φ0, span; combine::Symbol = :product)
    values = (ri, ro, φ0, span)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.Sector, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Sector}(DataModel.Sector, grids; combine)
end

function Annulus(ri, ro; combine::Symbol = :product)
    values = (ri, ro)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.Annulus, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Annulus}(DataModel.Annulus, grids; combine)
end

function Shell(t::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.Shell}(DataModel.Shell, (t,); combine)
end

function Polygon(points::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.Polygon}(DataModel.Polygon, (points,); combine)
end
