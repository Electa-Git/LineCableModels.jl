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

function DiskDefinition(r::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.DiskDefinition}(DataModel.DiskDefinition, (r,); combine)
end

function RectangleDefinition(w, h; combine::Symbol = :product)
    values = (w, h)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.RectangleDefinition, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.RectangleDefinition}(DataModel.RectangleDefinition, grids; combine)
end

function EllipseDefinition(a, b; combine::Symbol = :product)
    values = (a, b)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.EllipseDefinition, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.EllipseDefinition}(DataModel.EllipseDefinition, grids; combine)
end

function SectorDefinition(ri, ro, φ0, span; combine::Symbol = :product)
    values = (ri, ro, φ0, span)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.SectorDefinition, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.SectorDefinition}(DataModel.SectorDefinition, grids; combine)
end

function AnnulusDefinition(ri, ro; combine::Symbol = :product)
    values = (ri, ro)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        throw(MethodError(DataModel.AnnulusDefinition, values))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.AnnulusDefinition}(DataModel.AnnulusDefinition, grids; combine)
end

function ShellDefinition(t::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.ShellDefinition}(DataModel.ShellDefinition, (t,); combine)
end

function PolygonDefinition(points::Union{AbstractGrid, Gridspace}; combine::Symbol = :product)
    return Gridspace{DataModel.PolygonDefinition}(DataModel.PolygonDefinition, (points,); combine)
end
