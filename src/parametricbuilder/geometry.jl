"Construct one physical pose, or a finite space of poses."
function Pose2(;
        x = 0,
        y = 0,
        φ = 0,
        combine::Symbol = :product
)
    return _construction(DataModel.Pose2, DataModel.Pose2, (x, y, φ); combine)
end

"Declare a circular primitive, or a finite space of circular definitions."
function DiskDefinition(r; combine::Symbol = :product)
    return _construction(
        DataModel.DiskDefinition,
        DataModel.DiskDefinition,
        (r,);
        combine
    )
end

"Declare a rectangular primitive, or a finite space of rectangle definitions."
function RectangleDefinition(w, h; combine::Symbol = :product)
    return _construction(
        DataModel.RectangleDefinition,
        DataModel.RectangleDefinition,
        (w, h);
        combine
    )
end

"Declare an elliptical primitive, or a finite space of ellipse definitions."
function EllipseDefinition(a, b; combine::Symbol = :product)
    return _construction(
        DataModel.EllipseDefinition,
        DataModel.EllipseDefinition,
        (a, b);
        combine
    )
end

"Declare a sector primitive, or a finite space of sector definitions."
function SectorDefinition(ri, ro, φ0, span; combine::Symbol = :product)
    return _construction(
        DataModel.SectorDefinition,
        DataModel.SectorDefinition,
        (ri, ro, φ0, span);
        combine
    )
end

"Declare an annular primitive, or a finite space of annulus definitions."
function AnnulusDefinition(ri, ro; combine::Symbol = :product)
    return _construction(
        DataModel.AnnulusDefinition,
        DataModel.AnnulusDefinition,
        (ri, ro);
        combine
    )
end

"Declare a contextual shell, or a finite space of shell definitions."
function ShellDefinition(t; combine::Symbol = :product)
    return _construction(
        DataModel.ShellDefinition,
        DataModel.ShellDefinition,
        (t,);
        combine
    )
end

"Declare a polygon primitive, or a finite space of polygon definitions."
function PolygonDefinition(
        points::Union{AbstractGrid, Gridspace};
        combine::Symbol = :product
)
    return _construction(
        DataModel.PolygonDefinition,
        DataModel.PolygonDefinition,
        (points,);
        combine
    )
end
