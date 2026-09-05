"Construct one physical pose, or a finite space of poses."
function Pose2(;
        x = 0,
        y = 0,
        φ = 0,
        combine::Symbol = :product
)
    return parameterize(DataModel.Pose2, DataModel.Pose2, (x, y, φ); combine)
end

"Declare a circular primitive, or a finite space of circular primitives."
function Disk(r; combine::Symbol = :product)
    return parameterize(
        DataModel.Disk,
        DataModel.Disk,
        (r,);
        combine
    )
end

"Declare a rectangular primitive, or a finite space of rectangular primitives."
function Rectangle(w, h; combine::Symbol = :product)
    return parameterize(
        DataModel.Rectangle,
        DataModel.Rectangle,
        (w, h);
        combine
    )
end

"Declare an elliptical primitive, or a finite space of elliptical primitives."
function Ellipse(a, b; combine::Symbol = :product)
    return parameterize(
        DataModel.Ellipse,
        DataModel.Ellipse,
        (a, b);
        combine
    )
end

"Declare a filleted cable-sector primitive, or a finite space of sectors."
function Sector(;
        span,
        r_base,
        r_back,
        fillet = 0,
        combine::Symbol = :product
)
    return parameterize(
        DataModel.Sector,
        DataModel.Sector,
        (span, r_base, r_back, fillet);
        combine
    )
end

"Declare an annular primitive, or a finite space of annular primitives."
function Annulus(ri, ro; combine::Symbol = :product)
    return parameterize(
        DataModel.Annulus,
        DataModel.Annulus,
        (ri, ro);
        combine
    )
end

"Declare a contextual shell, or a finite space of contextual shells."
function Shell(t; combine::Symbol = :product)
    return parameterize(
        DataModel.Shell,
        DataModel.Shell,
        (t,);
        combine
    )
end

"Declare a polygon primitive, or a finite space of polygon primitives."
function Polygon(
        points::Union{AbstractGrid, Gridspace};
        combine::Symbol = :product
)
    return parameterize(
        DataModel.Polygon,
        DataModel.Polygon,
        (points,);
        combine
    )
end
