using LineCableModels.DataModel: Sector, Ellipse, Polygon, ArcPolygon, PerforatedShape
using LineCableModels.BaseGeometry: Curve, Segment, Arc

# Sample every concrete shape type away from the origin, at rotations that are no
# multiple of π/2. The bar contains two round holes. The course has six wires tangent
# to both circles of an annulus. The shell has three wires tangent to the inner circle
# alone. Every point and center of curvature lies within 0.05 m of the origin.
function shape_gallery()
    bar = Placement(0.005, 0.03, 0.25)
    course = Placement(-0.025, -0.01, 0.35)
    shell = Placement(0.03, 0.01, 1.9)
    half = Arc(Placement(0.025, -0.03, 0.3), 3.0e-3, π)
    return [Circle(2.0e-3, Placement(0.03, -0.02, 0.4)),
        Annulus(1.5e-3, 3.0e-3, Placement(-0.02, 0.01, 1.2)),
        Rectangle(4.0e-3, 2.0e-3, Placement(0.01, 0.025, 0.6)),
        Sector(0, 10.0e-3, 2π / 3, 0, 0, Placement(0.015, -0.01, 0.5)),
        Sector(5.0e-3, 10.0e-3, π / 2, 0, 1.0e-3, Placement(-0.01, -0.02, 2.3)),
        Ellipse(3.0e-3, 2.0e-3, Placement(0.02, 0.02, 0.7)),
        Polygon(
            [(0.0, 0.0), (4.0e-3, 0.0), (5.0e-3, 3.0e-3), (2.0e-3, 5.0e-3),
                (-1.0e-3, 3.0e-3)],
            Placement(-0.03, 0.005, 0.9)),
        ArcPolygon(Curve{Float64}[half, Segment(half(1), half(0))]),
        PerforatedShape(Rectangle(12.0e-3, 6.0e-3, bar),
            [Circle(1.0e-3, bar * Placement(3.0e-3, 0.5e-3, 0.8)),
                Circle(1.5e-3, bar * Placement(-2.5e-3, -0.5e-3, -0.6))]),
        PerforatedShape(Annulus(4.0e-3, 6.0e-3, course),
            [Circle(1.0e-3, course * Placement(0, 0, k * π / 3) *
                            Placement(5.0e-3, 0, 0.2))
             for k in 0:5]),
        PerforatedShape(Annulus(4.0e-3, 7.0e-3, shell),
            [Circle(1.0e-3, shell * Placement(0, 0, 2k * π / 3) *
                            Placement(5.0e-3, 0, 0.2))
             for k in 0:2])]
end
