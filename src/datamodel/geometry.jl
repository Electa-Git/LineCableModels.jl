"Return three cable-centre coordinates in a trefoil formation."
function trefoil_formation(x0::Real, y0::Real, radius::Real)
    radius > zero(radius) || throw(DomainError(radius, "radius must be positive"))
    x, y, r = promote(float(x0), float(y0), float(radius))
    distance = r / cosd(30)
    return (
        x,
        y + distance,
        x + distance * cosd(210),
        y + distance * sind(210),
        x + distance * cosd(330),
        y + distance * sind(330)
    )
end

"Return three cable-centre coordinates in a horizontal or vertical flat formation."
function flat_formation(x0::Real, y0::Real, spacing::Real; vertical::Bool = false)
    spacing > zero(spacing) ||
        throw(DomainError(spacing, "spacing must be positive"))
    x, y, distance = promote(float(x0), float(y0), float(spacing))
    return vertical ?
           (x, y, x, y - distance, x, y - 2 * distance) :
           (x, y, x + distance, y, x + 2 * distance, y)
end

"Return the outer radius of a materialised cable design."
function outer_radius(design::CableDesign)
    return support(boundary(design.geometry))
end
