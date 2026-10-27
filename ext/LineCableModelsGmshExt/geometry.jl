"""
$(TYPEDEF)

Store the shared points, curves, loops, and surfaces of a cable cross-section in the
caller's Gmsh session.

$(TYPEDFIELDS)
"""
struct Geometry{M <: FieldModel}
    "Model that supplies the regions, materials, placements, and point scale."
    model::M
    "Registered coordinate pairs \\[m\\] and their Gmsh point tags."
    points::Dict{NTuple{2, Float64}, Int}
    "Gmsh point tags and their registered coordinates \\[m\\]."
    positions::Dict{Int, NTuple{2, Float64}}
    "Coordinates in \\[m\\], grouped in cells whose side has length `model.scale` \\[m\\]."
    buckets::Dict{NTuple{2, Int}, Vector{NTuple{2, Float64}}}
    "Sorted endpoint tags and their Gmsh line tags."
    lines::Dict{NTuple{2, Int}, Int}
    "Circle centers and radii \\[m\\], sorted endpoint tags, and Gmsh arc tags."
    arcs::Dict{Tuple{NTuple{3, Float64}, NTuple{2, Int}}, Int}
    "First and last point tags of each stored curve."
    ends::Dict{Int, NTuple{2, Int}}
    "CAD curve of each Gmsh curve: coordinates and radii \\[m\\], angles \\[rad\\]."
    curves::Dict{Int, Curve{Float64}}
    "Angles in \\[rad\\] for cuts on each circle, keyed by center and radius in \\[m\\]."
    breaks::Dict{NTuple{3, Float64}, Vector{Float64}}
    "Coordinates in \\[m\\] of contacts, keyed by circle and angle in \\[rad\\]."
    contacts::Dict{Tuple{NTuple{3, Float64}, Float64}, NTuple{2, Float64}}
    "Counterclockwise and clockwise loop tags and oriented curves, keyed by unsigned tags."
    loops::Dict{Vector{Int}, @NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}}
    "Plane surface tags of each region in system order."
    surfaces::Vector{Vector{Int}}
    "Annular patch tags, region indices, and angular fractions \\[dimensionless\\]."
    patches::Vector{@NamedTuple{region::Int, surface::Int, arcs::NTuple{2, Int},
        spokes::NTuple{2, Int}, corners::NTuple{4, Int}, fraction::Float64}}
    "Counterclockwise and clockwise outline loops and oriented curves of each cable."
    outlines::Vector{@NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}}
    "Clockwise curve lists of the holes in each medium's main face."
    holes::@NamedTuple{air::Vector{Vector{Int}}, earth::Vector{Vector{Int}}}
    "Plane surface tags of the enclosed pockets of each medium."
    pockets::@NamedTuple{air::Vector{Int}, earth::Vector{Int}}
    "Terminal anchor coordinates and metal dimensions \\[m\\]."
    anchors::Vector{@NamedTuple{x::Float64, y::Float64, dimension::Float64}}
end

"""
$(TYPEDSIGNATURES)

Construct empty geometric registries for `model`. Allocate an empty vector of plane surfaces
for each region.
"""
function Geometry(model::FieldModel)
    return Geometry(model,
        Dict{NTuple{2, Float64}, Int}(),
        Dict{Int, NTuple{2, Float64}}(),
        Dict{NTuple{2, Int}, Vector{NTuple{2, Float64}}}(),
        Dict{NTuple{2, Int}, Int}(),
        Dict{Tuple{NTuple{3, Float64}, NTuple{2, Int}}, Int}(),
        Dict{Int, NTuple{2, Int}}(),
        Dict{Int, Curve{Float64}}(),
        Dict{NTuple{3, Float64}, Vector{Float64}}(),
        Dict{Tuple{NTuple{3, Float64}, Float64}, NTuple{2, Float64}}(),
        Dict{Vector{Int}, @NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}}(),
        [Int[] for _ in model.problem.system.regions],
        @NamedTuple{region::Int, surface::Int, arcs::NTuple{2, Int},
            spokes::NTuple{2, Int}, corners::NTuple{4, Int}, fraction::Float64}[],
        @NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}[],
        (air = Vector{Int}[], earth = Vector{Int}[]),
        (air = Int[], earth = Int[]),
        @NamedTuple{x::Float64, y::Float64, dimension::Float64}[])
end

"""
$(TYPEDSIGNATURES)

Return the defining `numbers`, converted to `Float64` and rounded to 15 significant digits.
Coordinates and radii are in \\[m\\].
"""
function key(numbers::Real...)
    return map(number -> round(Float64(number); sigdigits = 15), numbers)
end

"""
$(TYPEDSIGNATURES)

Return the registered coordinate pair within 64 machine spacings of `xy` \\[m\\], or its
rounded key. The spacing uses the largest absolute coordinate, with a minimum scale of 1 m.
"""
function key(geometry::Geometry, xy::Tuple{Real, Real})
    pair = key(xy...)
    haskey(geometry.points, pair) && return pair
    tolerance = 64eps(max(abs(Float64(xy[1])), abs(Float64(xy[2])), 1.0))
    scale = geometry.model.scale
    for x in
        floor(Int, (pair[1] - tolerance) / scale):floor(Int, (pair[1] + tolerance) / scale),
        y in
        floor(Int, (pair[2] - tolerance) / scale):floor(Int, (pair[2] + tolerance) / scale)

        for stored in get(geometry.buckets, (x, y), ())
            abs(stored[1] - pair[1]) <= tolerance &&
                abs(stored[2] - pair[2]) <= tolerance && return stored
        end
    end
    return pair
end

"""
$(TYPEDSIGNATURES)

Return `a` \\[rad\\] reduced to `[0, 2π)`, matched to a break of `circle` within
1e-12 rad. The circle's center and radius are in \\[m\\].
"""
function key(geometry::Geometry, circle::Circle, a::Real)
    normalized = mod(Float64(a), 2π)
    isapprox(normalized, 2π; rtol = 0, atol = 128eps(Float64)) && (normalized = 0.0)
    for stored in get(geometry.breaks, key(circle.at.x, circle.at.y, circle.r), Float64[])
        separation = abs(stored - normalized)
        min(separation, 2π - separation) <= 1e-12 && return stored
    end
    return normalized
end

"""
$(TYPEDSIGNATURES)

Return the Gmsh point tag at `xy` \\[m\\]. Reuse a registered coordinate pair when it
matches. New points use `geometry.model.scale` \\[m\\] as their provisional mesh size.
"""
function point!(geometry::Geometry, xy::Tuple{Real, Real})
    pair = key(geometry, xy)
    return get!(geometry.points, pair) do
        scale = geometry.model.scale
        cell = (floor(Int, pair[1] / scale), floor(Int, pair[2] / scale))
        push!(get!(Vector{NTuple{2, Float64}}, geometry.buckets, cell), pair)
        tag = gmsh.model.geo.addPoint(pair[1], pair[2], 0.0, Float64(scale))
        geometry.positions[tag] = pair
        return tag
    end
end

"""
$(TYPEDSIGNATURES)

Draw `segment`, whose endpoints are in \\[m\\], split at registered points on it.
Return oriented line tags. Reuse each sorted endpoint pair.
"""
function draw(geometry::Geometry, segment::Segment)
    a, b = segment(0), segment(1)
    p, q = point!(geometry, a), point!(geometry, b)
    dx, dy = b[1] - a[1], b[2] - a[2]
    distance = hypot(dx, dy)
    tolerance = 64eps(max(abs(a[1]), abs(a[2]), abs(b[1]), abs(b[2]), 1.0))
    ordered = Tuple{Float64, Int}[(0.0, p), (1.0, q)]
    scale = geometry.model.scale
    x_bins = range(floor(Int, (min(a[1], b[1]) - tolerance) / scale),
        floor(Int, (max(a[1], b[1]) + tolerance) / scale))
    y_bins = range(floor(Int, (min(a[2], b[2]) - tolerance) / scale),
        floor(Int, (max(a[2], b[2]) + tolerance) / scale))
    bins = length(x_bins) <= length(geometry.buckets) ÷ length(y_bins) ?
           Iterators.product(x_bins, y_bins) : keys(geometry.buckets)
    for (x, y) in bins
        x in x_bins && y in y_bins || continue
        for xy in get(geometry.buckets, (x, y), ())
            px, py = xy[1] - a[1], xy[2] - a[2]
            t = (px * dx + py * dy) / distance^2
            tolerance / distance < t < 1 - tolerance / distance || continue
            abs(px * dy - py * dx) <= tolerance * distance || continue
            tag = geometry.points[xy]
            tag in (p, q) && continue
            push!(ordered, (t, tag))
        end
    end
    sort!(ordered)
    tags = Int[]
    for i in 1:(length(ordered) - 1)
        start_tag, stop_tag = ordered[i][2], ordered[i + 1][2]
        pair = minmax(start_tag, stop_tag)
        tag = get!(geometry.lines, pair) do
            native = gmsh.model.geo.addLine(pair...)
            geometry.ends[native] = pair
            geometry.curves[native] = Segment(geometry.positions[pair[1]],
                geometry.positions[pair[2]])
            return native
        end
        push!(tags, start_tag == geometry.ends[tag][1] ? tag : -tag)
    end
    return tags
end

"""
$(TYPEDSIGNATURES)

Draw `arc`, whose center and radius are in \\[m\\] and signed span in \\[rad\\],
split at its circle's breaks. Create the center before the points along the arc in its
sense. Return oriented arc tags. Arc ends use exact registered contacts when present.
"""
function draw(geometry::Geometry, arc::Arc)
    cx, cy, radius = Float64(arc.at.x), Float64(arc.at.y), Float64(arc.r)
    circle = Circle(radius, Placement(cx, cy))
    cut!(geometry, circle)
    cut!(geometry, circle, arc.at.φ)
    cut!(geometry, circle, arc.at.φ + arc.span)
    circle_key = key(cx, cy, radius)
    sense = sign(Float64(arc.span))
    distance = abs(Float64(arc.span))
    distances = Float64[0.0, distance]
    for a in geometry.breaks[circle_key]
        offset = mod(sense * (a - Float64(arc.at.φ)), 2π)
        if 1e-12 < offset < distance - 1e-12 &&
           all(existing -> abs(existing - offset) > 1e-12, distances)
            push!(distances, offset)
        end
    end
    sort!(distances)
    angles = Float64(arc.at.φ) .+ sense .* distances
    center_tag = point!(geometry, (cx, cy))
    point_tags = Int[]
    for a in angles
        cosine, sine = cos(a), sin(a)
        isapprox(cosine, 0; rtol = 0, atol = 128eps(Float64)) && (cosine = 0.0)
        isapprox(sine, 0; rtol = 0, atol = 128eps(Float64)) && (sine = 0.0)
        xy = get(geometry.contacts, (circle_key, key(geometry, circle, a)),
            (cx + radius * cosine, cy + radius * sine))
        push!(point_tags, point!(geometry, xy))
    end
    tags = Int[]
    for i in 1:(length(point_tags) - 1)
        p, q = point_tags[i], point_tags[i + 1]
        p == q && continue
        tag = get!(geometry.arcs, (circle_key, minmax(p, q))) do
            native = gmsh.model.geo.addCircleArc(p, center_tag, q)
            geometry.ends[native] = (p, q)
            geometry.curves[native] = Arc(Placement(cx, cy, angles[i]), radius,
                angles[i + 1] - angles[i])
            return native
        end
        push!(tags, (p, q) == geometry.ends[tag] ? tag : -tag)
    end
    return tags
end

"""
$(TYPEDSIGNATURES)

Register the corners of `path` in path order, draw its curves, and return its
counterclockwise and clockwise Gmsh loops and oriented curves. Coordinates are in
\\[m\\] and angles in \\[rad\\]. Reuse the loop with the same unsigned curve tags.
"""
function draw(geometry::Geometry, path::AbstractVector{<:Curve})
    corners = [curve(0)
               for (previous, curve) in zip(circshift(path, 1), path)
               if previous isa Segment || curve isa Segment]
    for xy in corners
        point!(geometry, xy)
    end
    tags = reduce(vcat, [draw(geometry, curve) for curve in path]; init = Int[])
    return get!(geometry.loops, sort(abs.(tags))) do
        ccw = gmsh.model.geo.addCurveLoop(tags)
        cw = gmsh.model.geo.addCurveLoop(-reverse(tags))
        return (ccw = ccw, cw = cw, curves = tags)
    end
end

"""
$(TYPEDSIGNATURES)

Register the break `a` \\[rad\\] on `circle`, whose center and radius are in \\[m\\].
Match a break within 1e-12 rad and return its angle \\[rad\\].
"""
function cut!(geometry::Geometry, circle::Circle, a::Real)
    circle_key = key(circle.at.x, circle.at.y, circle.r)
    angles = get!(Vector{Float64}, geometry.breaks, circle_key)
    stored = key(geometry, circle, a)
    stored in angles || push!(angles, stored)
    return stored
end

"""
$(TYPEDSIGNATURES)

Register `a` \\[rad\\] on `circle`. Keep `xy` in \\[m\\] as the exact contact.
Keep the first contact recorded at the matched break. Return its angle \\[rad\\].
"""
function cut!(geometry::Geometry, circle::Circle, a::Real, xy::Tuple{Real, Real})
    stored = cut!(geometry, circle, a)
    get!(geometry.contacts, (key(circle.at.x, circle.at.y, circle.r), stored),
        (Float64(xy[1]), Float64(xy[2])))
    return stored
end

"""
$(TYPEDSIGNATURES)

Cut `circle`, whose center and radius are in \\[m\\], at its four quadrant angles.
Return `geometry`.
"""
function cut!(geometry::Geometry, circle::Circle)
    for a in (0.0, π / 2, π, 3π / 2)
        cut!(geometry, circle, a)
    end
    return geometry
end

"""
$(TYPEDSIGNATURES)

Cut the circle of `arc` at its quadrants and at the arc's endpoint angles \\[rad\\].
The center and radius are in \\[m\\]. Return `geometry`.
"""
function cut!(geometry::Geometry, arc::Arc)
    circle = Circle(arc.r, arc.at)
    cut!(geometry, circle)
    cut!(geometry, circle, arc.at.φ)
    cut!(geometry, circle, arc.at.φ + arc.span)
    return geometry
end

"""
$(TYPEDSIGNATURES)

Return `geometry` unchanged: a segment lies on no circle.
"""
function cut!(geometry::Geometry, ::Segment)
    return geometry
end

"""
$(TYPEDSIGNATURES)

Cut the circles of every curve in `contour`, whose coordinates and radii are in \\[m\\]
and angles in \\[rad\\]. Return `geometry`.
"""
function cut!(geometry::Geometry, contour::Contour)
    for path in (contour.outer, contour.holes...), curve in path

        cut!(geometry, curve)
    end
    return geometry
end

"""
$(TYPEDSIGNATURES)

Cut the circles of `shape` through its contour. Dimensions and coordinates are in
\\[m\\] and angles in \\[rad\\]. Return `geometry`.
"""
function cut!(geometry::Geometry, shape::AbstractShape)
    return cut!(geometry, Contour(shape))
end

"""
$(TYPEDSIGNATURES)

Cut the contour of `shape` and the contacts of its course of wires. Dimensions and
coordinates are in \\[m\\] and angles in \\[rad\\]. Contact points use the plan's
radial placements exactly. Return `geometry`.
"""
function cut!(geometry::Geometry, shape::DifferenceShape)
    cut!(geometry, Contour(shape))
    plan = compartments(shape)
    isempty(plan.holes) && return geometry
    cx, cy = plan.center
    inner = Circle(plan.inner_radius, Placement(cx, cy))
    outer = Circle(plan.outer_radius, Placement(cx, cy))
    for (a, wire) in zip(plan.angles, plan.holes)
        outer_xy = Placement(cx, cy, a)((plan.outer_radius, 0.0))
        inner_xy = Placement(cx, cy, a)((plan.inner_radius, 0.0))
        cut!(geometry, outer, a, outer_xy)
        cut!(geometry, inner, a, inner_xy)
        cut!(geometry, wire, a, outer_xy)
        cut!(geometry, wire, a + π, inner_xy)
        cut!(geometry, wire, a + π / 2)
        cut!(geometry, wire, a - π / 2)
    end
    return geometry
end

"""
$(TYPEDSIGNATURES)

Cut each pair of tangent or intersecting registered circles at their shared contacts.
Centers, radii, and contact coordinates are in \\[m\\]. Return `geometry`.
"""
function cut!(geometry::Geometry)
    circles = collect(keys(geometry.breaks))
    for i in eachindex(circles)
        cx, cy, r = circles[i]
        for j in (i + 1):length(circles)
            qx, qy, s = circles[j]
            dx, dy = qx - cx, qy - cy
            distance = hypot(dx, dy)
            scale = max(r, s, distance, 1.0)
            tolerance = 1e-12 * scale
            distance <= tolerance && continue
            distance > r + s + tolerance && continue
            distance < abs(r - s) - tolerance && continue
            along = (r^2 - s^2 + distance^2) / (2distance)
            height_squared = r^2 - along^2
            height_squared < -tolerance * scale && continue
            tangent = isapprox(distance, r + s; rtol = 0, atol = tolerance) ||
                      isapprox(distance, abs(r - s); rtol = 0, atol = tolerance)
            height = tangent ? 0.0 : sqrt(max(height_squared, 0.0))
            ux, uy = dx / distance, dy / distance
            bx, by = cx + along * ux, cy + along * uy
            contacts = height <= tolerance ? ((bx, by),) :
                       ((bx - height * uy, by + height * ux),
                (bx + height * uy, by - height * ux))
            p = Circle(r, Placement(cx, cy))
            q = Circle(s, Placement(qx, qy))
            for xy in contacts
                cut!(geometry, p, atan(xy[2] - cy, xy[1] - cx), xy)
                cut!(geometry, q, atan(xy[2] - qy, xy[1] - qx), xy)
            end
        end
    end
    return geometry
end

"""
$(TYPEDSIGNATURES)

Align the breaks of concentric `annuli`. Share their angles in \\[rad\\].
Centers and radii are in \\[m\\]. Return `geometry`.
"""
function cut!(geometry::Geometry, annuli::AbstractVector{<:Annulus})
    groups = Dict{NTuple{2, Float64}, Set{Float64}}()
    for annulus in annuli
        center = key(annulus.at.x, annulus.at.y)
        radii = get!(Set{Float64}, groups, center)
        union!(radii, key(annulus.ri, annulus.ro))
    end
    for (center, radii) in groups
        angles = sort!(unique(reduce(vcat,
            [geometry.breaks[(center..., r)] for r in radii]; init = Float64[])))
        for r in radii, a in angles

            cut!(geometry, Circle(r, Placement(center...)), a)
        end
    end
    return geometry
end

"""
$(TYPEDSIGNATURES)

Return the course plan of `shape`, whose dimensions and coordinates are in \\[m\\]
and angles in \\[rad\\]. A plan with no wires denotes a fill without a tangent course.
"""
function compartments(shape::DifferenceShape)
    return compartments(shape, shape.outer)
end

"""
$(TYPEDSIGNATURES)

Return a plan with no wires for a `shape` whose `outer` is not an annulus. Plan radii and
coordinates are in \\[m\\], and angles are in \\[rad\\].
"""
function compartments(::DifferenceShape, ::AbstractShape)
    return (center = (0.0, 0.0), inner_radius = 0.0, outer_radius = 0.0,
        shell_outer_radius = 0.0, angles = Float64[], holes = Circle{Float64}[],
        outer_spans = NTuple{2, Float64}[])
end

"""
$(TYPEDSIGNATURES)

Return the course of circular wires tangent to the inner circle of `outer`, or a plan
with no wires. Record the center and inner, contact, and shell radii \\[m\\], the wires
and their angles \\[rad\\], and the shell's free angular spans \\[rad\\]. Convert
the cutouts of `shape` to `Float64` and calculate their free spans.
"""
function compartments(shape::DifferenceShape, outer::Annulus)
    plan = invoke(compartments, Tuple{DifferenceShape, AbstractShape}, shape, outer)
    isempty(shape.holes) && return plan
    iszero(outer.ri) && return plan
    center = (Float64(outer.at.x), Float64(outer.at.y))
    wires = [hole for hole in shape.holes if hole isa Circle]
    cutouts = [hole for hole in shape.holes if hole isa AnnularSector]
    length(wires) + length(cutouts) == length(shape.holes) || return plan
    isempty(wires) && return plan
    entries = Tuple{Float64, Circle{Float64}}[]
    outer_contacts = Float64[]
    for wire in wires
        dx, dy = Float64(wire.at.x) - center[1], Float64(wire.at.y) - center[2]
        distance = hypot(dx, dy)
        radius = Float64(wire.r)
        tolerance = 1e-12 * max(Float64(outer.ro), distance, radius, 1.0)
        distance > tolerance || return plan
        isapprox(distance - radius, Float64(outer.ri); rtol = 0, atol = tolerance) ||
            return plan
        push!(outer_contacts, distance + radius)
        a = mod(atan(dy, dx), 2π)
        isapprox(a, 2π; rtol = 0, atol = 128eps(Float64)) && (a = 0.0)
        push!(entries, (a, Circle(radius, convert(Placement{Float64}, wire.at))))
    end
    computed_radius = first(outer_contacts)
    tolerance = 1e-12 * max(Float64(outer.ro), computed_radius, 1.0)
    for r in outer_contacts
        isapprox(r, computed_radius; rtol = 0, atol = tolerance) || return plan
    end
    computed_radius <= Float64(outer.ro) + tolerance || return plan
    contact_radius = if !isempty(cutouts)
        Float64(first(cutouts).ri)
    elseif isapprox(computed_radius, Float64(outer.ro); rtol = 0, atol = tolerance)
        Float64(outer.ro)
    else
        computed_radius
    end
    isapprox(contact_radius, computed_radius; rtol = 0, atol = tolerance) || return plan
    sectors = [AnnularSector(Float64(sector.ri), Float64(sector.ro), Float64(sector.span),
                   convert(Placement{Float64}, sector.at)) for sector in cutouts]
    for sector in sectors
        isapprox(sector.at.x, center[1]; rtol = 0, atol = tolerance) || return plan
        isapprox(sector.at.y, center[2]; rtol = 0, atol = tolerance) || return plan
        isapprox(sector.ri, contact_radius; rtol = 0, atol = tolerance) || return plan
        isapprox(sector.ro, Float64(outer.ro); rtol = 0, atol = tolerance) || return plan
        1e-12 < sector.span < 2π - 1e-12 || return plan
    end
    shell = Float64(outer.ro) - contact_radius > tolerance
    shell || isempty(sectors) || return plan
    spans = NTuple{2, Float64}[]
    if shell
        if isempty(sectors)
            push!(spans, (0.0, 2π))
        else
            occupied = NTuple{2, Float64}[]
            for sector in sectors
                a = mod(sector.at.φ - sector.span / 2, 2π)
                isapprox(a, 2π; rtol = 0, atol = 128eps(Float64)) && (a = 0.0)
                b = a + sector.span
                if b <= 2π + 1e-12
                    push!(occupied, (a, min(b, 2π)))
                else
                    push!(occupied, (a, 2π))
                    push!(occupied, (0.0, b - 2π))
                end
            end
            sort!(occupied; by = first)
            merged = NTuple{2, Float64}[]
            for interval in occupied
                if isempty(merged) || interval[1] > last(merged)[2] + 1e-12
                    push!(merged, interval)
                    continue
                end
                interval[1] < last(merged)[2] - 1e-12 && return plan
                merged[end] = (last(merged)[1], max(last(merged)[2], interval[2]))
            end
            cursor = 0.0
            for interval in merged
                interval[1] - cursor > 1e-12 &&
                    push!(spans, (cursor, interval[1] - cursor))
                cursor = max(cursor, interval[2])
            end
            2π - cursor > 1e-12 && push!(spans, (cursor, 2π - cursor))
        end
    end
    sort!(entries; by = first)
    if length(entries) > 1
        angles = first.(entries)
        separations = [mod(angles[mod1(i + 1, length(angles))] - angles[i], 2π)
                       for i in eachindex(angles)]
        minimum(separations) > 1e-12 || return plan
    end
    return (center = center, inner_radius = Float64(outer.ri),
        outer_radius = contact_radius, shell_outer_radius = Float64(outer.ro),
        angles = first.(entries), holes = last.(entries), outer_spans = spans)
end

"""
$(TYPEDSIGNATURES)

Build the plane surfaces of `contour`, whose coordinates and radii are in \\[m\\]
and angles in \\[rad\\]. If loops share a point, walk their material boundaries.
Otherwise add plane surfaces with clockwise holes. Return their tags.
"""
function faces!(geometry::Geometry, contour::Contour)
    outer = draw(geometry, contour.outer)
    holes = [draw(geometry, path) for path in contour.holes]
    isempty(holes) && return Int[gmsh.model.geo.addPlaneSurface([outer.ccw])]
    seen = Set{Int}()
    touching = false
    for loop in (outer, holes...)
        points = Set(p for edge in loop.curves for p in geometry.ends[abs(edge)])
        touching |= any(p -> p in seen, points)
        union!(seen, points)
    end
    if !touching
        return Int[gmsh.model.geo.addPlaneSurface([outer.ccw; [hole.cw for hole in holes]])]
    end
    surfaces, outside = walk!(geometry, (outer.curves, (-hole.curves for hole in holes)...))
    isempty(surfaces) && throw(ArgumentError("FEM material has no positive-area face"))
    isempty(outside) || throw(ArgumentError("uncontained FEM material hole"))
    return surfaces
end

"""
$(TYPEDSIGNATURES)

Build the plane surfaces of `shape` through its contour. Dimensions and coordinates are
in \\[m\\] and angles in \\[rad\\]. Return tags of plane surfaces.
"""
function faces!(geometry::Geometry, shape::AbstractShape)
    return faces!(geometry, Contour(shape))
end

"""
$(TYPEDSIGNATURES)

Build the compartments of `shape` when its wires are tangent to its inner circle, followed
by the shell's free sectors. Otherwise build its contour's faces. Dimensions and coordinates
are in \\[m\\] and angles in \\[rad\\]. Return tags of plane surfaces.
"""
function faces!(geometry::Geometry, shape::DifferenceShape)
    plan = compartments(shape)
    isempty(plan.holes) && return faces!(geometry, Contour(shape))
    surfaces = Int[]
    cx, cy = plan.center
    for i in eachindex(plan.holes)
        j = mod1(i + 1, length(plan.holes))
        a, b = plan.angles[i], plan.angles[j]
        span = mod(b - a, 2π)
        span <= 1e-12 && (span = 2π)
        wire, next_wire = plan.holes[i], plan.holes[j]
        outer_arc = Arc(Placement(cx, cy, a), plan.outer_radius, span)
        next_at = Placement(next_wire.at.x, next_wire.at.y, b)
        next_wire_arc = Arc(next_at, next_wire.r, -π)
        inner_arc = Arc(Placement(cx, cy, b), plan.inner_radius, -span)
        wire_at = Placement(wire.at.x, wire.at.y, a + π)
        wire_arc = Arc(wire_at, wire.r, -π)
        tags = vcat(draw(geometry, outer_arc), draw(geometry, next_wire_arc),
            draw(geometry, inner_arc), draw(geometry, wire_arc))
        members, outside = walk!(geometry, (tags,))
        isempty(members) && throw(ArgumentError("FEM material has no positive-area face"))
        isempty(outside) || throw(ArgumentError("uncontained FEM material hole"))
        append!(surfaces, members)
    end
    for (a, span) in plan.outer_spans
        sector = AnnularSector(plan.outer_radius, plan.shell_outer_radius, span,
            Placement(cx, cy, a + span / 2))
        append!(surfaces, faces!(geometry, Contour(sector)))
    end
    return surfaces
end

"""
$(TYPEDSIGNATURES)

Return the direction \\[rad\\] of oriented `edge` at distance `ds` \\[m\\] from
its first point. Read the stored CAD curve.
"""
function direction(geometry::Geometry, edge::Int, ds::Real)
    return direction(geometry, geometry.curves[abs(edge)], edge, ds)
end

"""
$(TYPEDSIGNATURES)

Return the direction \\[rad\\] of the stored `segment` along oriented `edge`.
The probe distance `ds` is in \\[m\\].
"""
function direction(::Geometry, segment::Segment, edge::Int, ds::Real)
    a, b = segment(0), segment(1)
    return atan(sign(edge) * (b[2] - a[2]), sign(edge) * (b[1] - a[1]))
end

"""
$(TYPEDSIGNATURES)

Return the direction \\[rad\\] of the chord from the first point of oriented `edge`
to distance `ds` \\[m\\] along its stored `arc`. With sense ``s`` and radius ``r``
\\[m\\], the chord direction is the tangent direction plus ``s\\,ds/(2r)``.
"""
function direction(geometry::Geometry, arc::Arc, edge::Int, ds::Real)
    cx, cy, radius = key(arc.at.x, arc.at.y, arc.r)
    p, q = geometry.ends[abs(edge)]
    a, b = geometry.positions[p], geometry.positions[q]
    orientation = sign((a[1] - cx) * (b[2] - cy) - (a[2] - cy) * (b[1] - cx))
    xy = edge > 0 ? a : b
    sense = sign(edge) * orientation
    return atan(sense * (xy[1] - cx), -sense * (xy[2] - cy)) + sense * ds / (2radius)
end

"""
$(TYPEDSIGNATURES)

Return whether `point` lies inside `polygon` by the even-odd rule. All coordinates are in
\\[m\\]. The last vertex connects to the first.
"""
function isinside(polygon::Vector{NTuple{2, Float64}}, point::NTuple{2, Float64})
    inside = false
    for i in eachindex(polygon)
        a, b = polygon[i], polygon[mod1(i + 1, length(polygon))]
        (a[2] > point[2]) == (b[2] > point[2]) && continue
        point[1] < a[1] + (point[2] - a[2]) * (b[1] - a[1]) / (b[2] - a[2]) &&
            (inside = !inside)
    end
    return inside
end

"""
$(TYPEDSIGNATURES)

Return the signed area \\[m²\\] of the closed polyline `points` \\[m\\]:

```math
A = \\frac{1}{2} \\sum_i
    [(x_i-x_0)(y_{i+1}-y_0) - (x_{i+1}-x_0)(y_i-y_0)],
```

where the last vertex connects to the first and ``(x_0,y_0)`` is the first vertex.
Counterclockwise paths have positive area.
"""
function area(points::Vector{NTuple{2, Float64}})
    origin = first(points)
    return sum(eachindex(points)) do i
        j = mod1(i + 1, length(points))
        return (points[i][1] - origin[1]) * (points[j][2] - origin[2]) -
               (points[j][1] - origin[1]) * (points[i][2] - origin[2])
    end / 2
end

"""
$(TYPEDSIGNATURES)

Walk oriented `boundaries`, with material on the left, into cycles. Return the plane
surface tags of positive cycles with their contained holes, and the curve lists of
uncontained negative cycles. Stored coordinates are in \\[m\\]. Compare incident
directions at a common distance from each contact.
"""
function walk!(geometry::Geometry, boundaries)
    multiplicities = Dict{Int, Int}()
    for tags in boundaries, edge in tags

        multiplicities[abs(edge)] = get(multiplicities, abs(edge), 0) + sign(edge)
    end
    all(n -> abs(n) <= 1, values(multiplicities)) ||
        throw(ArgumentError("overlapping boundaries in a FEM material face"))
    edges = sort!([tag * n for (tag, n) in multiplicities if !iszero(n)]; by = abs)
    samples = Dict{Int, Vector{NTuple{2, Float64}}}()
    outgoing = Dict{Int, Vector{Int}}()
    endpoints = Dict{Int, NTuple{2, Int}}()
    for edge in edges
        p, q = geometry.ends[abs(edge)]
        curve = geometry.curves[abs(edge)]
        coordinates = NTuple{2, Float64}[curve(t)
                                         for t in
                                             (0.0, 1e-6, ((1:15) ./ 16)..., 1 - 1e-6, 1.0)]
        coordinates[1], coordinates[end] = geometry.positions[p], geometry.positions[q]
        if edge < 0
            p, q = q, p
            reverse!(coordinates)
        end
        endpoints[edge] = (p, q)
        samples[edge] = coordinates
        push!(get!(Vector{Int}, outgoing, p), edge)
    end
    distances = Dict{Int, Float64}()
    for (p, q) in values(endpoints)
        a, b = geometry.positions[p], geometry.positions[q]
        ds = hypot(b[1] - a[1], b[2] - a[2]) / 64
        for tag in (p, q)
            distances[tag] = min(get(distances, tag, Inf), ds)
        end
    end
    successors = Dict{Int, Int}()
    for edge in edges
        ds = distances[endpoints[edge][2]]
        backwards = direction(geometry, -edge, ds)
        candidates = get(outgoing, endpoints[edge][2], Int[])
        isempty(candidates) && throw(ArgumentError("open FEM material boundary"))
        successors[edge] = candidates[argmin(map(candidates) do candidate
            return mod(backwards - direction(geometry, candidate, ds), 2π)
        end)]
    end
    length(unique(values(successors))) == length(edges) ||
        throw(ArgumentError("ambiguous FEM material boundary at a contact"))
    cycles = Vector{Int}[]
    polygons = Vector{NTuple{2, Float64}}[]
    remaining = Set(edges)
    for initial in edges
        initial in remaining || continue
        cycle = Int[]
        polygon = NTuple{2, Float64}[]
        edge = initial
        while edge in remaining
            delete!(remaining, edge)
            push!(cycle, edge)
            append!(polygon, samples[edge][1:(end - 1)])
            edge = successors[edge]
        end
        edge == initial || throw(ArgumentError("nonclosing FEM material face"))
        push!(cycles, cycle)
        push!(polygons, polygon)
    end
    areas = area.(polygons)
    outers = findall(>(0), areas)
    interiors = Dict(i => Int[] for i in outers)
    outside = Int[]
    for i in findall(<(0), areas)
        coordinates = samples[first(cycles[i])]
        middle = length(coordinates) ÷ 2
        probe = ((coordinates[middle][1] + coordinates[middle + 1][1]) / 2,
            (coordinates[middle][2] + coordinates[middle + 1][2]) / 2)
        parents = filter(owner -> isinside(polygons[owner], probe), outers)
        if isempty(parents)
            push!(outside, i)
            continue
        end
        owner = parents[argmin(areas[parents])]
        push!(interiors[owner], i)
    end
    loops = [i in outside ? 0 : gmsh.model.geo.addCurveLoop(cycle)
             for (i, cycle) in enumerate(cycles)]
    surfaces = [gmsh.model.geo.addPlaneSurface([loops[i]; loops[interiors[i]]])
                for i in outers]
    return surfaces, cycles[outside]
end

"""
$(TYPEDSIGNATURES)

Partition `annulus` into four-sided patches at its aligned breaks. Coordinates and radii
are in \\[m\\] and angles in \\[rad\\]. Record each patch's region `i`, surface,
arcs, spokes, corners, and fraction of a revolution. Return tags of plane surfaces.
"""
function patches!(geometry::Geometry, annulus::Annulus, i::Int)
    contour = Contour(annulus)
    outer = reduce(vcat, [draw(geometry, curve) for curve in contour.outer]; init = Int[])
    inner = reduce(vcat, [draw(geometry, curve) for curve in only(contour.holes)]; init = Int[])
    surfaces = Int[]
    for (oc, ic) in zip(outer, inner)
        a, b = geometry.ends[abs(oc)]
        p, q = oc > 0 ? (a, b) : (b, a)
        a, b = geometry.ends[abs(ic)]
        ip, iq = ic > 0 ? (a, b) : (b, a)
        start_tag = only(draw(geometry, Segment(geometry.positions[ip], geometry.positions[p])))
        stop_tag = only(draw(geometry, Segment(geometry.positions[iq], geometry.positions[q])))
        loop = gmsh.model.geo.addCurveLoop([oc, -stop_tag, -ic, start_tag])
        surface = gmsh.model.geo.addPlaneSurface([loop])
        u = geometry.positions[p] .- (annulus.at.x, annulus.at.y)
        v = geometry.positions[q] .- (annulus.at.x, annulus.at.y)
        a = atan(abs(u[1] * v[2] - u[2] * v[1]), u[1] * v[1] + u[2] * v[2])
        push!(geometry.patches,
            (region = i, surface = surface,
                arcs = (abs(oc), abs(ic)), spokes = (abs(start_tag), abs(stop_tag)),
                corners = (p, q, iq, ip), fraction = a / (2π)))
        push!(surfaces, surface)
    end
    return surfaces
end

"""
$(TYPEDSIGNATURES)

Build the geometry of `model` in an initialized Gmsh session under `options`. Add the model
named by the system's `id` and leave it in the session for the caller to remove. The native
`ThinLayerStrips` control selects patches for passive annuli. Coordinates and radii are
in \\[m\\] and angles in \\[rad\\]. Return the completed geometry with empty anchors.
"""
function build(::Type{Geometry}, model::FieldModel, options::ComputationOptions{GmshBackend})
    system = model.problem.system
    gmsh.model.add(system.id)
    geometry = Geometry(model)
    for region in system.regions
        cut!(geometry, region.shape)
    end
    contours = [Contour(resolve(at, design.outer))
                for (at, design) in zip(system.placements, system.designs)]
    for contour in contours
        cut!(geometry, contour)
    end
    cut!(geometry)
    annuli = findall(
        i -> system.regions[i].shape isa Annulus &&
             !model.materials[model.material[i]].conductor,
        eachindex(system.regions))
    strips = options.overrides.ThinLayerStrips == 1
    if strips
        cut!(geometry, Annulus[system.regions[i].shape for i in annuli])
    end
    for (i, region) in enumerate(system.regions)
        geometry.surfaces[i] = strips && i in annuli ?
                               patches!(geometry, region.shape, i) :
                               faces!(geometry, region.shape)
    end
    for contour in contours
        push!(geometry.outlines, draw(geometry, contour.outer))
    end
    for (in_air, holes, pockets) in ((true, geometry.holes.air, geometry.pockets.air),
        (false, geometry.holes.earth, geometry.pockets.earth))
        append!(holes,
            [-reverse(outline.curves)
             for (outline, air) in zip(geometry.outlines, model.air)
             if air == in_air])
        points = [Set(p for edge in tags for p in geometry.ends[abs(edge)])
                  for tags in holes]
        any(i -> any(j -> !isdisjoint(points[i], points[j]), 1:(i - 1)),
            eachindex(points)) || continue
        members, outside = walk!(geometry, holes)
        append!(pockets, members)
        empty!(holes)
        append!(holes, outside)
    end
    gmsh.model.geo.synchronize()
    return geometry
end
