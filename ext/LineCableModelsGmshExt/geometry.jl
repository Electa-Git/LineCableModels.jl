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
    "Ellipse centers and rotations \\[m, rad\\], semi-axes and normal offsets \\[m\\], sorted endpoint tags, and Gmsh tags of elliptic arcs and their offset splines."
    elliptic_arcs::Dict{Tuple{NTuple{6, Float64}, NTuple{2, Int}}, Int}
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
    "The strips of the conductors and the depth of each \\[m\\]."
    strips::Vector{@NamedTuple{region::Int, surface::Int, curves::NTuple{2, Vector{Int}},
        spokes::NTuple{2, Int}, corners::NTuple{4, Int}, depth::Float64}}
    "Counterclockwise and clockwise outline loops and oriented curves of each cable."
    outlines::Vector{@NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}}
    "Clockwise curve lists of the holes in each medium's main face."
    holes::@NamedTuple{air::Vector{Vector{Int}}, earth::Vector{Vector{Int}}}
    "Plane surface tags of the enclosed pockets of each medium."
    pockets::@NamedTuple{air::Vector{Int}, earth::Vector{Int}}
    "Terminal anchor coordinates and metal dimensions \\[m\\]."
    anchors::Vector{@NamedTuple{x::Float64, y::Float64, dimension::Float64}}
    "Segments per turn of the circle rule, for curves drawn through points \\[dimensionless\\]."
    segments::Int
end

"""
$(TYPEDSIGNATURES)

Construct empty geometric registries for `model`. Allocate an empty vector of plane surfaces
for each region. Set ``segments=segments(options.overrides.ConductorGeometryTolerance)``
for curves drawn through points.
"""
function Geometry(model::FieldModel, options::ComputationOptions{GmshBackend})
    return Geometry(model,
        Dict{NTuple{2, Float64}, Int}(),
        Dict{Int, NTuple{2, Float64}}(),
        Dict{NTuple{2, Int}, Vector{NTuple{2, Float64}}}(),
        Dict{NTuple{2, Int}, Int}(),
        Dict{Tuple{NTuple{3, Float64}, NTuple{2, Int}}, Int}(),
        Dict{Tuple{NTuple{6, Float64}, NTuple{2, Int}}, Int}(),
        Dict{Int, NTuple{2, Int}}(),
        Dict{Int, Curve{Float64}}(),
        Dict{NTuple{3, Float64}, Vector{Float64}}(),
        Dict{Tuple{NTuple{3, Float64}, Float64}, NTuple{2, Float64}}(),
        Dict{Vector{Int}, @NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}}(),
        [Int[] for _ in model.problem.system.regions],
        @NamedTuple{region::Int, surface::Int, arcs::NTuple{2, Int},
            spokes::NTuple{2, Int}, corners::NTuple{4, Int}, fraction::Float64}[],
        @NamedTuple{region::Int, surface::Int, curves::NTuple{2, Vector{Int}},
            spokes::NTuple{2, Int}, corners::NTuple{4, Int}, depth::Float64}[],
        @NamedTuple{ccw::Int, cw::Int, curves::Vector{Int}}[],
        (air = Vector{Int}[], earth = Vector{Int}[]),
        (air = Int[], earth = Int[]),
        @NamedTuple{x::Float64, y::Float64, dimension::Float64}[],
        segments(options.overrides.ConductorGeometryTolerance))
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
function cut!(geometry::Geometry, shape::PerforatedShape)
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
A pair is tangent when one of its circles is tangent to the other by `istangent`. Centers,
radii, and contact coordinates are in \\[m\\]. Return `geometry`.
"""
function cut!(geometry::Geometry)
    circles = collect(keys(geometry.breaks))
    for i in eachindex(circles)
        cx, cy, r = circles[i]
        for j in (i + 1):length(circles)
            qx, qy, s = circles[j]
            dx, dy = qx - cx, qy - cy
            distance = hypot(dx, dy)
            distance <= 1e-12 * max(r, s, distance, 1.0) && continue
            p = Circle(r, Placement(cx, cy))
            q = Circle(s, Placement(qx, qy))
            tangent = istangent(p, q) || istangent(q, p)
            tangent || abs(r - s) < distance < r + s || continue
            along = (r^2 - s^2 + distance^2) / (2distance)
            height = tangent ? 0.0 : sqrt(r^2 - along^2)
            ux, uy = dx / distance, dy / distance
            bx, by = cx + along * ux, cy + along * uy
            contacts = tangent ? ((bx, by),) :
                       ((bx - height * uy, by + height * ux),
                (bx + height * uy, by - height * ux))
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
function compartments(shape::PerforatedShape)
    return compartments(shape, shape.outer)
end

"""
$(TYPEDSIGNATURES)

Return a plan with no wires for a `shape` whose `outer` is not an annulus. Plan radii and
coordinates are in \\[m\\], and angles are in \\[rad\\].
"""
function compartments(::PerforatedShape, ::AbstractShape)
    return (center = (0.0, 0.0), inner_radius = 0.0, outer_radius = 0.0,
        shell_outer_radius = 0.0, angles = Float64[], holes = Circle{Float64}[],
        outer_spans = NTuple{2, Float64}[])
end

"""
$(TYPEDSIGNATURES)

Return the course of circular wires tangent to the inner circle of `outer`, or a plan
with no wires. Each wire touches the inner circle and the contact circle, and each cutout
touches the contact circle and the outer circle, by `istangent`. The contact circle is
the inner circle of the cutouts. Without cutouts, it is the outer circle when every wire
touches it, else the circle through the outer point of the first wire. Record the center
and inner, contact, and shell radii \\[m\\], the wires and their angles \\[rad\\], and the
shell's free angular spans \\[rad\\]. Convert the cutouts of `shape` to `Float64` and
calculate their free spans.
"""
function compartments(shape::PerforatedShape, outer::Annulus)
    plan = invoke(compartments, Tuple{PerforatedShape, AbstractShape}, shape, outer)
    isempty(shape.holes) && return plan
    iszero(outer.ri) && return plan
    center = (Float64(outer.at.x), Float64(outer.at.y))
    wires = [hole for hole in shape.holes if hole isa Circle]
    cutouts = [hole
               for hole in shape.holes
               if hole isa Sector && iszero(hole.apex) && iszero(nominal(hole.fillet))]
    length(wires) + length(cutouts) == length(shape.holes) || return plan
    isempty(wires) && return plan
    all(wire -> istangent(wire, Circle(outer.ri, outer.at)), wires) || return plan
    rim = Circle(outer.ro, outer.at)
    contact_radius = if !isempty(cutouts)
        Float64(first(cutouts).ri)
    elseif all(wire -> istangent(wire, rim), wires)
        Float64(outer.ro)
    else
        hypot(Float64(first(wires).at.x) - center[1],
            Float64(first(wires).at.y) - center[2]) + Float64(first(wires).r)
    end
    contact = Circle(contact_radius, outer.at)
    all(wire -> istangent(wire, contact), wires) || return plan
    all(sector -> istangent(sector, contact) && istangent(sector, rim), cutouts) ||
        return plan
    entries = Tuple{Float64, Circle{Float64}}[]
    for wire in wires
        dx, dy = Float64(wire.at.x) - center[1], Float64(wire.at.y) - center[2]
        a = mod(atan(dy, dx), 2π)
        isapprox(a, 2π; rtol = 0, atol = 128eps(Float64)) && (a = 0.0)
        push!(entries, (a, Circle(Float64(wire.r), convert(Placement{Float64}, wire.at))))
    end
    sectors = [Sector(
                   Float64(sector.ri), Float64(sector.ro), Float64(sector.span), 0.0, 0.0,
                   convert(Placement{Float64}, sector.at)) for sector in cutouts]
    tolerance = 1e-12 * max(Float64(outer.ro), contact_radius, 1.0)
    for sector in sectors
        isapprox(sector.at.x, center[1]; rtol = 0, atol = tolerance) || return plan
        isapprox(sector.at.y, center[2]; rtol = 0, atol = tolerance) || return plan
        1e-12 < sector.span < 2π - 1e-12 || return plan
    end
    shell = contact_radius != Float64(outer.ro)
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
function faces!(geometry::Geometry, shape::PerforatedShape)
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
        sector = Sector(plan.outer_radius, plan.shell_outer_radius, span, 0.0, 0.0,
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

Return the direction \\[rad\\] of the chord from the first point of oriented `edge` to
distance `ds` \\[m\\] along its stored elliptic `arc`. With sense ``s`` and curvature ``κ``
\\[1/m\\], the chord direction is the tangent direction, in ``[-π, π]``, plus ``s\\,ds\\,κ/2``.
The curvature is that of the offset arc at the first point,

```math
κ = \\frac{κ_0}{1 + o κ_0}, \\qquad κ_0 = \\frac{ab}{(a^2 \\sin^2 θ + b^2 \\cos^2 θ)^{3/2}},
```

where ``a`` and ``b`` \\[m\\] are the semi-axes, ``o`` \\[m\\] the offset and ``θ`` \\[rad\\]
the parametric angle of the point.
"""
function direction(geometry::Geometry, arc::EllipticArc, edge::Int, ds::Real)
    p, q = geometry.ends[abs(edge)]
    xy = geometry.positions[edge > 0 ? p : q]
    start, finish = arc(0), arc(1)
    forward = hypot(xy[1] - start[1], xy[2] - start[2]) <=
              hypot(xy[1] - finish[1], xy[2] - finish[2])
    t = forward ? 0 : 1
    θ = arc.φ0 + t * arc.span
    κ₀ = arc.a * arc.b / (arc.a^2 * sin(θ)^2 + arc.b^2 * cos(θ)^2)^(3 / 2)
    κ = κ₀ / (1 + arc.offset * κ₀)
    sense = forward ? sign(arc.span) : -sign(arc.span)
    return rem2pi(direction(arc, t) + (forward ? 0 : π), RoundNearest) + sense * ds * κ / 2
end

"""
$(TYPEDSIGNATURES)

Walk oriented `boundaries`, with material on the left, into cycles. Return the tag of
the Gmsh plane surface of each positive cycle with its contained holes, and the curve
lists of uncontained negative cycles. Stored coordinates are in \\[m\\]. Compare incident
directions at a common distance from each contact.

A cycle's path takes the stored curve of each edge, reversed where the edge is negative.
The sign of the path's `area` is the cycle's orientation. A negative cycle lies in the
smallest positive cycle whose region contains the middle of its first curve.
"""
function walk!(geometry::Geometry, boundaries)
    multiplicities = Dict{Int, Int}()
    for tags in boundaries, edge in tags

        multiplicities[abs(edge)] = get(multiplicities, abs(edge), 0) + sign(edge)
    end
    all(n -> abs(n) <= 1, values(multiplicities)) ||
        throw(ArgumentError("overlapping boundaries in a FEM material face"))
    edges = sort!([tag * n for (tag, n) in multiplicities if !iszero(n)]; by = abs)
    outgoing = Dict{Int, Vector{Int}}()
    endpoints = Dict{Int, NTuple{2, Int}}()
    for edge in edges
        p, q = geometry.ends[abs(edge)]
        if edge < 0
            p, q = q, p
        end
        endpoints[edge] = (p, q)
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
    paths = Vector{Curve{Float64}}[]
    remaining = Set(edges)
    for initial in edges
        initial in remaining || continue
        cycle = Int[]
        path = Curve{Float64}[]
        edge = initial
        while edge in remaining
            delete!(remaining, edge)
            push!(cycle, edge)
            push!(path, edge > 0 ? geometry.curves[edge] : reverse(geometry.curves[-edge]))
            edge = successors[edge]
        end
        edge == initial || throw(ArgumentError("nonclosing FEM material face"))
        push!(cycles, cycle)
        push!(paths, path)
    end
    areas = area.(paths)
    outers = findall(>(0), areas)
    contours = Dict(i => Contour(paths[i], Vector{Curve{Float64}}[]) for i in outers)
    interiors = Dict(i => Int[] for i in outers)
    outside = Int[]
    for i in findall(<(0), areas)
        probe = first(paths[i])(1 / 2)
        parents = filter(owner -> probe ∈ contours[owner], outers)
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
        push!(geometry.patches,
            (region = i, surface = surface,
                arcs = (abs(oc), abs(ic)), spokes = (abs(start_tag), abs(stop_tag)),
                corners = (p, q, iq, ip),
                fraction = turn(geometry.curves[abs(oc)]).fraction))
        push!(surfaces, surface)
    end
    return surfaces
end

"""
$(TYPEDSIGNATURES)

Partition the conducting region `i` of `shape` into four-sided strips and a core. The strips
lie between the contour and its `parallel` at the depth ``e`` \\[m\\] inside the metal, the
depth of its region's [`section`](@ref). Each curve of each path of the contour bounds one
strip. Its outer side is the curve's drawn pieces, and its inner side the drawn pieces of
its parallel. A spoke joins each junction of the contour to the junction of the parallels.

Each strip's loop runs counterclockwise. The parallels bound the core. Record each strip's
region, surface, side curves, spokes, corners and depth ``e``. Return the tags of the strips
and of the core.

# Errors

- Throws `DomainError` at a corner whose interior angle exceeds π, where a path of the
  parallel gains a curve.
"""
function strips!(geometry::Geometry, shape::AbstractShape, i::Int, depth::Real)
    contour = Contour(shape)
    inner = parallel(contour, -depth)
    for (path, innerpath) in
        zip((contour.outer, contour.holes...), (inner.outer, inner.holes...))
        length(innerpath) > length(path) &&
            throw(DomainError(shape,
                "region $i has a reflex corner: the inward parallel of its outline " *
                "gains a curve there, so the outline has no strips"))
    end
    surfaces = Int[]
    loops = Int[]
    for (k, (path, innerpath)) in
        enumerate(zip((contour.outer, contour.holes...), (inner.outer, inner.holes...)))
        tags = Int[]
        for j in eachindex(path)
            outer = draw(geometry, path[j])
            pieces = draw(geometry, innerpath[j])
            a, b = geometry.ends[abs(first(outer))]
            p = first(outer) > 0 ? a : b
            a, b = geometry.ends[abs(last(outer))]
            q = last(outer) > 0 ? b : a
            a, b = geometry.ends[abs(first(pieces))]
            p′ = first(pieces) > 0 ? a : b
            a, b = geometry.ends[abs(last(pieces))]
            q′ = last(pieces) > 0 ? b : a
            s = only(draw(geometry, Segment(geometry.positions[p], geometry.positions[p′])))
            t = only(draw(geometry, Segment(geometry.positions[q], geometry.positions[q′])))
            loop = k == 1 ? [outer; t; -reverse(pieces); -s] :
                   [s; pieces; -t; -reverse(outer)]
            corners = k == 1 ? (p, q, q′, p′) : (p, p′, q′, q)
            surface = gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop(loop)])
            push!(geometry.strips,
                (region = i, surface = surface, curves = (abs.(outer), abs.(pieces)),
                    spokes = (s, t), corners = corners, depth = depth))
            push!(surfaces, surface)
            append!(tags, pieces)
        end
        push!(loops, gmsh.model.geo.addCurveLoop(k == 1 ? tags : -reverse(tags)))
    end
    push!(surfaces, gmsh.model.geo.addPlaneSurface(loops))
    return surfaces
end

"""
$(TYPEDSIGNATURES)

Build the geometry of `model` in an initialized Gmsh session under `options`. Add the model
named by the system's `id` and leave it in the session for the caller to remove. The native
`ThinLayerStrips` control selects patches for passive annuli. Coordinates and radii are
in \\[m\\] and angles in \\[rad\\]. Carry ``segments=segments(τ)`` from the options
to the drawing. Return the completed geometry with terminal anchors.
"""
function build(::Type{Geometry}, model::FieldModel, options::ComputationOptions{GmshBackend})
    system = model.problem.system
    gmsh.model.add(system.id)
    geometry = Geometry(model, options)
    for region in system.regions
        cut!(geometry, region.shape)
    end
    contours = [Contour(resolve(at, design.outer))
                for (at, design) in zip(system.placements, system.designs)]
    for contour in contours
        cut!(geometry, contour)
    end
    cut!(geometry)
    annuli = findall(model.layer)
    strips = options.overrides.ThinLayerStrips == 1
    if strips
        cut!(geometry, Annulus[system.regions[i].shape for i in annuli])
    end
    conductor = [model.materials[g].conductor for g in model.material]
    for (i, region) in enumerate(system.regions)
        depth = section(model, i, options.overrides).depth
        geometry.surfaces[i] = strips && model.layer[i] ?
                               patches!(geometry, region.shape, i) :
                               conductor[i] && depth > 0 ?
                               strips!(geometry, region.shape, i, depth) :
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
    append!(geometry.anchors, [anchor(geometry, t) for t in eachindex(system.terminals)])
    return geometry
end

"""
$(TYPEDSIGNATURES)

Find the measurement anchor of terminal `t` in `geometry`. Probe upward from each region's
lowest point. Its metal dimension sets the vertical and horizontal offsets. Accept the
first point that lies inside the metal for horizontal shifts of 0 and 1 % of the smallest
metal dimension of the terminal.

# Returns

- `x`, `y`: anchor coordinates \\[m\\].
- `dimension`: the smallest metal dimension of the terminal \\[m\\].

# Errors

- Throws `ArgumentError` when no probe lies inside the terminal's metal.
"""
function anchor(geometry::Geometry, t::Integer)
    system = geometry.model.problem.system
    indices = system.terminals[t].regions
    dimensions = [dimension(system.regions[i].shape)
                  for i in indices]
    minimum_dimension = minimum(dimensions)
    for (i, width) in zip(indices, dimensions)
        low = lowest(Contour(system.regions[i].shape))
        for dy in (0.125, 0.25, 0.5, 0.75), dx in (0.0, 0.125, -0.125, 0.25, -0.25)

            xy = (low[1] + dx * width, low[2] + dy * width)
            if any(geometry.surfaces[i]) do surface
                return all(
                    shift -> gmsh.model.isInside(2, surface,
                        [xy[1] + shift * minimum_dimension, xy[2], 0.0]) > 0,
                    (0.0, 0.01))
            end
                return (x = Float64(xy[1]), y = Float64(xy[2]),
                    dimension = Float64(minimum_dimension))
            end
        end
    end
    return throw(ArgumentError("cannot locate an interior metal anchor for terminal $t"))
end

"""
$(TYPEDSIGNATURES)

Return the sorted unique boundary curve tags of the plane `surfaces`. Gmsh removes curves
shared by two surfaces.
"""
function curves(surfaces::AbstractVector{<:Integer})
    isempty(surfaces) && return Int[]
    boundaries = gmsh.model.getBoundary([(2, surface) for surface in surfaces],
        true, false, false)
    tags = [Int(tag) for (dim, tag) in boundaries if dim == 1]
    return sort!(unique(tags))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of region `i` of `model` under the native `controls`. The
parameters are a named tuple of five numbers:

- `target`, the target of the element size in the region \\[m\\], before `MeshSizeFactor`.
- `width`, the width of the law of its conductor \\[m\\].
- `divisions`, the number of parts of that width in the bulk cap.
- `fraction`, the part of that width that the graded layer can take.
- `depth`, the depth of its strips \\[m\\], zero without strips.

The construction of the region selects the target of a strand or of a region laid alone.
`ConductorThicknessElements` gives the divisions of an annular region. A patched passive
annulus, a region of `model.layer` under `ThinLayerStrips` 1, takes at least the target
``2π r_o/n`` \\[m\\]. That is the arc of one of its ``n = segments(τ)`` segments on its
outer circle ``r_o``, for ``τ`` the `ConductorGeometryTolerance` \\[m\\].
"""
function section(model::FieldModel, i::Integer, controls::NamedTuple)
    region = model.problem.system.regions[i]
    parameters = section(
        region.source.primitive, region.shape, model.construction[i], controls)
    model.layer[i] && controls.ThinLayerStrips == 1 || return parameters
    bound = 2π * region.shape.ro / segments(controls.ConductorGeometryTolerance)
    return merge(parameters, (target = max(parameters.target, bound),))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the perforated primitive `x`: the parameters of its outer,
with the resolved `shape`.
"""
function section(x::PerforatedShape, shape::AbstractShape,
        construction::Union{Strand, Solid}, controls::NamedTuple)
    return section(x.outer, shape, construction, controls)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of an offset by its resolved `shape`.
"""
function section(::Offset, shape::AbstractShape, construction::Union{Strand, Solid},
        controls::NamedTuple)
    return section(shape, shape, construction, controls)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of `offset` with the perforated resolved `shape`, through the
outer of `shape`.
"""
function section(
        offset::Offset, shape::PerforatedShape, construction::Union{Strand, Solid},
        controls::NamedTuple)
    return section(offset, shape, shape.outer, construction, controls)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of an offset fill whose container is the annulus `outer`: the
parameters of the annulus.
"""
function section(::Offset, ::PerforatedShape, outer::Annulus,
        construction::Union{Strand, Solid}, controls::NamedTuple)
    return section(outer, outer, construction, controls)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of an offset fill whose container is the layer `outer`: the
parameters of the layer.
"""
function section(offset::Offset, ::PerforatedShape, outer::PerforatedShape,
        construction::Union{Strand, Solid}, controls::NamedTuple)
    return section(offset, outer, construction, controls)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the round strand `x`, a circle of radius ``r`` \\[m\\]. The
target is ``r``, with the width ``r``, five divisions, the fraction ``0.5`` and no strips.
"""
function section(x::Circle, ::AbstractShape, ::Strand, ::NamedTuple)
    return (target = x.r, width = x.r, divisions = 5, fraction = 0.5, depth = zero(x.r))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the round region `x` laid alone, a circle. They are those of
a round strand, with the target ``r/5``.
"""
function section(x::Circle, ::AbstractShape, ::Solid, ::NamedTuple)
    return (target = x.r / 5, width = x.r, divisions = 5, fraction = 0.5, depth = zero(x.r))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the round strand `x`, an ellipse of radius ``r`` \\[m\\]. The
target is ``r``, with the width ``r``, five divisions, the fraction ``0.5`` and no strips.
The radius ``r`` is half the `dimension` of the ellipse: the depth at which its inward
parallel first degenerates, as for a circle.
"""
function section(x::Ellipse, ::AbstractShape, ::Strand, ::NamedTuple)
    r = dimension(x) / 2
    return (target = r, width = r, divisions = 5, fraction = 0.5, depth = zero(r))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the round region `x` laid alone, an ellipse. They are those
of a round strand, with the target ``r/5``.
"""
function section(x::Ellipse, ::AbstractShape, ::Solid, ::NamedTuple)
    r = dimension(x) / 2
    return (target = r / 5, width = r, divisions = 5, fraction = 0.5, depth = zero(r))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the annular region `x`, with the radii ``r_i`` and ``r_o``
\\[m\\] and the thickness ``w = r_o - r_i``. The target is ``w/2``, with the width ``w``,
`ConductorThicknessElements` divisions, the fraction ``0.45`` and no strips.
"""
function section(x::Annulus, ::AbstractShape, ::Union{Strand, Solid}, controls::NamedTuple)
    w = x.ro - x.ri
    return (target = w / 2, width = w,
        divisions = controls.ConductorThicknessElements, fraction = 0.45, depth = zero(w))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the layer `shape`, whose thickness ``w = dimension(shape)``
\\[m\\] is the same at every placement. The target is ``w/2``, with the width ``w``,
`ConductorThicknessElements` divisions, the fraction ``0.45`` and no strips.
"""
function section(::Offset, shape::PerforatedShape, ::AbstractShape, ::Union{Strand, Solid},
        controls::NamedTuple)
    w = dimension(shape)
    return (target = w / 2, width = w,
        divisions = controls.ConductorThicknessElements, fraction = 0.45, depth = zero(w))
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the strip strand `x`, a sector with the sides ``a`` and ``b``
\\[m\\] of the rectangle that it bends. The target is ``\\min(a, b)/2``. Its resolved
`shape` gives the width ``d = dimension(shape)`` \\[m\\], with five divisions, the fraction
``0.45`` and strips at the depth ``0.45 d``. The sides are

```math
a = r_o - r_i, \\qquad b = \\varphi (r_o + r_i) / 2,
```

where ``r_o`` and ``r_i`` are the outer and inner radii \\[m\\] and ``\\varphi`` is the
angular span \\[rad\\].
"""
function section(x::Sector, shape::AbstractShape, ::Strand, ::NamedTuple)
    a, b = x.ro - x.ri, x.span * (x.ro + x.ri) / 2
    width = dimension(shape)
    fraction = 0.45
    return (target = min(a, b) / 2, width = width,
        divisions = 5, fraction = fraction, depth = fraction * width)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the strip region `x` laid alone, a sector. They are those of
a strip strand, with the target ``\\min(a, b)/5``.
"""
function section(x::Sector, shape::AbstractShape, ::Solid, ::NamedTuple)
    a, b = x.ro - x.ri, x.span * (x.ro + x.ri) / 2
    width = dimension(shape)
    fraction = 0.45
    return (target = min(a, b) / 5, width = width,
        divisions = 5, fraction = fraction, depth = fraction * width)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the strip strand `x`, a rectangle, a polygon or an arc polygon
of `dimension` ``d_x`` \\[m\\]. A rectangle's `dimension` is its smaller side. The target is
``d_x/2``. Its resolved `shape` gives the width ``d = dimension(shape)`` \\[m\\], with five
divisions, the fraction ``0.45`` and strips at the depth ``0.45 d``.
"""
function section(x::Union{Rectangle, Polygon, ArcPolygon}, shape::AbstractShape, ::Strand,
        ::NamedTuple)
    d = dimension(x)
    width = dimension(shape)
    fraction = 0.45
    return (target = d / 2, width = width,
        divisions = 5, fraction = fraction, depth = fraction * width)
end

"""
$(TYPEDSIGNATURES)

Return the mesh parameters of the strip region `x` laid alone, a rectangle, a polygon or an
arc polygon. They are those of a strip strand, with the target ``d_x/5``.
"""
function section(x::Union{Rectangle, Polygon, ArcPolygon}, shape::AbstractShape, ::Solid,
        ::NamedTuple)
    d = dimension(x)
    width = dimension(shape)
    fraction = 0.45
    return (target = d / 5, width = width,
        divisions = 5, fraction = fraction, depth = fraction * width)
end

"""
$(TYPEDSIGNATURES)

Write `geometry` to `path` as the frequency-independent physical CAD read by the native
sources. Write Gmsh's unrolled geometry to `path * "_unrolled"`, replace each point's size
with `MeshFine`, append the native lists and physical groups, and remove the unrolled file.
Coordinates, radii, thicknesses, arc lengths, and anchor dimensions are in \\[m\\].
Angular fractions are dimensionless, and the angular match is in \\[rad\\].
Numbers use 17 significant digits. Each patch lists the pairs of circular arcs that touch
inside one another at its arcs' ends, the larger first, with both fractions.
"""
function Base.write(path::AbstractString, geometry::Geometry)
    gmsh.model.geo.synchronize()
    model = geometry.model
    system = model.problem.system
    material_surfaces = [Int[] for _ in model.materials]
    for (i, g) in enumerate(model.material)
        append!(material_surfaces[g], geometry.surfaces[i])
    end
    terminal_surfaces = [reduce(vcat, geometry.surfaces[terminal.regions]; init = Int[])
                         for terminal in system.terminals]
    cable_curves = [abs.(outline.curves) for outline in geometry.outlines]
    layers = [(shape = region.shape, curves = curves(geometry.surfaces[i]))
              for (i, region) in enumerate(system.regions)
              if model.layer[i]]
    conductor_boundaries = [curves(geometry.surfaces[i])
                            for i in eachindex(system.regions)
                            if model.materials[model.material[i]].conductor]
    conductor_curves = reduce(vcat, conductor_boundaries; init = Int[])
    layer_boundaries = reduce(vcat, [layer.curves for layer in layers]; init = Int[])
    layer_curves = sort!(setdiff(unique(layer_boundaries), conductor_curves))
    boundaries = map(layer_curves) do edge
        adjacent = filter(layer -> edge in layer.curves, layers)
        shape = first(adjacent).shape
        p = first(geometry.ends[edge])
        xy = geometry.positions[p]
        radius = hypot(xy[1] - shape.at.x, xy[2] - shape.at.y)
        thickness = minimum(dimension(layer.shape) for layer in adjacent)
        points = unique(last.(gmsh.model.getBoundary([(1, edge)], false, false, false)))
        cables = [k - 1 for (k, tags) in enumerate(cable_curves) if edge in tags]
        return (radius = radius, thickness = thickness, points = points,
            cables = cables, fraction = turn(geometry.curves[abs(edge)]).fraction)
    end
    incident = Dict{Int, Vector{Int}}()
    for ((circle, (p, q)), tag) in geometry.arcs
        push!(get!(incident, p, Int[]), tag)
        push!(get!(incident, q, Int[]), tag)
    end
    tangencies = map(geometry.patches) do patch
        entries = NTuple{4, Float64}[]
        for A in patch.arcs, e in geometry.ends[A], N in get(incident, e, Int[])
            N == A && continue
            a = first(geometry.ends[A]) == e ? A : -A
            n = first(geometry.ends[N]) == e ? N : -N
            ds = turn(geometry.curves[A]).length
            δa = direction(geometry, a, ds) - direction(geometry, a, 0)
            δn = direction(geometry, n, ds) - direction(geometry, n, 0)
            Δ = direction(geometry, n, 0) - direction(geometry, a, 0)
            cos(Δ) > 0 && abs(sin(Δ)) <= tolerance(one(sin(Δ))) || continue
            sign(δa) * sign(δn) > 0 || continue
            L, S = abs(δa) < abs(δn) ? (A, N) : (N, A)
            push!(entries, (L, S, turn(geometry.curves[L]).fraction,
                turn(geometry.curves[S]).fraction))
        end
        return entries
    end
    raw = path * "_unrolled"
    gmsh.write(raw)
    open(path, "w") do io
        println(io, "// Frequency-independent physical CAD; numerical sizing is native.")
        for line in eachline(raw)
            startswith(line, "Transfinite") && continue
            startswith(line, "Recombine") && continue
            occursin(r"^cl__\d+ =", line) && continue
            matched = match(r"^Point\((\d+)\)", line)
            if matched !== nothing
                tag = parse(Int, matched[1])
                xyz = gmsh.model.getValue(0, tag, Float64[])
                numbers = [@sprintf("%.17g", Float64(number)) for number in xyz]
                println(io, "Point(", tag, ") = {", join(numbers, ", "), ", MeshFine};")
            else
                println(io, line)
            end
        end
        boundary_radii = [boundary.radius for boundary in boundaries]
        boundary_thicknesses = [boundary.thickness for boundary in boundaries]
        boundary_fractions = [boundary.fraction for boundary in boundaries]
        for (name, entries) in (("FEMPolygonCurves", Int[]),
            ("FEMPolygonFirstPoints", Int[]), ("FEMPolygonLastPoints", Int[]),
            ("FEMPolygonLengths", Int[]), ("FEMPolygonAreaScales", Int[]),
            ("FEMMaterialSurfaces", reduce(vcat, material_surfaces; init = Int[])),
            ("FEMCableCurves", unique(reduce(vcat, cable_curves; init = Int[]))),
            ("FEMLayerBoundaryCurves", layer_curves),
            ("FEMLayerBoundaryRadii", boundary_radii),
            ("FEMLayerBoundaryThicknesses", boundary_thicknesses),
            ("FEMLayerBoundaryFractions", boundary_fractions))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, name, "() = {", join(numbers, ", "), "};")
        end
        for (i, patch) in enumerate(geometry.patches)
            for (name, entries) in (("Arcs", patch.arcs), ("Spokes", patch.spokes),
                ("Corners", patch.corners))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "FEMAnnular", name, "~{", i - 1, "}() = {",
                    join(numbers, ", "), "};")
            end
            println(io, "FEMAnnularSurface~{", i - 1, "} = ", patch.surface, ";")
            println(io, "FEMAnnularRegion~{", i - 1, "} = ", patch.region - 1, ";")
            shape = system.regions[patch.region].shape
            for (name, number) in (("Fraction", patch.fraction), ("Radius", shape.ro),
                ("Thickness", dimension(shape)))
                println(io, "FEMAnnular", name, "~{", i - 1, "} = ",
                    @sprintf("%.17g", Float64(number)), ";")
            end
            for (name, entries) in (
                ("Contacts", [number for pair in tangencies[i] for number in pair[1:2]]),
                ("ContactFractions",
                    [number for pair in tangencies[i] for number in pair[3:4]]))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "FEMAnnular", name, "~{", i - 1, "}() = {",
                    join(numbers, ", "), "};")
            end
        end
        println(io, "FEMAnnularPatches = ", length(geometry.patches), ";")
        println(io, "FEMAngularMatch = ", @sprintf("%.17g", tolerance(1.0)), ";")
        for (k, strip) in enumerate(geometry.strips)
            println(io, "FEMStripRegion~{", k - 1, "} = ", strip.region - 1, ";")
            println(io, "FEMStripSurface~{", k - 1, "} = ", strip.surface, ";")
            for (name, entries) in
                (("OuterCurves", strip.curves[1]), ("InnerCurves", strip.curves[2]),
                ("InnerLengths",
                    [turn(geometry.curves[c]).length for c in strip.curves[2]]),
                ("Spokes", strip.spokes), ("Corners", strip.corners))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(
                    io, "FEMStrip", name, "~{", k - 1, "}() = {", join(numbers, ", "), "};")
            end
        end
        println(io, "FEMStrips = ", length(geometry.strips), ";")
        for (i, boundary) in enumerate(boundaries)
            lists = (("Points", boundary.points), ("Cables", boundary.cables))
            for (name, entries) in lists
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "FEMLayerBoundary", name, "~{", i - 1, "}() = {",
                    join(numbers, ", "), "};")
            end
        end
        for (medium, holes, pockets) in (("Air", geometry.holes.air, geometry.pockets.air),
            ("Earth", geometry.holes.earth, geometry.pockets.earth))
            println(io, "FEM", medium, "Holes() = {};")
            for tags in holes
                println(io, "FEMHole = newll; Curve Loop(FEMHole) = {", join(tags, ","), "};")
                println(io, "FEM", medium, "Holes() += {FEMHole};")
            end
            numbers = [@sprintf("%.17g", Float64(number)) for number in pockets]
            println(io, "FEM", medium, "Pockets() = {", join(numbers, ", "), "};")
        end
        for (k, tags) in enumerate(cable_curves)
            numbers = [@sprintf("%.17g", Float64(number)) for number in tags]
            println(io, "FEMCableCurves~{", k - 1, "}() = {", join(numbers, ", "), "};")
            name = @sprintf("LCM/cable_contour/%04d/%s", k, system.designs[k].id)
            escaped = replace(name, '\\' => "\\\\", '"' => "\\\"")
            println(io, "Physical Curve(\"", escaped, "\",", TAGS.cable_contour_base + k,
                ") = {", join(tags, ","), "};")
        end
        for (i, surfaces) in enumerate(geometry.surfaces)
            edges = curves(surfaces)
            dimtags = [(1, edge) for edge in edges]
            tags = unique(last.(gmsh.model.getBoundary(dimtags, false, false, false)))
            lists = (("Surfaces", surfaces), ("Curves", edges), ("Points", tags))
            for (name, entries) in lists
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "FEMRegion", name, "~{", i - 1, "}() = {",
                    join(numbers, ", "), "};")
            end
            material = model.material[i]
            conductor = model.materials[material].conductor
            println(io, "FEMPassive~{", i - 1, "} = ", Int(!conductor), ";")
            println(io, "FEMConductor~{", i - 1, "} = ", Int(conductor), ";")
            conductor || continue
            println(io, "FEMConductorMaterial~{", i - 1, "} = ", material, ";")
            owned = [strip for strip in geometry.strips if strip.region == i]
            println(io, "FEMConductorStrips~{", i - 1, "} = ", length(owned), ";")
            if !isempty(owned)
                println(io, "FEMStripDepth~{", i - 1, "} = ",
                    @sprintf("%.17g", first(owned).depth), ";")
            end
            arcs = [turn(geometry.curves[abs(edge)]) for edge in edges]
            thicknesses = map(edges) do edge
                thickness = Inf
                for layer in layers
                    edge in layer.curves || continue
                    thickness = min(thickness, dimension(layer.shape))
                end
                return isfinite(thickness) ? thickness : 0.0
            end
            for (name, entries) in (("ArcLengths", [arc.length for arc in arcs]),
                ("ArcFractions", [arc.fraction for arc in arcs]),
                ("LayerThickness", thicknesses), (
                    "BulkSurfaces", setdiff(surfaces, [strip.surface for strip in owned])))
                numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
                println(io, "FEMConductor", name, "~{", i - 1, "}() = {",
                    join(numbers, ", "), "};")
            end
        end
        for (g, material) in enumerate(model.materials)
            escaped = replace(material.name, '\\' => "\\\\", '"' => "\\\"")
            println(io, "Physical Surface(\"", escaped, "\",", material.tag, ") = {",
                join(material_surfaces[g], ","), "};")
        end
        conductor = Int[]
        passive = Int[]
        for (g, material) in enumerate(model.materials)
            append!(material.conductor ? conductor : passive, material_surfaces[g])
        end
        for (tag, name, entities) in
            ((TAGS.conductors, "LCM/domain/conductors", conductor),
            (TAGS.passive, "LCM/domain/passive_cable_media", passive))
            isempty(entities) && continue
            println(io, "Physical Surface(\"", name, "\",", tag, ") = {",
                join(entities, ","), "};")
        end
        for (t, surfaces) in enumerate(terminal_surfaces)
            groups = (
                ("Surface", TAGS.terminal_base + t,
                    @sprintf("LCM/terminal/%04d/%s", t, model.label[t]), surfaces),
                ("Curve", TAGS.terminal_contour_base + t,
                    @sprintf("LCM/terminal_contour/%04d", t), curves(surfaces)))
            for (entity, tag, name, entities) in groups
                escaped = replace(name, '\\' => "\\\\", '"' => "\\\"")
                println(io, "Physical ", entity, "(\"", escaped, "\",", tag, ") = {",
                    join(entities, ","), "};")
            end
        end
        for (name, entries) in (("X", [p.x for p in geometry.anchors]),
            ("Y", [p.y for p in geometry.anchors]),
            ("MetalDimension", [p.dimension for p in geometry.anchors]))
            numbers = [@sprintf("%.17g", Float64(number)) for number in entries]
            println(io, "FEMReceiver", name, "() = {", join(numbers, ", "), "};")
        end
        return
    end
    rm(raw)
    return
end

"""
$(TYPEDSIGNATURES)

Draw the elliptic curve ``at(a cos θ,b sin θ)`` in \\[m\\] with Gmsh ellipse arcs.
Split its parametric span \\[rad\\] into intervals at most ``π/2``. Reuse the elliptic
arc of the same ellipse and endpoints from `elliptic_arcs`. Return its stored tag with
the requested orientation.

For a nonzero normal offset ``d``, draw a spline through ``p(θ)+dn(θ)``.
Use ``m=max(2,⌈segments⋅|Δα|/(2π)⌉)`` intervals at equal normal angles
``α_j=α_1+jΔα/m`` from the normal ``α_1`` at the piece's start, each the point of greatest
projection on ``α_j``, by [`apsides`](@ref). Angles are in \\[rad\\] and the offset in
\\[m\\].
"""
function draw(geometry::Geometry, arc::EllipticArc)
    tags = Int[]
    n = ceil(Int, abs(arc.span) / (π / 2))
    center = point!(geometry, (arc.at.x, arc.at.y))
    major = point!(geometry, arc.at(arc.a >= arc.b ? (arc.a, 0) : (0, arc.b)))
    for k in 0:(n - 1)
        curve = EllipticArc(
            arc.at, arc.a, arc.b, arc.φ0 + k * arc.span / n, arc.span / n, arc.offset)
        p, q = point!(geometry, curve(0)), point!(geometry, curve(1))
        tag = get!(geometry.elliptic_arcs,
            (key(curve.at.x, curve.at.y, curve.at.φ, curve.a, curve.b, curve.offset),
                minmax(p, q))) do
            start_tag, stop_tag = curve.span > 0 ? (p, q) : (q, p)
            if iszero(curve.offset)
                native = gmsh.model.geo.addEllipseArc(start_tag, center, major, stop_tag)
            else
                α₁ = direction(curve, 0) - sign(curve.span) * π / 2
                Δα = sign(curve.span) * 2π * turn(curve).fraction
                m = max(2, ceil(Int, geometry.segments * turn(curve).fraction))
                points = [p]
                for j in 1:(m - 1)
                    push!(points, point!(geometry, curve(apsides(curve, α₁ + j * Δα / m).near)))
                end
                push!(points, q)
                native = gmsh.model.geo.addSpline(curve.span > 0 ? points : reverse(points))
            end
            geometry.ends[native] = (start_tag, stop_tag)
            geometry.curves[native] = EllipticArc(curve.at, curve.a, curve.b,
                curve.span > 0 ? curve.φ0 : curve.φ0 + curve.span, abs(curve.span), curve.offset)
            return native
        end
        push!(tags, geometry.ends[tag] == (p, q) ? tag : -tag)
    end
    return tags
end

"""
$(TYPEDSIGNATURES)

Return `geometry` unchanged for an elliptic arc.
"""
function cut!(geometry::Geometry, ::EllipticArc)
    return geometry
end
