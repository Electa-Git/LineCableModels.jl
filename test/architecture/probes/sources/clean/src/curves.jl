# One member per method, by a comprehension and by the predicates of selection calls.
segments(curves) = [curve for curve in curves if curve isa Segment]
arcs(curves) = filter(curve -> curve isa Arc, curves)
function located(curves)
    return (findfirst(curve -> curve isa EllipticArc, curves),
        findlast(curve -> curve isa EllipticArc, curves))
end
