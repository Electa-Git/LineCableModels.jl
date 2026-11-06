function partition(curves)
    return [curve for curve in curves if curve isa Union{Segment, Arc}]
end
partition(curves, selected) = [curve for curve in curves if curve isa EllipticArc]
