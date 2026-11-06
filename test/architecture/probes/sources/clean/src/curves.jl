function partition(curves)
    return ([curve for curve in curves if curve isa Segment],
        [curve for curve in curves if curve isa Arc],
        [curve for curve in curves if curve isa EllipticArc])
end
