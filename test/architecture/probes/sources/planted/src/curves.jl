function partition(curves)
    return [curve for curve in curves if curve isa Union{Segment, Arc}]
end
partition(curves, selected) = [curve for curve in curves if curve isa EllipticArc]
# Every member tested in one method: an exhaustive switch is still a switch.
function partition(curves, selected, rest)
    return ([curve for curve in curves if curve isa Segment],
        [curve for curve in curves if curve isa Arc],
        [curve for curve in curves if curve isa EllipticArc])
end
# The predicates of two selection calls test two members.
function selected(curves)
    return (filter(curve -> curve isa Segment, curves),
        findall(curve -> curve isa Arc, curves))
end
# One member tested beside its negation.
function split(curves)
    return (findall(curve -> curve isa Arc, curves),
        findall(other -> !(other isa Arc), curves))
end
