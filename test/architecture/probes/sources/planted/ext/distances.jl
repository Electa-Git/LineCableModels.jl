function section(x::Annulus, shape::AbstractShape, controls::NamedTuple)
    w = x.ro - x.ri
    return w
end
