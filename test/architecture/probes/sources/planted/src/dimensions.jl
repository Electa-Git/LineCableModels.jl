dimension(shape::Rectangle) = 1
dimension(curve::Clothoid) = 1
dimension(shape::AbstractShape) = 1
dimension(::Sector) = 1
dimension(shape::Circle{T}) where {T} = 1
dimension(shape::Union{Annulus, Ellipse}) = 1
