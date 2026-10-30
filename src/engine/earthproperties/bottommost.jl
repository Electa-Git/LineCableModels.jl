defaults(::Type{Formula{:bottommost}}) = (;)

function (::Formula{:bottommost})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <:
                                                                                     Real}
    return last(layers)
end
