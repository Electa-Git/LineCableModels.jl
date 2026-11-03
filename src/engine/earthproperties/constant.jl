defaults(::Type{Formula{:constant}}) = (;)

function (::Formula{:constant})(layer::EarthLayer{T}, f::T) where {T <: Real}
    return layer
end
