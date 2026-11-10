defaults(::Type{Formula{:bottommost}}) = (;)

function (::Formula{:bottommost})(layers::Tuple{Vararg{EarthLayer{T}}}, f::T) where {T <:
                                                                                     Real}
    deepest = last(layers)
    return EquivalentEarthLayer{T}(deepest.rho, deepest.eps_r, deepest.mu_r)
end
