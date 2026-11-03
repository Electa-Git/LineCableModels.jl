defaults(::Type{Formula{:lossless}}) = (;)

function (::Formula{:lossless})(material::Material{T}, f::T) where {T <: Real}
    return im * 2π * f * ε₀ * material.eps_r
end
