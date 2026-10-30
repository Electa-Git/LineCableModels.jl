defaults(::Type{Formula{:lossless}}) = (;)

function (::Formula{:lossless})(material::Material{T}, f::T) where {T <: Real}
    return im * 2π * f * vacuum_permittivity * material.eps_r
end
