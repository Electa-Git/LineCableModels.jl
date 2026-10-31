defaults(::Type{Formula{:lossy}}) = (;)

function (::Formula{:lossy})(material::Material{T}, f::T) where {T <: Real}
    return complex(
        inv(material.rho) +
        2π * f * vacuum_permittivity * material.eps_r * material.tan_delta,
        2π * f * vacuum_permittivity * material.eps_r)
end
