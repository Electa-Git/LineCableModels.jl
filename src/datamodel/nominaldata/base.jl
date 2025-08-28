import Base: eltype

# Scalar-type query
eltype(::NominalData{T}) where {T} = T
eltype(::Type{NominalData{T}}) where {T} = T