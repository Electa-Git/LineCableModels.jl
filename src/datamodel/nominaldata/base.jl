
# Scalar-type query
Base.eltype(::NominalData{T}) where {T} = T
Base.eltype(::Type{NominalData{T}}) where {T} = T