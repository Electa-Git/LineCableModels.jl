import Base: eltype

eltype(::ConductorGroup{T}) where {T} = T
eltype(::Type{ConductorGroup{T}}) where {T} = T