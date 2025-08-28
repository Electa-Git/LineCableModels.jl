import Base: eltype

eltype(::InsulatorGroup{T}) where {T} = T
eltype(::Type{InsulatorGroup{T}}) where {T} = T