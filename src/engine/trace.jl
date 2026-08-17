"""
$(TYPEDEF)

Inspection result containing the public line parameters and selected completed
primitive matrices from one EMT calculation.

$(TYPEDFIELDS)
"""
struct EMTTrace{R, T <: Real}
    result::R
    frequencies::Vector{T}
    phase_map::Vector{Int}
    cable_map::Vector{Int}
    Zin::Array{Complex{T}, 3}
    Pin::Array{Complex{T}, 3}
    Zg::Array{Complex{T}, 3}
    Pg::Array{Complex{T}, 3}
    Z::Array{Complex{T}, 3}
    P::Array{Complex{T}, 3}
end

Base.eltype(::EMTTrace{<:Any, T}) where {T} = T

function Base.show(io::IO, trace::EMTTrace)
    print(io, "EMTTrace(", length(trace.frequencies), " frequencies, ",
        length(trace.phase_map), " primitive conductors)")
end
