"""
$(TYPEDEF)

Inspection result containing the public line parameters and selected completed
primitive matrices from one analytical line-parameter calculation.

$(TYPEDFIELDS)
"""
struct LineParametersTrace{R <: LineParameters, T <: Real} <: AbstractProblemResult
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

Base.eltype(::LineParametersTrace{<:LineParameters, T}) where {T} = T

basis(trace::LineParametersTrace) = basis(trace.result)
observe(trace::LineParametersTrace, selector, args...) = observe(trace.result, selector, args...)
observables(::Type{<:LineParametersTrace}) = observables(LineParameters)

function Base.show(io::IO, trace::LineParametersTrace)
    print(io, "LineParametersTrace(", length(trace.frequencies), " frequencies, ",
        length(trace.phase_map), " primitive conductors)")
end
