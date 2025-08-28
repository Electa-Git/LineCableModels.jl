import Base: show

# Pretty printing with uncertainty information if present
function show(io::IO, params::LineParameters{T}) where {T}
    n_cond, _, n_freq = size(params.Z)
    print(io, "LineParameters with $(n_cond) conductors at $(n_freq) frequencies")
    if T <: Complex{Measurement{Float64}}
        print(io, " (with uncertainties)")
    end
end