module ValuesProbeClean
struct Tight{T <: Real}
    x::T
end
struct Named
    name::String
    Named(name::String) = new(name)
end
struct Holder{T <: Real}
    values::Vector{T}
end
struct Grid
    cells::Vector{Float64}
    Grid(cells::Vector{Float64}) = new(cells)
end
Tight(x::Integer, y::Integer) = Tight(x + y)
function weighted(values; weight = 2)
    total = sum(values)
    return weight * total
end
scaled(values, k) = map(value -> k * value, values)
function summed(values)
    total = sum(values)
    return total
end
end
