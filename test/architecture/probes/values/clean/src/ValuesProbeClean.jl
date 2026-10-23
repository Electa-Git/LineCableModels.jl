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
Tight(x::Integer, y::Integer) = Tight(x + y)
scaled(values, k) = map(value -> k * value, values)
function summed(values)
    total = sum(values)
    return total
end
end
