module ValuesProbePlanted
struct Loose
    x::Any
end
struct Wide
    x::Real
end
struct Open{T}
    x::T
end
struct Plain
    x::Float64
end
Plain(x, y) = Plain(x + y)
struct Nested
    x::Vector{Real}
    y::Dict{Symbol, Any}
    Nested(x::Vector{Real}, y::Dict{Symbol, Any}) = new(x, y)
end
struct Either
    x::Union{Float64, Int}
    Either(x::Union{Float64, Int}) = new(x)
end
struct Blend
    x::Vector{Union{Integer, Float64}}
    Blend(x::Vector{Union{Integer, Float64}}) = new(x)
end
function boxed_keyword(values; start = 0)
    total = start
    foreach(value -> (total += value), values)
    return total
end
function boxed(values)
    total = 0
    foreach(value -> (total += value), values)
    return total
end
end
