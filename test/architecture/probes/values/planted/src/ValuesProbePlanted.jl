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
function boxed(values)
    total = 0
    foreach(value -> (total += value), values)
    return total
end
end
