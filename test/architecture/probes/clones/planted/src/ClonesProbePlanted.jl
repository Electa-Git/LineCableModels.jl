module ClonesProbePlanted

module Commons
abstract type AbstractFormulation end
end

# One computation written twice, its body at the floor: caught in both directions.
function spread_left(values, weights)
    total = sum(weights)
    mean = sum(values .* weights) / total
    deviation = sqrt(sum(weights .* abs2.(values .- mean)) / total)
    skew = sum(weights .* (values .- mean) .^ 3) / (total * deviation^3)
    low, high = mean - 2 * deviation, mean + 2 * deviation
    return (low, high, skew, high - low)
end

function spread_right(values, weights)
    total = sum(weights)
    mean = sum(values .* weights) / total
    deviation = sqrt(sum(weights .* abs2.(values .- mean)) / total)
    skew = sum(weights .* (values .- mean) .^ 3) / (total * deviation^3)
    low, high = mean - 2 * deviation, mean + 2 * deviation
    return (low, high, skew, high - low)
end

struct Box{T <: Real}
    longest::T
    middle::T
    shortest::T
    volume::T
    area::T
    edges::T
    aspect::T
end

struct Crate{T <: Real}
    longest::T
    middle::T
    shortest::T
    volume::T
    area::T
    edges::T
    aspect::T
end

# The constructors of two declared types share one body, which the guard admits.
function Box(; width = 1.0, height = 2.0, depth = 3.0, scale = 1.0)
    w, h, d = promote(scale * width, scale * height, scale * depth)
    w > 0 && h > 0 && d > 0 || throw(DomainError((w, h, d), "a side is not positive"))
    longest, shortest = max(w, h, d), min(w, h, d)
    return Box(longest, w + h + d - longest - shortest, shortest, w * h * d,
        2 * (w * h + h * d + d * w), 4 * (w + h + d), longest / shortest)
end

function Crate(; width = 1.0, height = 2.0, depth = 3.0, scale = 1.0)
    w, h, d = promote(scale * width, scale * height, scale * depth)
    w > 0 && h > 0 && d > 0 || throw(DomainError((w, h, d), "a side is not positive"))
    longest, shortest = max(w, h, d), min(w, h, d)
    return Crate(longest, w + h + d - longest - shortest, shortest, w * h * d,
        2 * (w * h + h * d + d * w), 4 * (w + h + d), longest / shortest)
end

# The body of a constructor in a function that constructs nothing: caught against both.
function sorted_sides(; width = 1.0, height = 2.0, depth = 3.0, scale = 1.0)
    w, h, d = promote(scale * width, scale * height, scale * depth)
    w > 0 && h > 0 && d > 0 || throw(DomainError((w, h, d), "a side is not positive"))
    longest, shortest = max(w, h, d), min(w, h, d)
    return tuple(longest, w + h + d - longest - shortest, shortest, w * h * d,
        2 * (w * h + h * d + d * w), 4 * (w + h + d), longest / shortest)
end

struct Formula{ID} <: Commons.AbstractFormulation end

# The methods of a formula share one body and the function inside it, which the guard admits.
function (formula::Formula{:left})(x)
    integrand = function (u)
        z = x * cos(u) + sqrt(abs(1 - x^2)) * sin(u)
        damping = exp(-z^2 / 2) / (1 + z^2)
        phase = 3 * u - x * sin(2 * u) + atan(z, 1 + x^2)
        return damping * (cos(phase) - x * sin(phase) / (2 + cos(u)))
    end
    return sum(integrand(2π * (k - 1 / 2) / 64) for k in 1:64) * 2π / 64
end

function (formula::Formula{:right})(x)
    integrand = function (u)
        z = x * cos(u) + sqrt(abs(1 - x^2)) * sin(u)
        damping = exp(-z^2 / 2) / (1 + z^2)
        phase = 3 * u - x * sin(2 * u) + atan(z, 1 + x^2)
        return damping * (cos(phase) - x * sin(phase) / (2 + cos(u)))
    end
    return sum(integrand(2π * (k - 1 / 2) / 64) for k in 1:64) * 2π / 64
end

end
