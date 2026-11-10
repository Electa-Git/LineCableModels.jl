module ClonesProbeClean

# One computation written twice, its body one token under the floor: not a reference.
function bounds_left(values)
    low, high = extrema(values)
    middle = (low + high) / 2
    half = (high - low) / 2
    below = count(<(middle), values)
    share = below / length(values)
    spread = sum(abs2, values .- middle) / length(values)
    return (middle - 1.5 * half, middle + 1.5 * half,
        below, share, 1 - share, half / middle, spread)
end

function bounds_right(values)
    low, high = extrema(values)
    middle = (low + high) / 2
    half = (high - low) / 2
    below = count(<(middle), values)
    share = below / length(values)
    spread = sum(abs2, values .- middle) / length(values)
    return (middle - 1.5 * half, middle + 1.5 * half,
        below, share, 1 - share, half / middle, spread)
end

# A computation above the floor, written once.
function simpson(f, a, b, n)
    h = (b - a) / (2 * n)
    odd = sum(f(a + (2 * k - 1) * h) for k in 1:n)
    even = sum(f(a + 2 * k * h) for k in 1:(n - 1))
    estimate = h / 3 * (f(a) + 4 * odd + 2 * even + f(b))
    trapezoid = h * (f(a) / 2 + odd + even + f(b) / 2)
    return estimate, abs(estimate - trapezoid)
end

end
