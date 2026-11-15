function cleared(A, B)
    for image in (A, B)
        fill!(image, 0)
    end
end
function positive(r_in, r_ex)
    for r in (r_in, r_ex)
        r >= 0 || throw(DomainError(r))
        iszero(r) && return false
        r == 0 && return false
    end
    return true
end
function copied(A, B, reflections)
    for (field, image) in enumerate((A, B))
        copyto!(image, view(reflections,:,:,field))
    end
end
truncated(x, order) = [x^n for n in 0:2 if n <= order]
function skipped(xs)
    for (i, x) in enumerate(xs)
        i == 1 && continue
        x > 0 || continue
    end
end
function diagonal(M, n)
    for i in 1:n, j in 1:n

        i == j && continue
        M[i, j] = 0
    end
end
function oriented(contour, d)
    for (i, path) in enumerate((contour.outer, contour.holes...))
        i == 1 || reverse!(path)
    end
end
function paired(A, B, xs)
    for (image, x) in zip((A, B), xs)
        x > 0 && fill!(image, x)
        isempty(image) && continue
    end
end
picked(xs) = [xs[pass] for pass in 1:2]
const WAVES = (sin, cos)
indexed(x) = [(sin, cos)[pass](x) + WAVES[pass](x) for pass in 1:2]
function renewed(xs)
    for pass in 1:2
        pass = xs[pass]
        pass == 1 && return true
    end
    return false
end
shadowed(xs) = [count(k -> k == 1, xs) for k in 1:2]
function applied(f, g, A, B)
    for (formula, image) in ((f, A), (g, B))
        fill!(image, formula(1))
    end
end
function limited(xs, limit)
    for (k, scale) in ((1, 2.0), (2, 0.5))
        scale * xs[k] > limit && return true
        k <= length(xs) || continue
    end
    return false
end
term(n, x) = n == 0 ? log(x) : x^n / n
series(x, order) = term(0, x) + sum(n -> term(n, x), 1:order)
scaled(x; factor = 2) = factor * x
doubled(x) = scaled(x; factor = 4)
clipped(x, bound) = x > bound ? bound : x
capped(x) = clipped(x, 1)
halves(xs) = map(x -> x > 0 ? x / 2 : x, xs)
signs(a, b) = map(x -> x > 0, (a, b))
tripled(f) = ntuple(k -> f(k), 3)
function grown(contour)
    paths = [contour.outer]
    append!(paths, contour.holes)
    for (i, path) in enumerate(paths)
        i == 1 || reverse!(path)
    end
end
countdown(x, n) = (n -= 1; n == 0 ? x : -x)
started(x) = countdown(x, 1)
