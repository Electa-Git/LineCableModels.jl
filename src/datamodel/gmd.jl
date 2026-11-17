"""
$(TYPEDSIGNATURES)

Return the geometric mean distance of `A` and `B` in \\[m\\].
The closed forms for [`Circle`](@ref) and [`Annulus`](@ref) follow Maxwell,
§692 [maxwell1881treatise](@cite). For a disk of radius ``r`` and an annulus
of inner radius ``a`` and outer radius ``b``, their own terms are

```math
\\log g_\\mathrm{disk} = \\log r - \\frac14,\\qquad
\\log g_\\mathrm{annulus} = \\log b -
\\left[\\frac{a^4\\log(b/a)}{(b^2-a^2)^2}-\\frac{3a^2-b^2}{4(b^2-a^2)}\\right].
```

For ``s = (b-a)/(b+a)`` with `nominal(s) ≤ 1/3`, the annulus term uses the equivalent series:

```math
\\delta = \\log b - \\log g = \\frac18 \\sum_{n=3}^{\\infty} c_n s^{n-2},
\\qquad c_3 = \\frac{16}{3},
```

```math
c_n = \\begin{cases}
-4\\left(\\frac{1}{n-1}+\\frac{1}{n-3}\\right), & n \\ge 4 \\text{ even}, \\\\
\\frac{1}{n}+\\frac{6}{n-2}+\\frac{1}{n-4}, & n \\ge 5 \\text{ odd}.
\\end{cases}
```

The partial sum `δ` accumulates terms ``\\delta_n = c_n s_n/8``, with ``s_n = s^{n-2}``.
After `abs(nominal(δ_n)) < eps(nominal(δ))`, summation includes two further terms.
[`gmd`](@ref) evaluates the closed form for larger `nominal(s)`.

For an annulus and any region contained in its hole,

```math
\\log g = \\frac{b^2\\log b-a^2\\log a}{b^2-a^2}-\\frac12.
```

For two disjoint disks or annuli, ``g`` is the distance between their centers.
Radii and distances are in \\[m\\]. The own terms apply to coincident regions.
The disk expression also holds for an annulus with zero inner radius.
For other pairs, [`gmd`](@ref) evaluates the two contours.
"""
function gmd(A::Union{Annulus, Circle}, B::Union{Annulus, Circle})
    a, b = r_in(A), r_ex(A)
    ai, bi = r_in(B), r_ex(B)
    d=hypot(A.at.x-B.at.x, A.at.y-B.at.y)
    tol=tolerance(max(nominal(b), nominal(bi)))
    if nominal(d)<=tol && abs(nominal(a-ai))<=tol && abs(nominal(b-bi))<=tol
        if iszero(nominal(a))
            r = b
            return exp(log(r)-1/4)
        end
        s = (b - a) / (b + a)
        if nominal(s) <= 1 / 3
            s_n = s
            δ_n = (2 / 3) * s_n
            δ = δ_n
            n = 3
            remaining = -1
            while remaining != 0
                n += 1
                s_n *= s
                c_n = iseven(n) ? -4 * (1 / (n - 1) + 1 / (n - 3)) :
                      1 / n + 6 / (n - 2) + 1 / (n - 4)
                δ_n = (c_n / 8) * s_n
                δ += δ_n
                if remaining < 0 && abs(nominal(δ_n)) < eps(nominal(δ))
                    remaining = 2
                elseif remaining >= 0
                    remaining -= 1
                end
            end
            return exp(log(b) - δ)
        end
        return exp(log(b)-(a^4/(b^2-a^2)^2*log(b/a)-(3a^2-b^2)/(4*(b^2-a^2))))
    elseif nominal(d+bi)<=nominal(a)+tol
        return exp((b^2*log(b)-a^2*log(a))/(b^2-a^2)-1/2)
    elseif nominal(d+b)<=nominal(ai)+tol
        return exp((bi^2*log(bi)-ai^2*log(ai))/(bi^2-ai^2)-1/2)
    elseif nominal(d)>=nominal(b+bi)-tol
        return d
    end
    gmd(Contour(A), Contour(B))
end
function gmd(A::Union{Annulus, Circle}, B::AbstractShape)
    a, b = r_in(A), r_ex(A)
    D = Contour(B)
    iszero(nominal(a)) && return gmd(Contour(A), D)
    center=(nominal(A.at.x), nominal(A.at.y))
    far=maximum(D.outer) do curve
        outline=map(nominal, curve)
        hypot((outline(apsides(outline, center).far) .- center)...)
    end
    if far<=nominal(a)+tolerance(nominal(b))
        return exp((b^2*log(b)-a^2*log(a))/(b^2-a^2)-1/2)
    end
    gmd(Contour(A), D)
end
gmd(A::AbstractShape, B::Union{Annulus, Circle}) = gmd(B, A)
gmd(A::AbstractShape, B::AbstractShape) = gmd(Contour(A), Contour(B))

"""
$(TYPEDSIGNATURES)

Return the geometric mean radius ``g`` \\[m\\] of a conductor of cross-section `S` and
relative permeability `μ` \\[dimensionless\\]. A circle or an annulus of outer radius ``b``
\\[m\\] has

```math
\\log g = \\log b - \\mu (\\log b - \\log \\operatorname{gmd}(S,S)).
```

The difference in parentheses is the internal term of the tube, equal to ``1/4`` for a disk.
Any other cross-section takes the disk's convention for the permeability:

```math
\\log g = \\log \\operatorname{gmd}(S,S) - \\frac{\\mu - 1}{4}.
```

The radius depends on the cross-section alone, wherever it lies. At ``\\mu = 1`` it is the
[`gmd`](@ref) of `S` with itself.
"""
function gmr(S::Union{Circle, Annulus}, μ::Real)
    b = r_ex(S)
    return exp(log(b) - μ * (log(b) - log(gmd(S, S))))
end
gmr(S::AbstractShape, μ::Real) = exp(log(gmd(S, S)) - (μ - 1) / 4)
