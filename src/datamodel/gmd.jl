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
    C, D=Contour(A), Contour(B)
    a=isempty(C.holes) ? zero(first(C.outer).r) : first(first(C.holes)).r
    b=first(C.outer).r
    ai=isempty(D.holes) ? zero(first(D.outer).r) : first(first(D.holes)).r
    bi=first(D.outer).r
    d=hypot(A.at.x-B.at.x, A.at.y-B.at.y)
    tol=tolerance(max(nominal(b), nominal(bi)))
    if nominal(d)<=tol && abs(nominal(a-ai))<=tol && abs(nominal(b-bi))<=tol
        if iszero(nominal(a))
            r = b
            return exp(log(r)-1/4)
        end
        return exp(log(b)-(a^4/(b^2-a^2)^2*log(b/a)-(3a^2-b^2)/(4*(b^2-a^2))))
    elseif nominal(d+bi)<=nominal(a)+tol
        return exp((b^2*log(b)-a^2*log(a))/(b^2-a^2)-1/2)
    elseif nominal(d+b)<=nominal(ai)+tol
        return exp((bi^2*log(bi)-ai^2*log(ai))/(bi^2-ai^2)-1/2)
    elseif nominal(d)>=nominal(b+bi)-tol
        return d
    end
    gmd(C, D)
end
function gmd(A::Union{Annulus, Circle}, B::AbstractShape)
    C, D=Contour(A), Contour(B)
    isempty(C.holes) && return gmd(C, D)
    a=first(first(C.holes)).r
    b=first(C.outer).r
    center=(nominal(A.at.x), nominal(A.at.y))
    far=maximum(D.outer) do curve
        outline=map(nominal, curve)
        hypot((outline(apsides(outline, center).far) .- center)...)
    end
    if far<=nominal(a)+tolerance(nominal(b))
        return exp((b^2*log(b)-a^2*log(a))/(b^2-a^2)-1/2)
    end
    gmd(C, D)
end
gmd(A::AbstractShape, B::Union{Annulus, Circle}) = gmd(B, A)
gmd(A::AbstractShape, B::AbstractShape) = gmd(Contour(A), Contour(B))
