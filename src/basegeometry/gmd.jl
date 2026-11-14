"""
$(TYPEDSIGNATURES)

Return the derivative `(dx/dt, dy/dt)` of the curve position, in \\[m\\],
with respect to the dimensionless parameter `t`.
For a segment this is the endpoint difference. For a circular arc it is the
rotated vector ``r\\,span(-\\sin θ,\\cos θ)``, where ``θ=t\\,span``.
For an elliptic arc it is the rotated vector

```math
span\\left(1+\\frac{d ab}{ρ^3}\\right)(-a\\sin θ,b\\cos θ),
\\qquad ρ=\\sqrt{a^2\\sin^2 θ+b^2\\cos^2 θ},
```

where ``θ=φ0+t\\,span``, ``a,b`` are semi-axes \\[m\\], ``d`` is the normal
offset \\[m\\], and ``span`` and ``φ0`` are angles \\[rad\\].
"""
function derivative(curve::Segment, t::Real)
    return curve.q .- curve.p
end

function derivative(curve::Arc, t::Real)
    θ = t * curve.span
    dx, dy = (-curve.r * curve.span * sin(θ), curve.r * curve.span * cos(θ))
    cosine, sine = cos(curve.at.φ), sin(curve.at.φ)
    return (cosine * dx - sine * dy, sine * dx + cosine * dy)
end

function derivative(curve::EllipticArc, t::Real)
    θ = curve.φ0 + t * curve.span
    ρ = hypot(curve.a * sin(θ), curve.b * cos(θ))
    factor = curve.span * (1 + curve.offset * curve.a * curve.b / ρ^3)
    dx, dy = (-curve.a * sin(θ) * factor, curve.b * cos(θ) * factor)
    cosine, sine = cos(curve.at.φ), sin(curve.at.φ)
    return (cosine * dx - sine * dy, sine * dx + cosine * dy)
end

"""
$(TYPEDSIGNATURES)

Return the geometric mean distance of the regions bounded by `A` and `B`, in \\[m\\].
A region paired with itself gives its geometric mean radius. Grover's definition
[grover1981inductance](@cite) and the double contour identity give

```math
\\log g = \\frac{1}{|A||B|}\\int_A\\int_B \\log|x-y|\\,dA\\,dB
= -\\frac{1}{|A||B|}\\oint_{\\partial A}\\oint_{\\partial B}
H(|x-y|)(n_x\\cdot n_y)\\,ds_x\\,ds_y,\\qquad
H(r)=\\frac{r^2}{4}(\\log r-1).
```

Here ``|A|`` and ``|B|`` are areas in \\[m²\\], ``x,y`` are positions in \\[m\\],
and ``n_x,n_y`` are dimensionless outward unit normals. Distances inside logarithms
are numerical values in meters. A hole's path takes the region's outward normal,
opposite to the right-hand normal of its counterclockwise path.
Quadrature splits at interior nearest projections of endpoints and of each outer quadrature point.

If disjoint disks about the centroids contain the regions, the area average follows

```math
\\log g=\\log|D|+\\operatorname{Re}\\sum_{k=1}^{K}
\\frac{(-1)^{k+1}}{kD^k}\\sum_{j=0}^{k}\\binom{k}{j}
\\mathbb E[\\omega^j]\\mathbb E[(-\\zeta)^{k-j}],
\\qquad \\int_A\\zeta^j\\,dA=\\frac{1}{2i}\\oint_{\\partial A}
\\overline{\\zeta}\\zeta^j\\,d\\zeta.
```

Here ``D`` separates the centroids, and ``\\zeta,\\omega`` measure positions
from the centroids of `A` and `B`, all in \\[m\\]. The contour rule includes signed
hole paths.

The absolute tolerance on ``\\log g`` is ``τ=10^{-6}``.
Three FEM meshes of the bare and plated bars gave an inductance error bar of
0.35 nH/m in about 2900 nH/m. An error ``τ`` moves the inductance by at most
``\\mu_0 τ/(2\\pi) = 2\\times10^{-4}`` nH/m, where ``\\mu_0`` is vacuum permeability.

The calculation selects each regime by its computed remainder or quadrature error estimate.
The series tail uses half this target. The moments use the remainder through
an adaptive vector rule. The normalized moment error ``δ`` satisfies
``(2δ+δ^2)\\sum_{k=1}^K χ^k/k=τ/2``. Here ``χ`` is the dimensionless ratio of the sum
of the disks' radii to the distance between their centroids.
The series applies when its cost estimate does not exceed
the initial nested adaptive rule's cost estimate.

Each curve pair receives ``τ|A||B|/N``, split equally among its pieces.
For two straight pieces with midpoint separation ``ℓ`` and summed half-length
``S<ℓ``, the ellipse semiaxis is ``α=(1+ℓ/S)/2``, and
``R=α+\\sqrt{α^2-1}``. A one-dimensional Gauss rule integrates the outer
curve, with the inner integral in closed form. The expression
``8M R^{-2n}/(1-R^{-1})`` bounds its error, where ``M`` bounds the transformed kernel on
that ellipse.

The bound uses ``|\\log(1+z)| ≤ -\\log(1-|z|)`` for
``|z|<1``, which includes the complex logarithm's argument. Its order is
the least integer that meets the piece budget.
This rule applies when it costs no more than the initial 31-point adaptive rule.

Other pieces try Gauss–Kronrod orders 3, 4, 5, 6, 7, and 15 in that order.
The first five use the initial panels only. A curve pair accepts the first order
whose summed estimates meet its budget. Order 15 allows adaptive subdivision.
The order is a runtime integer shared by the inner and outer rules.

The relative tolerance is zero. The nested rule splits the piece budget equally
between its inner and outer integrals.

When either curve has zero turn fraction, the inner integral follows its
straight line between the callable endpoints.
Its value is the difference of the primitive

```math
F(u,h)=(u^3/24+h^2u/8)\\log(u^2+h^2)-u^3/9-5h^2u/12
+h^3\\arctan(u/h)/6,
```

where ``u`` is distance along the segment from the perpendicular foot and ``h`` is
its distance from the outer point, both in \\[m\\]. The calculation uses
the limits directly at ``h=0`` and ``u=h=0``.

# Notes

If adaptive order 15 exceeds a curve pair's budget, a warning reports the pair's
error estimate and budget.
"""
function gmd(A::Contour{TA}, B::Contour{TB}) where {TA, TB}
    τ = 1e-6
    cfg = (order = 15, rtol = 0.0)
    T = promote_type(TA, TB)
    mass = area(A) * area(B)
    counts = map((A, B)) do contour
        sum(length, (contour.outer, contour.holes...))
    end
    count = prod(counts)
    curved = map((A, B)) do contour
        sum(curve -> iszero(nominal(turn(curve).fraction)) ? 0 : 1,
            (curve for path in (contour.outer, contour.holes...) for curve in path))
    end
    nearcost = 31 * (count - prod(curved)) + 31^2 * prod(curved)
    centers = (centroid(A), centroid(B))
    displacement = complex((centers[2] .- centers[1])...)
    radii = map((A, B), centers) do contour, center
        point = nominal.(center)
        reach = zero(first(point))
        for path in (contour.outer, contour.holes...), curve in path

            outline = map(nominal, curve)
            reach = max(reach, hypot((outline(apsides(outline, point).far) .- point)...))
        end
        reach
    end
    separation = abs(nominal(displacement))
    χ = sum(radii) / separation
    order = 1
    if χ < 1
        while χ^(order + 1) / ((order + 1) * (1 - χ)) > τ / 2 &&
            31 * sum(counts) + order * (order + 1) / 2 <= nearcost
            order += 1
        end
    end
    if χ < 1 && χ^(order + 1) / ((order + 1) * (1 - χ)) <= τ / 2 &&
       31 * sum(counts) + order * (order + 1) / 2 <= nearcost
        degree = order
        harmonic = sum(χ^k / k for k in 1:degree)
        target = (τ / (2 * harmonic)) / (sqrt(1 + τ / (2 * harmonic)) + 1)
        moments = map((A, B), centers, radii, counts) do contour, center, radius, terms
            local values, value
            values = zeros(Complex{T}, degree + 1)
            for (index, path) in enumerate((contour.outer, contour.holes...)), curve in path

                local orientation = index == 1 ? 1 : -1
                # The barriers reduce QuadGK compilation in flatten. On Julia 1.12.7,
                # x86_64 Linux, three fresh processes per variant gave these medians.
                # coaxial_rectangular_strands: first call 267.355 → 11.697 s, warm 158.587 → 160.619 ms.
                # milliken_solid_segments: first call 278.112 → 7.479 s, warm 36.103 → 36.212 ms.
                # Each warm sample is the minimum of ten calls.
                value,
                _ = Base.inferencebarrier(QuadGK.quadgk!)(
                    (powers,
                        u) -> begin
                        point = complex((curve(u) .- center)...)
                        tangent = complex(derivative(curve, u)...)
                        term = orientation * conj(point) * tangent / (2im)
                        scaled = point / radius
                        for k in 0:degree
                            powers[k + 1] = term
                            term *= scaled
                        end
                        powers
                    end, zeros(Complex{T}, degree + 1), 0.0,
                    1.0; order = cfg.order, rtol = cfg.rtol,
                    atol = target * abs(nominal(area(contour))) / terms,
                    norm = x -> maximum(y -> abs(nominal(y)), x))::Tuple{
                    Vector{Complex{T}}, Float64}
                values .+= value
            end
            values ./ area(contour) .* (radius / displacement) .^ (0:degree)
        end
        logarithm = log(abs(displacement))
        for k in 1:order
            value = zero(Complex{T})
            coefficient = 1.0
            for j in 0:k
                value += coefficient * moments[2][j + 1] * moments[1][k - j + 1] *
                         (isodd(k - j) ? -1 : 1)
                coefficient *= (k - j) / (j + 1)
            end
            logarithm += (isodd(k) ? 1 : -1) * real(value) / k
        end
        return exp(logarithm)
    end
    allocated = (QuadGK.alloc_segbuf(Float64, T, Float64; size = 1000),
        QuadGK.alloc_segbuf(Float64, T, Float64; size = 1000))
    atol = τ * abs(nominal(mass)) / count
    total = zero(T)
    for (i, path) in enumerate((A.outer, A.holes...)),
        (j, other) in enumerate((B.outer, B.holes...))

        orientation = (i == 1 ? 1 : -1) * (j == 1 ? 1 : -1)
        for firstcurve in path, secondcurve in other

            let (outercurve,
                    innercurve) = iszero(nominal(turn(firstcurve).fraction)) &&
                                  nominal(turn(secondcurve).fraction) > 0 ?
                                  (secondcurve, firstcurve) : (firstcurve, secondcurve)
                outernominal,
                innernominal = map(nominal, outercurve), map(nominal, innercurve)
                ss,
                tt = map(((outernominal, innernominal),
                    (innernominal, outernominal))) do (curve, othercurve)
                    ts = [0.0, 1.0]
                    for pt in (othercurve(0), othercurve(1))
                        v = apsides(curve, pt).near
                        0 < v < 1 && push!(ts, v)
                    end
                    sort!(unique!(ts))
                end
                subtotal = zero(T)
                for nearorder in (3, 4, 5, 6, 7, 15)
                    subtotal = zero(T)
                    pairerror = 0.0
                    for (u0, u1) in zip(ss[1:(end - 1)], ss[2:end]),
                        (v0, v1) in zip(tt[1:(end - 1)], tt[2:end])

                        p0, p1,
                        pm = outernominal(u0), outernominal(u1), outernominal((u0 + u1) / 2)
                        q0, q1,
                        qm = innernominal(v0), innernominal(v1), innernominal((v0 + v1) / 2)
                        budget = atol / ((length(ss) - 1) * (length(tt) - 1))
                        ℓ = hypot((pm .- qm)...)
                        halfsum = (hypot((p1 .- p0)...) + hypot((q1 .- q0)...)) / 2
                        bound = Inf
                        n = 1
                        if iszero(nominal(turn(outernominal).fraction)) &&
                           iszero(nominal(turn(innernominal).fraction)) &&
                           ℓ > halfsum
                            semiaxis = (1 + ℓ / halfsum) / 2
                            ellipse = semiaxis + sqrt(semiaxis^2 - 1)
                            radius = halfsum * semiaxis
                            magnitude = (ℓ + radius)^2 / 4 *
                                        (abs(log(ℓ)) -
                                         log1p(-radius / ℓ) +
                                         1) *
                                        hypot((p1 .- p0)...) * hypot((q1 .- q0)...)
                            needed = 8 * magnitude <= budget * (1 - inv(ellipse)) ? 1.0 :
                                     log(8 * magnitude / (budget * (1 - inv(ellipse)))) /
                                     (2 * log(ellipse))
                            if needed <= 31
                                n = max(1, ceil(Int, needed))
                                bound = 8 * magnitude * ellipse^(-2n) / (1 - inv(ellipse))
                            end
                        end
                        far = bound <= budget
                        nodes, weights = if far
                            abscissas, coefficients = QuadGK.gauss(n)
                            ((abscissas .+ 1) ./ 2, coefficients ./ 2)
                        else
                            (Float64[], Float64[])
                        end
                        if iszero(nominal(turn(innercurve).fraction))
                            integral = (far ? map : Base.inferencebarrier(quadgk))(
                                u -> begin
                                    x = outercurve(u)
                                    tangent = derivative(outercurve, u)
                                    endpoints = innercurve(1) .- innercurve(0)
                                    innerlength = hypot(endpoints...)
                                    ex, ey = endpoints ./ innerlength
                                    start = innercurve(0) .+ v0 .* endpoints .- x
                                    h = abs(start[1] * ey - start[2] * ex)
                                    along0 = start[1] * ex + start[2] * ey
                                    along1 = along0 + (v1 - v0) * innerlength
                                    values = map((along0, along1)) do z
                                        r2 = z^2 + h^2
                                        if iszero(nominal(r2))
                                            zero(r2 * z)
                                        elseif iszero(nominal(h))
                                            z^3 * (log(r2) / 24 - 1 / 9)
                                        else
                                            (z^3 / 24 + h^2 * z / 8) * log(r2) -
                                            z^3 / 9 -
                                            5h^2 * z / 12 + h^3 / 6 * atan(z / h)
                                        end
                                    end
                                    (values[2] - values[1]) * dot(tangent, (ex, ey))
                                end,
                                (far ? (u0 .+ (u1 - u0) .* nodes,) : (u0, u1))...;
                                (far ? (;) :
                                 (order = nearorder,
                                    maxevals = nearorder == 15 ? 10^7 : 0,
                                    rtol = cfg.rtol, atol = budget,
                                    segbuf = allocated[2], norm = abs ∘ nominal))...)
                            if far
                                subtotal += (u1 - u0) * dot(weights, integral)
                                gapestimate = bound
                            else
                                value, gap = integral::Tuple{T, Float64}
                                subtotal += value
                                gapestimate = gap
                            end
                        else
                            innererror = [0.0]
                            value,
                            gap = Base.inferencebarrier(quadgk)(
                                u -> begin
                                    vnear = apsides(innernominal, outernominal(u)).near
                                    ends = v0 < vnear < v1 ? (v0, vnear, v1) : (v0, v1)
                                    innervalue,
                                    innergap = Base.inferencebarrier(quadgk)(
                                        v -> begin
                                            x, y = outercurve(u), innercurve(v)
                                            r2 = (x[1] - y[1])^2 + (x[2] - y[2])^2
                                            H = iszero(nominal(r2)) ? zero(r2) :
                                                r2 * (log(r2) / 2 - 1) / 4
                                            H *
                                            dot(derivative(outercurve, u), derivative(innercurve, v))
                                        end,
                                        ends...;
                                        order = nearorder,
                                        maxevals = nearorder == 15 ? 10^7 : 0,
                                        rtol = cfg.rtol,
                                        atol = budget / (2 * (u1 - u0)), segbuf = allocated[1],
                                        norm = abs ∘ nominal)::Tuple{T, Float64}
                                    innererror[1] = max(innererror[1], innergap)
                                    innervalue
                                end,
                                u0,
                                u1;
                                order = nearorder,
                                maxevals = nearorder == 15 ? 10^7 : 0,
                                rtol = cfg.rtol,
                                atol = budget / 2, segbuf = allocated[2],
                                norm = abs ∘ nominal)::Tuple{T, Float64}
                            subtotal += value
                            gapestimate = gap + (u1 - u0) * innererror[1]
                        end
                        pairerror += gapestimate
                    end
                    if pairerror <= atol
                        break
                    end
                    if nearorder == 15
                        @warn "QuadGK error estimate exceeds the curve pair budget" estimate=pairerror budget=atol
                    end
                end
                total += orientation * subtotal
            end
        end
    end
    return exp(-total / mass)
end
