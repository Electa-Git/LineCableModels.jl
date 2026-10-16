# Integrated finite-face charge element for the explicit boundary shunt model.
# Coordinates are in m. Green functions use charge/(2pi*epsilon0).
# The buffers of the boundary shunt formula hold the dense arrays of one solve.

"""
$(TYPEDSIGNATURES)

Evaluate the Fourier coefficients of the electrostatic Green function of the concentric
dielectric layers of `g`, grounded at the inner anchor and at the reference, for a target and
a source in any two layers. The radial basis ``r^{\\pm m}``, with ``1,\\log r`` for the zero
mode, is the classical cylindrical Laplace basis. See [Schelkunoff1934](@cite),
Eqs. (122)-(124), p. 573, in the small-radius-to-wavelength limit of cylindrical fields. Here
it is used for the electrostatic potential. The layered-kernel strategy is also supported by
[Bernal1997](@cite), Sections II-III, whose planar Prony/matrix-pencil kernel is not used here.

With ``t=\\log r``, layer ``k`` between the radii ``\\rho_{k-1}<\\rho_k`` with relative
permittivity ``\\varepsilon_k`` satisfies ``u''-m^2u=0`` for the mode ``m``. The solution
``u_L`` vanishes at the inner anchor and ``u_R`` at the reference, and both keep ``u`` and
``\\varepsilon u'`` continuous at every interface. In layer ``k``,
``u_L\\propto (r/\\rho_{k-1})^m+\\Gamma^L_k(\\rho_{k-1}/r)^m`` and
``u_R\\propto(\\rho_k/r)^m+\\Gamma^R_k(r/\\rho_k)^m``, with ``\\Gamma^L_1=\\Gamma^R_L=-1`` at
grounded metal. Matching at ``\\rho_k`` gives

```math
\\Gamma^L_{k+1}=\\frac{\\varepsilon_{k+1}(1+\\gamma)-\\varepsilon_k(1-\\gamma)}
{\\varepsilon_{k+1}(1+\\gamma)+\\varepsilon_k(1-\\gamma)},\\qquad
\\gamma=\\Gamma^L_k(\\rho_{k-1}/\\rho_k)^{2m},
```

and the mirror relation for ``\\Gamma^R_k``, the reflection form of the load transformation of
[Sunde1968](@cite), Section 1.6, Eqs. (1.39)-(1.42), p. 15. Every ``|\\Gamma|\\le 1``, so no
mode overflows. With the invariant Wronskian
``W_m=\\varepsilon(u_L'u_R-u_Lu_R')``, the mode of the Green function is
``H_m=u_L(r_<)u_R(r_>)/W_m``, and the kernel is ``H_0+2\\sum_m H_m\\cos m\\Delta\\theta``. For
the inner point in layer ``p`` and the outer point in layer ``q``,

```math
2H_m=\\frac{\\prod_{k=p}^{q-1}T_k}{\\varepsilon_p m(1-\\Gamma^L_p\\Gamma^R_p
(\\rho_{p-1}/\\rho_p)^{2m})}\\left(\\frac{r_<}{r_>}\\right)^m
\\left(1+\\Gamma^L_p\\frac{\\rho_{p-1}^{2m}}{r_<^{2m}}\\right)
\\left(1+\\Gamma^R_q\\frac{r_>^{2m}}{\\rho_q^{2m}}\\right),\\qquad
T_k=\\frac{2\\varepsilon_k}{\\varepsilon_k(1+\\gamma)+\\varepsilon_{k+1}(1-\\gamma)},
```

with ``\\gamma=\\Gamma^R_{k+1}(\\rho_k/\\rho_{k+1})^{2m}``. The zero mode is

```math
H_0(r,s)=\\frac{R_<(R_{\\mathrm{total}}-R_>)}{R_{\\mathrm{total}}},\\qquad
R(r)=\\int_{\\rho_0}^r\\frac{dr'}{\\varepsilon(r')r'},
```

where ``R_<`` and ``R_>`` are the smaller and larger of ``R(r),R(s)``. The log-radius
normalization follows the coaxial capacitance ``C=2\\pi\\varepsilon/\\log(r_o/r_i)`` \\[F/m\\].
See [Sunde1968](@cite), Section 1.5, Eq. (1.20), p. 11. The kernel is reciprocal before
discretization.

The logarithms are extracted before the modes are truncated: the direct logarithm in the
source's layer and the transmitted one in another layer, with the static coefficient
``\\prod_k 2\\varepsilon_k/(\\varepsilon_k+\\varepsilon_{k+1})/\\varepsilon_p``, and, for a
target and a source in one layer, its two nearest images with the static contrasts
``(\\varepsilon_p-\\varepsilon_j)/(\\varepsilon_p+\\varepsilon_j)`` of the adjacent layers,
also used in [Campione2018](@cite), Section 2, Eq. (3). That paper uses a local planar
approximation to the braid geometry. It supports these dielectric image factors, not this
concentric Green function or its radial recursion. A source on an interface sees the mean of its
two permittivities. Only the remaining reflections and transmissions are
truncated.

# Arguments

- `g`: the dielectric layers' radii \\[m\\] and positive relative permittivities
  \\[dimensionless\\].
- `modes`: number of retained positive Fourier modes.

# Returns

- `radii`: the interfaces ``\\rho_0,\\ldots,\\rho_L`` \\[m\\].
- `epsilon`: the layers' relative permittivities.
- `R`: ``R`` at each interface and `Rtotal` at the reference \\[dimensionless\\].
- `inner`, `outer`: the static image contrasts of each layer \\[dimensionless\\].
- `direct`: the static coefficient of the extracted logarithm, by target and source layer.
- `modal`: the remaining coefficients, by mode, by the target's inner or outer factor
  ``(\\rho_{i-1}/\\bar z)^m`` or ``(z/\\rho_i)^m``, by the source's factor of the same kind,
  and by target and source layer. The kernel multiplies charge per unit length divided by
  ``2\\pi\\varepsilon_0`` \\[V\\] to give potential \\[V\\].
"""
function _shunt_kernel_coefficients(g, modes)
    layers = g.layers
    L = length(layers)
    radii = [first(layers).ri; [layer.ro for layer in layers]]
    epsilon = [layer.epsilon for layer in layers]
    R = zeros(L + 1)
    for (k, layer) in pairs(layers)
        R[k + 1] = R[k] + log(layer.ro/layer.ri)/layer.epsilon
    end
    inner = [k == 1 ? -1.0 : (epsilon[k]-epsilon[k-1])/(epsilon[k]+epsilon[k-1]) for k in 1:L]
    outer = [k == L ? -1.0 : (epsilon[k]-epsilon[k+1])/(epsilon[k]+epsilon[k+1]) for k in 1:L]
    direct = [prod(k -> 2epsilon[k]/(epsilon[k]+epsilon[k+1]), min(i, j):(max(i, j)-1);
                  init = 1.0)/epsilon[min(i, j)] for i in 1:L, j in 1:L]
    modal = zeros(modes, 2, 2, L, L)
    left, right, transmission = zeros(L), zeros(L), zeros(L)
    for m in 1:modes
        left[1] = -1
        for k in 2:L
            gamma = left[k-1]*(radii[k-1]/radii[k])^(2m)
            left[k] = (epsilon[k]*(1+gamma)-epsilon[k-1]*(1-gamma)) /
                      (epsilon[k]*(1+gamma)+epsilon[k-1]*(1-gamma))
        end
        right[L] = -1
        for k in (L - 1):-1:1
            gamma = right[k+1]*(radii[k+1]/radii[k+2])^(2m)
            denominator = epsilon[k]*(1+gamma)+epsilon[k+1]*(1-gamma)
            right[k] = (epsilon[k]*(1+gamma)-epsilon[k+1]*(1-gamma))/denominator
            transmission[k] = 2epsilon[k]/denominator
        end
        for j in 1:L, i in 1:L
            p, q = minmax(i, j)
            X = left[p]*right[p]*(radii[p]/radii[p+1])^(2m)
            if i == j
                weight = epsilon[p]*m
                modal[m, 1, 1, i, j] = (left[p]/(1-X)-inner[p])/weight
                modal[m, 2, 2, i, j] = (right[p]/(1-X)-outer[p])/weight
                cross = left[p]*right[p]*(radii[p]/radii[p+1])^m/((1-X)*weight)
                modal[m, 1, 2, i, j] = cross
                modal[m, 2, 1, i, j] = cross
            else
                scale = prod(view(transmission, p:(q - 1)))/(epsilon[p]*m*(1-X))
                transmitted = (scale-direct[i, j]/m)*(radii[p+1]/radii[q])^m
                outward = scale*right[q]*(radii[p+1]/radii[q+1])^m
                inward = scale*left[p]*(radii[p]/radii[q])^m
                both = scale*left[p]*right[q]*(radii[p]/radii[q+1])^m
                # The point in layer p takes its outer factor and the point in layer q its
                # inner factor in the transmitted term.
                near, far = i < j ? (2, 1) : (1, 2)
                modal[m, near, far, i, j] = transmitted
                modal[m, 2, 2, i, j] = outward
                modal[m, 1, 1, i, j] = inward
                modal[m, far, near, i, j] = both
            end
        end
    end
    return (; radii, epsilon, R, Rtotal = R[L + 1], inner, outer, direct, modal)
end

# The zero-mode, extracted logarithmic and nearest-image terms of the Green function for a
# target `z` in layer `i` and a source in layer `j`. The remaining Fourier modes enter by the
# batched products of `_shunt_kernel_matrix!`. With `split_images`, the logarithms of the
# distances to the source and its images, which the tape faces integrate, are left out.
function _shunt_kernel(z, i, source, j, k; split_images = false)
    r, s = abs(z), abs(source)
    Rr = k.R[i] + log(r/k.radii[i])/k.epsilon[i]
    Rs = k.R[j] + log(s/k.radii[j])/k.epsilon[j]
    result = min(Rr, Rs)*(k.Rtotal-max(Rr, Rs))/k.Rtotal
    if split_images
        result += k.direct[i, j]*log(max(r, s))
        # The nearest images belong to a target and a source in one layer.
        i == j && (result += (k.inner[i]*log(s)+k.outer[i]*(2log(k.radii[i+1])-log(r))) /
                             k.epsilon[i])
    else
        result += k.direct[i, j]*log(max(r, s)/abs(z-source))
        if i == j
            result -= k.inner[i]/k.epsilon[i]*log(abs(1-k.radii[i]^2/(conj(z)*source)))
            result -= k.outer[i]/k.epsilon[i]*log(abs(1-z*conj(source)/k.radii[i+1]^2))
        end
    end
    result
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the arrays that solving `group` uses, sized by `plan`. The `N` sectors
of the rings have `r` targets and `c` sources each. The border consists of the other `mb`
targets and `nb` charge coefficients:

- `ring_kernel` holds the potentials of the rings' sources at the first sector's targets,
  `coupling_kernel` at the border targets, and `border` the potentials of the
  border coefficients at every target. `column_norms` holds the norms that equilibrate
  their columns.
- `excitation` holds the excitation of each terminal at the targets, and `coefficients`
  the equilibrated charge coefficients of the solve.
- `phases` holds the unitary Fourier transform over the ring. `ring_factors` holds the
  ring's Fourier blocks and their QR factors, `ring_border` the border columns in the ring's
  Fourier rows, and `coupling` the ring's Fourier columns in the border rows.
- `ring_residual`, `border_gram`, `gram_product`, `coupling_product`, `border_reduced` and
  `border_system` hold the eliminated border system and its factors.
- `ring_rhs`, `ring_solution`, `border_product`, `border_rhs`, `border_solution`,
  `residual_rhs` and `system_rhs` hold the right-hand sides and solutions of one solve.
- `moments` and `segments` hold the logarithmic moments of one tape face and their QuadGK
  segments.
- `target_inner`, `target_outer`, `source_inner` and `source_outer` hold the ratios
  ``\\rho_{k-1}/\\bar z`` and ``z/\\rho_k`` of the targets and of the sources or face nodes
  to the interfaces of their layer ``k``. The `*_power` arrays hold their running powers.
- `target_modes` and `source_modes` hold the Fourier factors of at most 64 modes at a time.
- `face_kernel` holds the Green function from the targets to the nodes of one face, and
  `face_layers` the layer of those nodes.
- `plus_collocation`, `minus_collocation`, `plus_excitation` and `minus_excitation` hold
  perturbed blocks of at most 192 targets, when `plan.derivatives` directions follow.
- `drive`, `stationarity`, `residual`, `border_normal` and `coefficient_derivative` hold the
  right-hand sides of the capacitance derivative, the residual of the nominal least-squares
  solve and the derivative of the coefficients.
"""
function initialize_buffers(::Formula{:boundary}, ::Type{T}, group::DielectricGroup, plan,
        buffers) where {T}
    (; rows, columns, ports, order, sources, modes, quadrature, derivatives) = plan
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    mb, nb = rows - N*r, columns - N*c
    # Only a ring eliminates the border rows. Without one, they are the whole system.
    kept = iszero(N) ? 0 : mb
    top = min(N*(r - c), nb)
    kernel = max(length(sources), quadrature)
    factors = 4min(64, modes)
    perturbed = iszero(derivatives) ? 0 : rows
    block = min(192, perturbed)
    tangent = iszero(derivatives) ? 0 : columns
    return merge(buffers,
        (ring_kernel = Matrix{T}(undef, r, N*c),
            coupling_kernel = Matrix{T}(undef, kept, N*c),
            border = Matrix{T}(undef, rows, nb),
            column_norms = Vector{T}(undef, columns),
            excitation = Matrix{T}(undef, rows, ports),
            coefficients = Matrix{T}(undef, columns, ports),
            phases = Matrix{Complex{T}}(undef, N, N),
            ring_factors = Array{Complex{T}}(undef, r, c, N),
            ring_border = Array{Complex{T}}(undef, r, nb, N),
            coupling = Matrix{Complex{T}}(undef, kept, N*c),
            ring_residual = Matrix{Complex{T}}(undef, N*(r - c), nb),
            border_gram = Matrix{T}(undef, kept, kept),
            gram_product = Matrix{Complex{T}}(undef, kept, kept),
            coupling_product = Matrix{Complex{T}}(undef, kept, nb),
            border_reduced = Matrix{T}(undef, kept, nb),
            border_system = Matrix{T}(undef, 2top + mb, nb),
            ring_rhs = Array{Complex{T}}(undef, r, ports, N),
            ring_solution = Array{Complex{T}}(undef, c, ports, N),
            border_product = Matrix{Complex{T}}(undef, kept, ports),
            border_rhs = Matrix{T}(undef, mb, ports),
            border_solution = Matrix{T}(undef, nb, ports),
            residual_rhs = Matrix{Complex{T}}(undef, N*(r - c), ports),
            system_rhs = Matrix{T}(undef, 2top + mb, ports),
            moments = Vector{T}(undef, order+1),
            segments = alloc_segbuf(T, Vector{T}, T; size = 32),
            target_inner = Vector{Complex{T}}(undef, rows),
            target_outer = Vector{Complex{T}}(undef, rows),
            source_inner = Vector{Complex{T}}(undef, kernel),
            source_outer = Vector{Complex{T}}(undef, kernel),
            target_inner_power = Vector{Complex{T}}(undef, rows),
            target_outer_power = Vector{Complex{T}}(undef, rows),
            source_inner_power = Vector{Complex{T}}(undef, kernel),
            source_outer_power = Vector{Complex{T}}(undef, kernel),
            target_modes = Matrix{T}(undef, rows, factors),
            source_modes = Matrix{T}(undef, kernel, factors),
            face_kernel = Matrix{T}(undef, rows, quadrature),
            face_layers = Vector{Int}(undef, quadrature),
            plus_collocation = Matrix{T}(undef, block, columns),
            minus_collocation = Matrix{T}(undef, block, columns),
            plus_excitation = Matrix{T}(undef, block, ports),
            minus_excitation = Matrix{T}(undef, block, ports),
            drive = Matrix{T}(undef, perturbed, ports),
            stationarity = Matrix{T}(undef, tangent, ports),
            residual = Matrix{T}(undef, perturbed, ports),
            border_normal = Matrix{Complex{T}}(undef, iszero(tangent) ? 0 : nb, ports),
            coefficient_derivative = Matrix{T}(undef, tangent, ports)))
end

function _shunt_kernel_matrix!(matrix, targets, target_layers, sources, source_layers, k;
        split_images = false, buffers)
    # Same Green function, batched into small Fourier blocks for BLAS. No dense
    # N-by-modes array: the products use only 4*64 columns per geometric boundary set.
    @inbounds for j in eachindex(sources), i in eachindex(targets)

        matrix[i, j] = _shunt_kernel(
            targets[i], target_layers[i], sources[j], source_layers[j], k; split_images)
    end
    nr, nc = length(targets), length(sources)
    ta, tb = @view(buffers.target_inner[1:nr]), @view(buffers.target_outer[1:nr])
    sa, sb = @view(buffers.source_inner[1:nc]), @view(buffers.source_outer[1:nc])
    tpa = @view buffers.target_inner_power[1:nr]
    tpb = @view buffers.target_outer_power[1:nr]
    spa = @view buffers.source_inner_power[1:nc]
    spb = @view buffers.source_outer_power[1:nc]
    # Each point's Fourier factors refer to the interfaces of its own layer.
    @inbounds for i in eachindex(targets)
        ta[i] = k.radii[target_layers[i]]/conj(targets[i])
        tb[i] = targets[i]/k.radii[target_layers[i]+1]
    end
    @inbounds for j in eachindex(sources)
        sa[j] = k.radii[source_layers[j]]/conj(sources[j])
        sb[j] = sources[j]/k.radii[source_layers[j]+1]
    end
    U, V = @view(buffers.target_modes[1:nr, :]), @view(buffers.source_modes[1:nc, :])
    modes = size(k.modal, 1)
    # The targets of one layer that follow each other share their coefficients.
    start = 1
    while start <= nr
        layer = target_layers[start]
        stop = start
        while stop < nr && target_layers[stop + 1] == layer
            stop += 1
        end
        rows = start:stop
        fill!(@view(tpa[rows]), 1)
        fill!(@view(tpb[rows]), 1)
        fill!(spa, 1)
        fill!(spb, 1)
        for first_mode in 1:64:modes
            block = first_mode:min(first_mode + 63, modes)
            for (column, m) in enumerate(block)
                @views tpa[rows] .*= ta[rows]
                @views tpb[rows] .*= tb[rows]
                spa .*= sa
                spb .*= sb
                @inbounds for i in rows
                    U[i - start + 1, 4column - 3]=real(tpa[i])
                    U[i - start + 1, 4column - 2]=imag(tpa[i])
                    U[i - start + 1, 4column - 1]=real(tpb[i])
                    U[i - start + 1, 4column]=imag(tpb[i])
                end
                @inbounds for j in eachindex(sources)
                    source = source_layers[j]
                    a=k.modal[m, 1, 1, layer, source]*spa[j] +
                      k.modal[m, 1, 2, layer, source]*spb[j]
                    b=k.modal[m, 2, 1, layer, source]*spa[j] +
                      k.modal[m, 2, 2, layer, source]*spb[j]
                    V[j, 4column - 3]=real(a)
                    V[j, 4column - 2]=imag(a)
                    V[j, 4column - 1]=real(b)
                    V[j, 4column]=imag(b)
                end
            end
            @views mul!(matrix[rows, :], U[1:length(rows), 1:4length(block)],
                transpose(V[:, 1:4length(block)]), 1.0, 1.0)
        end
        start = stop + 1
    end
    matrix
end

# The potential of the inner anchor at unit voltage, relative to the reference, at a point `z`
# in `layer`.
function _shunt_core_voltage(z, layer, k)
    1-(k.R[layer]+log(abs(z)/k.radii[layer])/k.epsilon[layer])/k.Rtotal
end

# Point `j` of `count` on the circle at `fraction` of a wire's radius, in the wire's own
# polar frame. A rotation that maps one ring wire onto the next maps its points too.
# Round wires keep auxiliary sources inside. The tape has no source contour or inset: its
# physical faces are integrated below.
function _shunt_wire_point(wire, fraction, j, count)
    complex(wire.x, wire.y) + fraction*wire.r*cis(atan(wire.y, wire.x) + 2pi*j/count)
end

"""
$(TYPEDSIGNATURES)

Order the wires of `g` for collocation and find its rings. The wires of a ring share one
terminal, one layer, one radius and one pitch circle, they are equally spaced, and each
exposes the same targets. A ring has at least two wires. A rotation by ``2\\pi/n``, with ``n``
the greatest common divisor of the ring sizes, maps every ring onto itself in the concentric
layers. The rings split into ``n`` sectors that it permutes cyclically, each with ``N_i/n``
consecutive wires of ring ``i``. The sectors come first, then the other wires. A sector needs
at least as many targets as sources, or no wire is reduced.

# Returns

- `wires`: the wire indices in collocation order.
- `ring`: the number of sectors, zero without a ring.
- `sector`: the number of wires in a sector.
- `exposed`: for each wire in that order, the indices of its `2level.wire` surface points
  that lie outside every other wire of its terminal.
"""
function _shunt_layout(g, level)
    count = 2level.wire
    exposed = map(g.wires) do wire
        tolerance = _shunt_tol(g.layers[wire.layer].ro)
        filter(0:(count - 1)) do j
            point = _shunt_wire_point(wire, 1, j, count)
            # Contact points and buried arcs do not expose dielectric to this terminal.
            all(other -> other === wire || other.terminal != wire.terminal ||
                abs(point-complex(other.x, other.y))-other.r > tolerance, g.wires)
        end
    end
    angles = [mod(atan(wire.y, wire.x), 2pi) for wire in g.wires]
    # Wires of one terminal and layer with one radius on one pitch circle, in angular order.
    rings = Vector{Int}[]
    for index in sortperm(angles)
        wire = g.wires[index]
        tolerance = _shunt_tol(g.layers[wire.layer].ro)
        found = findfirst(rings) do candidate
            lead = g.wires[first(candidate)]
            lead.terminal == wire.terminal && lead.layer == wire.layer &&
                abs(lead.r - wire.r) <= tolerance &&
                abs(hypot(lead.x, lead.y) - hypot(wire.x, wire.y)) <= tolerance
        end
        found === nothing ? push!(rings, [index]) : push!(rings[found], index)
    end
    filter!(rings) do ring
        lead = g.wires[first(ring)]
        tolerance = _shunt_tol(g.layers[lead.layer].ro)
        length(ring) > 1 && all(pairs(ring)) do (position, index)
            spacing = angles[index] - angles[first(ring)] - 2pi*(position - 1)/length(ring)
            abs(hypot(lead.x, lead.y)*spacing) <= tolerance &&
                exposed[index] == exposed[first(ring)]
        end
    end
    n = isempty(rings) ? 0 : gcd(length.(rings))
    sector = sum(ring -> length(ring) ÷ n, rings; init = 0)
    reduced = n > 0 &&
              sum(ring -> length(ring) ÷ n * length(exposed[first(ring)]), rings) >=
              sector*level.wire
    wires = Int[]
    if reduced
        for k in 0:(n - 1), ring in rings
            share = length(ring) ÷ n
            append!(wires, view(ring, (k*share + 1):((k + 1)*share)))
        end
    end
    append!(wires, setdiff(eachindex(g.wires), wires))
    return (; wires, ring = reduced ? n : 0, sector = reduced ? sector : 0,
        exposed = exposed[wires])
end

# The coefficients (A, B, C) of the three-term recurrence of the Jacobi polynomials,
# P_{n+1} = (A t + B) P_n - C P_{n-1}, for n = 1, ..., p - 1.
function _shunt_jacobi_recurrence(alpha, beta, p)
    s = alpha+beta
    return map(1:(p - 1)) do n
        ((2n+s+1)*(2n+s+2)/(2(n+1)*(n+s+1)),
            (alpha^2-beta^2)*(2n+s+1)/(2(n+1)*(n+s+1)*(2n+s)),
            (n+alpha)*(n+beta)*(2n+s+2)/((n+1)*(n+s+1)*(2n+s)))
    end
end

function _shunt_jacobi!(values, t, alpha, beta, recurrence)
    values[1] = 1
    length(values) == 1 && return values
    values[2] = ((alpha+beta+2)*t+alpha-beta)/2
    for (n, (A, B, C)) in enumerate(recurrence)
        values[n + 2] = (A*t+B)*values[n + 1]-C*values[n]
    end
    values
end

function _shunt_jacobi(t, alpha, beta, p)
    _shunt_jacobi!(Vector{Float64}(undef, p+1), t, alpha, beta,
        _shunt_jacobi_recurrence(alpha, beta, p))
end

function _shunt_gauss_jacobi(n, alpha, beta)
    s = alpha+beta
    diagonal = [(beta-alpha)/(s+2);
                [(beta^2-alpha^2)/((2k+s)*(2k+s+2)) for k in 1:(n - 1)]]
    off = [2/(s+2)*sqrt((1+alpha)*(1+beta)/(s+3));
           [2/(2k+s)*sqrt(k*(k+alpha)*(k+beta)*(k+s) /
                          ((2k+s-1)*(2k+s+1))) for k in 2:(n - 1)]]
    decomposition = eigen(SymTridiagonal(diagonal, off))
    # Normalized measure w(t)dt / integral(w). No arclength factor belongs here.
    decomposition.values, vec(decomposition.vectors[1, :] .^ 2)
end

"""
$(TYPEDSIGNATURES)

Find the leading electrostatic exponent at a right-angle metal corner on a
dielectric interface. The layer of the metal occupies angle ``\\pi/2`` and the adjacent
dielectric angle ``\\pi``. Continuity of potential and normal displacement,
with constant potential on the two metal faces, gives this specialization:

```math
\\varepsilon_l\\cos(\\nu\\pi/2)\\sin(\\nu\\pi)
+\\varepsilon_o\\sin(\\nu\\pi/2)\\cos(\\nu\\pi)=0.
```

The smallest root ``0<\\nu<1`` gives potential relative to the metal proportional
to ``\\rho^\\nu`` and surface charge density proportional to ``\\rho^{\\nu-1}``,
where ``\\rho`` is distance from the corner \\[m\\]. The admissible branch
obeys Meixner's finite-energy edge condition [Meixner1972](@cite).
Material-dependent metal-dielectric wedge singularities, rather than a
universal homogeneous exponent, are treated by [VanBladel1985](@cite).
For equal permittivities this equation recovers ``\\nu=2/3``.

# Arguments

- `epsilon_layer`: relative permittivity of the metal's layer \\[dimensionless\\].
- `epsilon_outer`: adjacent dielectric relative permittivity \\[dimensionless\\].

# Returns

- Leading exponent ``\\nu`` \\[dimensionless\\], found by bisection excluding
  the trivial zero root. This is the specified two-dielectric corner, not a
  general wedge solver.

# Errors

Nonpositive permittivities raise `BoundarySolveError`. For positive permittivities the
determinant is positive near zero and equals ``-\\varepsilon_o`` at one, so the root is
always bracketed.
"""
function _shunt_junction_exponent(epsilon_layer, epsilon_outer)
    epsilon_layer > 0 && epsilon_outer > 0 || throw(BoundarySolveError(:unsupported,
        (; epsilon_layer, epsilon_outer), "boundary shunt requires positive permittivities"))
    # Layer sector pi/2, outer dielectric sector pi. Determinant has no cotangent
    # poles. Exclude the trivial nu=0 root and bracket the first positive root.
    determinant(nu) = epsilon_layer*cospi(nu/2)*sinpi(nu) +
                      epsilon_outer*sinpi(nu/2)*cospi(nu)
    left, right = 1e-10, 1.0
    for _ in 1:60
        mid = (left+right)/2
        if determinant(mid) > 0
            left = mid
        else
            right = mid
        end
    end
    (left+right)/2
end

# A tape has two circular arc faces and two radial end faces. The record of a face
# holds its geometry, its endpoint exponents, quadrature, weighted Jacobi values, the
# Jacobi recurrence coefficients and the layer of its tape.
const ShuntArc = NamedTuple{(:radius, :phi, :span), NTuple{3, Float64}}
const ShuntEnd = NamedTuple{(:phi, :mid, :half), NTuple{3, Float64}}
const ShuntFace{S} = NamedTuple{
    (:shape, :alpha, :beta, :p, :nodes, :weights, :weighted, :beta_norm, :terminal, :points,
        :recurrence, :layer),
    Tuple{S, Float64, Float64, Int, Vector{Float64}, Vector{Float64}, Matrix{Float64},
        Float64, Int, Vector{ComplexF64}, Vector{NTuple{3, Float64}}, Int}}
# One tape has an inner arc face, an outer arc face and two end faces, in that order.
const ShuntTapeFaces = Tuple{ShuntFace{ShuntArc}, ShuntFace{ShuntArc},
    ShuntFace{ShuntEnd}, ShuntFace{ShuntEnd}}

_shunt_face_point(face::ShuntArc, t) = face.radius*cis(face.phi+face.span*t/2)
_shunt_face_point(face::ShuntEnd, t) = (face.mid+face.half*t)*cis(face.phi)

"""
$(TYPEDSIGNATURES)

Build whole-face charge expansions on the tape's two circular arcs and two
radial ends, retaining its physical thickness. The finite-face precedent is
[Bernal1997](@cite), Section IV, Eq. (10), which uses Maxwell-weighted
Chebyshev functions. This implementation instead uses material-dependent
Jacobi weights, with exponents from `_shunt_junction_exponent` at dielectric
interfaces and ``\\nu=2/3`` at homogeneous right-angle corners:

```math
d\\widehat q(t)=\\sum_{n=0}^{p} a_n P_n^{(\\alpha,\\beta)}(t)\\,d\\mu(t),
\\qquad
d\\mu(t)=\\frac{(1-t)^\\alpha(1+t)^\\beta\\,dt}
{2^{\\alpha+\\beta+1}\\mathrm{B}(\\alpha+1,\\beta+1)},
\\qquad \\alpha=\\nu_+-1,\\quad\\beta=\\nu_--1.
```

Here ``t\\in[-1,1]`` is the face parameter, ``\\nu_\\pm`` are its endpoint
exponents, and ``P_n`` are Jacobi polynomials. The normalized charge measure
``d\\widehat q=dq/(2\\pi\\varepsilon_0)`` and coefficients ``a_n`` have units
\\[V\\]. ``dq`` is charge per cable length \\[C/m\\]. Each face's total charge
is ``2\\pi\\varepsilon_0 a_0``. Do not multiply this measure by a second
arclength Jacobian. Encoding known singularities in the approximation space
is also supported by [Classen2011](@cite), Sections 2 and 4. Their FIT/DG
discretizations are not the boundary-integral element implemented here.

# Arguments

- `g`: the dielectric layers' radii \\[m\\] and relative permittivities
  \\[dimensionless\\].
- `s`: tape geometry, terminal index and layer. The inner and outer radii are in
  \\[m\\], and the center angle and angular span are in \\[rad\\].
- `p`: highest Jacobi degree on each face.
- `quadrature`: number of normalized Gauss–Jacobi nodes per face.

# Returns

- The 4 face records of the tape, with physical coordinates \\[m\\], endpoint exponents,
  quadrature weights and weighted polynomial values.
"""
function _shunt_tape_faces(g, s, p, quadrature)
    # A metal face on a dielectric interface meets that dielectric at its corners.
    layer = g.layers[s.layer]
    at_inner = s.layer > 1 && abs(s.ri - layer.ri) <= _shunt_tol(layer.ro)
    at_outer = s.layer < length(g.layers) && abs(s.ro - layer.ro) <= _shunt_tol(layer.ro)
    nu_inner = at_inner ? _shunt_junction_exponent(layer.epsilon, g.layers[s.layer-1].epsilon) :
               2/3
    nu_outer = at_outer ? _shunt_junction_exponent(layer.epsilon, g.layers[s.layer+1].epsilon) :
               2/3
    arc(radius, nu) = (; shape = ShuntArc((radius, s.phi, s.span)), alpha = nu-1, beta = nu-1)
    radial(sign) = (; shape = ShuntEnd((s.phi+sign*s.span/2, (s.ri+s.ro)/2, (s.ro-s.ri)/2)),
        alpha = nu_outer-1, beta = nu_inner-1)
    descriptions = (arc(s.ri, nu_inner), arc(s.ro, nu_outer), radial(-1), radial(1))
    map(descriptions) do face
        nodes, weights = _shunt_gauss_jacobi(quadrature, face.alpha, face.beta)
        polys = transpose(reduce(hcat, [_shunt_jacobi(t, face.alpha, face.beta, p)
                                        for t in nodes]))
        weighted = weights .* polys
        loggamma = SpecialFunctions.loggamma
        beta_norm = exp(loggamma(face.alpha+1)+loggamma(face.beta+1) -
                        loggamma(face.alpha+face.beta+2))
        merge(face,
            (; p, nodes, weights, weighted, beta_norm, terminal = s.terminal,
                points = _shunt_face_point.(Ref(face.shape), nodes),
                recurrence = _shunt_jacobi_recurrence(face.alpha, face.beta, p),
                layer = s.layer))
    end
end

function _shunt_face_projection(z, face::ShuntArc)
    radius = abs(z)
    angle_difference = atan(sin(angle(z)-face.phi), cos(angle(z)-face.phi))
    u = 2angle_difference/face.span
    delta = abs(radius-face.radius)/(face.radius*face.span/2)
    u, delta
end

function _shunt_face_projection(z, face::ShuntEnd)
    rotated = z*cis(-face.phi)
    u = (real(rotated)-face.mid)/face.half
    delta = abs(imag(rotated))/face.half
    u, delta
end

# Distance from `z` to the face point at the cosine-parameter offset `difference`.
function _shunt_face_distance(face::ShuntArc, z, difference)
    # Exact opposed-arc distance, also valid in the self limit.
    hypot(abs(z)-face.radius, 2sqrt(abs(z)*face.radius)*sin(face.span*difference/4))
end
_shunt_face_distance(face::ShuntEnd, z, difference) =
    hypot(face.half*difference, imag(z*cis(-face.phi)))

"""
$(TYPEDSIGNATURES)

Integrate the logarithmic kernel against the normalized face charge modes:

```math
M_n(z)=\\int_{-1}^{1} P_n^{(\\alpha,\\beta)}(t)
\\log\\!\\left(\\frac{|z-\\zeta(t)|}{1\\,\\mathrm{m}}\\right)d\\mu(t).
```

The face position is ``\\zeta(t)`` \\[m\\]. ``d\\mu`` is the normalized measure
in `_shunt_tape_faces`. Ordinary quadrature can lose accuracy close to a
source boundary [HelsingOjala2008](@cite), Section 1. Here direct and image
logarithms are integrated separately from the smooth Green remainder.
Nearby targets use ``t=\\cos\\theta`` and adaptive Gauss–Kronrod integration,
split at the projected singularity or near peak. Distant targets use the
face's Gauss–Jacobi rule. This is not Helsing-Ojala's rational-quadrature
algorithm, and that reference does not prescribe the corner basis.

# Arguments

- `z`: target coordinate in the complex cross-sectional plane \\[m\\].
- `face`: physical face geometry, normalized weights and Jacobi degree.

# Keywords

- `rtol`, `atol`, `maxevals`: QuadGK's relative and absolute tolerances for the
  dimensionless moment vector, and its evaluation budget.
- `result`: in-place moment vector, also returned.
- `segments`: reusable QuadGK segment buffer.

# Returns

- Logarithmic moments \\[dimensionless\\].
- QuadGK's error estimate. The far-field fixed rule returns zero because it
  supplies no estimate, not because its error is known to vanish.

# Errors

An unresolved zero-distance evaluation or a nonfinite moment raises
`BoundarySolveError`. The logarithm is not replaced by an arbitrary finite floor.
"""
function _shunt_log_moments(z, face; rtol, atol, maxevals, result, segments)
    u, delta = _shunt_face_projection(z, face.shape)
    if hypot(max(abs(u)-1, 0), delta) > 0.2
        fill!(result, 0)
        @inbounds for j in eachindex(face.points)
            value = log(abs(face.points[j]-z))
            for n in eachindex(result)
                result[n] += face.weighted[j, n]*value
            end
        end
        return result, 0.0
    end
    # Weighted logarithmic product integration: compute the log moments of the
    # Jacobi modes, independently of the smooth Green remainder. t=cos(theta)
    # absorbs endpoint weights. Split at the logarithmic singularity or near peak.
    theta0 = acos(clamp(u, -1, 1))
    splits = [0.0, theta0, pi]
    if delta > 0
        width = min(sqrt(delta), delta/max(sin(theta0), eps(Float64)))
        for sign in (-1, 1), factor in (1, 4)

            push!(splits, clamp(theta0+sign*factor*width, 0, pi))
        end
    end
    sort!(unique!(splits))
    function integrand!(output, theta)
        t = cos(theta)
        difference = abs(u) <= 1 ?
                     -2sin((theta+theta0)/2)*sin((theta-theta0)/2) : t-u
        distance = _shunt_face_distance(face.shape, z, difference)
        # Roundoff may identify the endpoint only after its contribution is
        # below integration accuracy. Never manufacture a finite log(0) floor.
        distance > 0 || throw(BoundarySolveError(:quadrature,
            (; point = z, theta), "unresolved logarithmic integration point"))
        weight = sin(theta/2)^(2face.alpha+1)*cos(theta/2)^(2face.beta+1)/face.beta_norm
        _shunt_jacobi!(output, t, face.alpha, face.beta, face.recurrence)
        output .*= log(distance)*weight
    end
    # Reuse quadrature work vectors. Modal integration must not allocate a new
    # polynomial vector at each of its thousands of function evaluations.
    integral,
    estimate = quadgk!(integrand!, result, splits;
        segbuf = segments, rtol, atol, order = max(7, cld(face.p+1, 2)), maxevals)
    all(isfinite, integral) || throw(BoundarySolveError(:nonfinite,
        (; point = z), "nonfinite logarithmic integral"))
    return integral, estimate
end

# Add the integrated charge columns of every tape face to `columns`, for targets in
# `target_layers`. A target in the face's layer integrates the direct logarithm and the two
# nearest images, and a target in another layer the transmitted logarithm. On an interface, an
# image is the target itself. Count the moments whose QuadGK estimate exceeds the requested
# target.
function _shunt_tape_columns!(columns, targets, target_layers, faces, k; integration, buffers)
    estimated_error = 0.0
    unmet = 0
    offset = 0
    for tape in faces, face in tape
        block = @view columns[:, (offset + 1):(offset + face.p + 1)]
        kernel = @view buffers.face_kernel[1:length(targets), 1:length(face.points)]
        layers = fill!(@view(buffers.face_layers[1:length(face.points)]), face.layer)
        _shunt_kernel_matrix!(
            kernel, targets, target_layers, face.points, layers, k; split_images = true, buffers)
        mul!(block, kernel, face.weighted)
        for (i, z) in pairs(targets)
            layer = target_layers[i]
            same = layer == face.layer
            inner, outer = k.radii[layer]^2/conj(z), k.radii[layer + 1]^2/conj(z)
            inner_same = abs(inner-z) <= _shunt_tol(k.radii[layer + 1])
            outer_same = abs(outer-z) <= _shunt_tol(k.radii[layer + 1])
            inner_image = same ? k.inner[layer]/k.epsilon[layer] : 0.0
            outer_image = same ? k.outer[layer]/k.epsilon[layer] : 0.0
            direct = k.direct[layer, face.layer] + (inner_same ? inner_image : 0.0) +
                     (outer_same ? outer_image : 0.0)
            for (point, coefficient) in ((inner, inner_same ? 0.0 : inner_image),
                (outer, outer_same ? 0.0 : outer_image), (z, direct))
                coefficient == 0 && continue
                moments,
                err = _shunt_log_moments(point, face; integration...,
                    result = buffers.moments, segments = buffers.segments)
                @inbounds for n in eachindex(moments)
                    block[i, n] -= coefficient*moments[n]
                end
                estimated_error = max(estimated_error, abs(coefficient)*err)
                unmet += !(err <= max(integration.rtol*norm(moments), integration.atol))
            end
        end
        offset += face.p+1
    end
    (; columns, log_moment_error = estimated_error, unmet)
end

"""
Record one boundary solve. `log_moment_error` is the largest weighted QuadGK error
estimate of the logarithmic moments. `condition` is the 1-norm condition estimate of the
equilibrated collocation system's triangular factors, the larger of the Fourier blocks
together and the border. `ring` is the number of Fourier classes that the ring reduction
separates, zero without a ring. `unknowns` and `equations` count the charge coefficients and
collocation equations, and `level` is the resolution.
"""
const InternalShuntDiagnostic = NamedTuple{
    (:log_moment_error, :condition, :ring, :unknowns, :equations, :level),
    Tuple{Float64, Float64, Int, Int, Int, typeof(DEFAULT_RESOLUTION)}}

# The fixed geometry of one solve in the wire order of `layout`: tape faces, collocation
# targets with their terminals and layers, auxiliary sources with their layers and the
# Fourier coefficients, with the array sizes they give. The sectors' targets and sources come
# first, wire by wire, then the other wires' and the tape faces'.
function _shunt_plan(selected::Formula{:boundary}, g, layout, derivatives = 0)
    (; resolution, integration) = selected.options.data
    level = resolution
    faces = ShuntTapeFaces[_shunt_tape_faces(g, tape, level.order, level.quadrature)
                           for tape in g.tapes]
    points, terminals, point_layers = ComplexF64[], Int[], Int[]
    sources, source_layers = ComplexF64[], Int[]
    for (index, exposed) in zip(layout.wires, layout.exposed)
        wire = g.wires[index]
        append!(points, _shunt_wire_point(wire, 1, j, 2level.wire) for j in exposed)
        append!(terminals, fill(wire.terminal, length(exposed)))
        append!(point_layers, fill(wire.layer, length(exposed)))
        append!(sources,
            _shunt_wire_point(wire, 0.75, j, level.wire) for j in 0:(level.wire - 1))
        append!(source_layers, fill(wire.layer, level.wire))
    end
    nface = max(24, 4(level.order+1))
    parameters = -cos.(pi .* ((1:nface) .- 0.5) ./ nface)
    for tape in faces, face in tape
        append!(points, [_shunt_face_point(face.shape, t) for t in parameters])
        append!(terminals, fill(face.terminal, nface))
        append!(point_layers, fill(face.layer, nface))
    end
    k = _shunt_kernel_coefficients(g, level.modes)
    columns = length(sources)+4(level.order+1)*length(faces)
    ports = max(maximum(wire -> wire.terminal, g.wires; init = 1),
        maximum(tape -> tape.terminal, g.tapes; init = 1))
    ring = layout.ring
    return (; faces, points, terminals, point_layers, sources, source_layers, k, level,
        integration, layout, ring,
        block_rows = sum(length, view(layout.exposed, 1:layout.sector); init = 0),
        block_columns = layout.sector*level.wire,
        rows = length(points), columns, ports, modes = level.modes,
        quadrature = isempty(faces) ? 0 : level.quadrature, order = level.order, derivatives)
end

function _shunt_rhs!(rhs, points, point_layers, terminals, k)
    fill!(rhs, 0)
    @inbounds for i in eachindex(points)
        rhs[i, 1] = -_shunt_core_voltage(points[i], point_layers[i], k)
        rhs[i, terminals[i]] = 1
    end
    return rhs
end

function _shunt_matrix!(matrix, points, point_layers, sources, source_layers, faces, k;
        integration, buffers)
    count = length(sources)
    _shunt_kernel_matrix!(
        @view(matrix[:, 1:count]), points, point_layers, sources, source_layers, k; buffers)
    return _shunt_tape_columns!(@view(matrix[:, (count + 1):end]), points, point_layers,
        faces, k; integration, buffers)
end

# Total charge of each terminal and the charge moment of the core potential, per face.
function _shunt_charge!(charge, offset, face::ShuntFace{ShuntArc}, k)
    charge[face.terminal, offset + 1] = 1
    # The core potential is constant on an arc: only the uniform mode has a moment.
    charge[1, offset + 1] = -_shunt_core_voltage(first(face.points), face.layer, k)
    return charge
end

function _shunt_charge!(charge, offset, face::ShuntFace{ShuntEnd}, k)
    charge[face.terminal, offset + 1] = 1
    @views mul!(charge[1:1, (offset + 1):(offset + face.p + 1)],
        transpose(_shunt_core_voltage.(face.points, face.layer, Ref(k))), face.weighted,
        -1.0, 0.0)
    return charge
end

function _shunt_charge_map(g, plan)
    (; sources, source_layers, faces, level, k) = plan
    charge = zeros(plan.ports, plan.columns)
    @inbounds for (position, index) in pairs(plan.layout.wires)
        for j in ((position - 1) * level.wire + 1):(position * level.wire)
            charge[1, j] = -_shunt_core_voltage(sources[j], source_layers[j], k)
            charge[g.wires[index].terminal, j] = 1
        end
    end
    offset = length(sources)
    for tape in faces, face in tape
        _shunt_charge!(charge, offset, face, k)
        offset += face.p+1
    end
    return charge
end

"""
$(TYPEDSIGNATURES)

Compute the lossless terminal capacitance of exposed wires and finite tapes
inside a layered circular reference shield. Round-wire auxiliary sources and
corner-weighted integrated tape charges share a column-equilibrated least-
squares solve. A ring of identical, equally spaced wires on one terminal separates into one
block per Fourier index, and the tapes and other wires form a border that is eliminated
exactly. Charge is extracted by exact total-charge and core-potential
moments, not by symmetrizing the answer. The layout and plan fix the geometry of the solve,
the buffers of the boundary shunt formula hold its arrays, and `directions` adds the arrays
of the capacitance derivative.

# Arguments

- `selected`: the boundary shunt formula. Its options set the resolution and the
  logarithmic-moment quadrature controls. The reference settings are not an
  accuracy guarantee for arbitrary geometry or weak terminal couplings.
- `g`: local radii and positions \\[m\\] and relative dielectric permittivities.
- `directions=nothing`: physical input perturbations, when derivatives follow.

# Returns

- Local capacitance `C` \\[F/m\\], the numerical record, the factors of the solve, the
  charge coefficients and charge map, and the workspace of the solve.

# Notes

The layered Green function and whole-face charge approach has precedent in
[Bernal1997](@cite), but this solve uses oversampled boundary collocation and
unpivoted QR, not that paper's Galerkin system. The circular kernel, wire-source
placement at `0.75` times each wire radius, Jacobi adaptation and resolution are
implementation choices. An accuracy bound for this combined discretization is not
established by the cited methods.

# Errors

Fewer equations than charge coefficients, degenerate charge columns and nonfinite
capacitance are `BoundarySolveError`, not geometry `DomainError`. Monte Carlo must not
condition its samples on them. An unmet quadrature target is reported as a warning, once
per solve. The condition estimate is recorded, without a judgement.
"""
function _shunt_capacitance(selected::Formula{:boundary}, group, directions = nothing)
    layout = _shunt_layout(group, selected.options.data.resolution)
    plan = _shunt_plan(selected, group, layout,
        directions === nothing ? 0 : size(directions, 2))
    buffers = initialize_buffers(selected, Float64, group, plan, (;))
    workspace = (; input = group, plan, buffers)
    return merge(_shunt_capacitance!(selected, workspace), (; workspace))
end

# Assemble the collocation system in the buffers of `workspace` and solve it.
function _shunt_capacitance!(selected::Formula{:boundary}, workspace)
    g, plan, buffers = workspace.input, workspace.plan, workspace.buffers
    (; faces, points, terminals, point_layers, sources, source_layers, k, integration) = plan
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    m, n = plan.rows, plan.columns
    m >= n || throw(BoundarySolveError(:rank, (; rows = m, columns = n),
        "fewer boundary equations than charge coefficients"))
    (; ring_kernel, coupling_kernel, border, column_norms) = buffers
    kept = size(coupling_kernel, 1)
    # The rings need one block row of their kernel, and their columns at the border rows.
    ring, ring_layers = @view(sources[1:(N * c)]), @view(source_layers[1:(N * c)])
    _shunt_kernel_matrix!(ring_kernel, @view(points[1:r]), @view(point_layers[1:r]), ring,
        ring_layers, k; buffers)
    border_rows = (N * r + 1):(N * r + kept)
    _shunt_kernel_matrix!(coupling_kernel, @view(points[border_rows]),
        @view(point_layers[border_rows]), ring, ring_layers, k; buffers)
    tape = _shunt_matrix!(border, points, point_layers, @view(sources[(N * c + 1):end]),
        @view(source_layers[(N * c + 1):end]), faces, k; integration, buffers)
    tape.unmet == 0 ||
        @warn "Boundary logarithmic quadrature target was not met" integrals=tape.unmet log_moment_error=tape.log_moment_error rtol=integration.rtol atol=integration.atol maxevals=integration.maxevals
    # A ring column's norm over the ring rows is the same in every sector, which keeps the
    # equilibrated blocks circulant.
    for j in 1:c
        scale = norm(@view ring_kernel[:, j:c:end])
        isfinite(scale) && scale > 0 || throw(BoundarySolveError(
            :rank, (; column = j), "degenerate charge column"))
        column_norms[j:c:(N * c)] .= scale
        @views ring_kernel[:, j:c:end] ./= scale
        @views coupling_kernel[:, j:c:end] ./= scale
    end
    for j in axes(border, 2)
        scale = norm(@view border[:, j])
        isfinite(scale) && scale > 0 || throw(BoundarySolveError(
            :rank, (; column = N*c+j), "degenerate charge column"))
        column_norms[N * c + j] = scale
        @views border[:, j] ./= scale
    end
    factored = _shunt_factor!(plan, buffers)
    rhs = _shunt_rhs!(buffers.excitation, points, point_layers, terminals, k)
    coefficients = _shunt_solve!(buffers.coefficients, factored, plan, buffers, rhs) ./
                   column_norms
    charge = _shunt_charge_map(g, plan)
    C = (2pi*vacuum_permittivity(Float64)) .* (charge*coefficients)
    C[1, 1] += 2pi*vacuum_permittivity(Float64)/k.Rtotal
    all(isfinite, C) ||
        throw(BoundarySolveError(:nonfinite, (;), "nonfinite terminal capacitance"))
    diagnostic = InternalShuntDiagnostic((
        tape.log_moment_error, factored.condition, N, n, m, plan.level))
    return (; C, diagnostic, factored, coefficients, charge)
end

"""
$(TYPEDSIGNATURES)

Factor the equilibrated collocation system that `buffers` holds for `plan`. A unitary
discrete Fourier transform over the ring separates its circulant part into one block
``\\Lambda_p=\\sum_l B_l e^{2\\pi i p l/N}`` per Fourier index ``p``, where ``B_l`` is the
potential of wire ``l``'s sources at the first wire's targets. Unpivoted QR factors each
block, ``\\Lambda_p=Q_p R_p``. The border rows enter through ``M=G R^{-1}`` and the Cholesky
factor ``L`` of ``K=I+\\Re(M M^H)``, where ``G`` holds the ring's Fourier columns at the border
rows. The border columns enter through one QR factor of the blocks' residual rows stacked on
``L^{-1}(H-\\Re(M S))``, where ``H`` and ``S`` hold the border columns at the border rows and
at the blocks' leading rows. Without a ring, that factor is the QR factor of the whole system.

# Returns

- The factors that `_shunt_solve!` and `_shunt_normal!` use, and `condition`, the larger of
  the 1-norm condition estimates of the Fourier blocks together and of the border factor.
"""
function _shunt_factor!(plan, buffers)
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    nb = plan.columns - N*c
    (; ring_kernel, coupling_kernel, border, phases, ring_factors, ring_border, coupling,
        ring_residual, border_gram, gram_product, coupling_product, border_reduced,
        border_system) = buffers
    kept = size(coupling, 1)
    for l in 0:(N - 1), p in 0:(N - 1)
        phases[p + 1, l + 1] = cispi(2mod(p*l, N)/N)/sqrt(N)
    end
    fill!(ring_factors, 0)
    fill!(ring_border, 0)
    fill!(coupling, 0)
    for l in 0:(N - 1), p in 0:(N - 1)
        phase = phases[p + 1, l + 1]
        block, rows = (l * c + 1):((l + 1) * c), (l * r + 1):((l + 1) * r)
        @views ring_factors[:, :, p + 1] .+= (sqrt(N)*phase) .* ring_kernel[:, block]
        @views ring_border[:, :, p + 1] .+= conj(phase) .* border[rows, :]
        @views coupling[:, (p * c + 1):((p + 1) * c)] .+= phase .* coupling_kernel[:, block]
    end
    factors = [qr!(@view ring_factors[:, :, p]) for p in 1:N]
    largest, inverse = 0.0, 0.0
    for p in 1:N
        lmul!(factors[p].Q', @view ring_border[:, :, p])
        R = UpperTriangular(@view ring_factors[1:c, :, p])
        magnitude = opnorm(R, 1)
        largest = max(largest, magnitude)
        inverse = max(inverse, cond(R, 1)/magnitude)
        rdiv!(@view(coupling[:, ((p - 1) * c + 1):(p * c)]), R)
    end
    # Without a ring or without border rows, K is empty.
    mul!(gram_product, coupling, coupling')
    border_gram .= real.(gram_product)
    for i in 1:kept
        border_gram[i, i] += 1
    end
    gram = cholesky!(Hermitian(border_gram))
    fill!(coupling_product, 0)
    for p in 1:N
        @views mul!(coupling_product, coupling[:, ((p - 1) * c + 1):(p * c)],
            ring_border[1:c, :, p], true, true)
    end
    @views border_reduced .= border[(N * r + 1):(N * r + kept), :] .- real.(coupling_product)
    ldiv!(gram.L, border_reduced)
    for p in 1:N
        @views ring_residual[((p - 1) * (r - c) + 1):(p * (r - c)), :] .=
            ring_border[(c + 1):r, :, p]
    end
    residual = qr!(ring_residual)
    top = min(N*(r - c), nb)
    R = residual.R
    @views border_system[1:top, :] .= real.(R)
    @views border_system[(top + 1):(2top), :] .= imag.(R)
    @views border_system[(2top + 1):end, :] .= iszero(N) ? border : border_reduced
    system = qr!(border_system)
    condition = max(largest*inverse, iszero(nb) ? 1.0 :
                                     cond(UpperTriangular(@view system.factors[1:nb, 1:nb]), 1))
    return (; factors, gram, residual, top, system, condition)
end

"""
$(TYPEDSIGNATURES)

Write the equilibrated coefficients of the least-squares solution for the right-hand sides
`rhs` into `x`. With ``d_p`` the leading rows of ``Q_p^H`` times the transformed ring rows,
``w=L^{-1}(b_b-\\Re\\sum_p M_p d_p)`` and the blocks' residual rows give the border
coefficients ``t`` by the border factor. Then ``v=L^{-H}(w-L^{-1}\\tilde H t)``, each block's
coefficients are ``y_p=R_p^{-1}(d_p+M_p^H v-S_p t)``, and the inverse transform returns the
ring's coefficients.
"""
function _shunt_solve!(x, factored, plan, buffers, rhs)
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    nb = plan.columns - N*c
    (; factors, gram, residual, top, system) = factored
    (; phases, ring_factors, ring_border, coupling, border_reduced, ring_rhs, ring_solution,
        border_product, border_rhs, border_solution, residual_rhs, system_rhs) = buffers
    kept = size(coupling, 1)
    fill!(ring_rhs, 0)
    for q in 0:(N - 1), p in 0:(N - 1)
        @views ring_rhs[:, :, p + 1] .+= conj(phases[p + 1, q + 1]) .*
                                         rhs[(q * r + 1):((q + 1) * r), :]
    end
    for p in 1:N
        lmul!(factors[p].Q', @view ring_rhs[:, :, p])
    end
    w = border_rhs
    w .= @view rhs[(N * r + 1):end, :]
    fill!(border_product, 0)
    for p in 1:N
        @views mul!(border_product, coupling[:, ((p - 1) * c + 1):(p * c)],
            ring_rhs[1:c, :, p], true, true)
    end
    @views w[1:kept, :] .-= real.(border_product)
    @views ldiv!(gram.L, w[1:kept, :])
    for p in 1:N
        @views residual_rhs[((p - 1) * (r - c) + 1):(p * (r - c)), :] .=
            ring_rhs[(c + 1):r, :, p]
    end
    lmul!(residual.Q', residual_rhs)
    @views system_rhs[1:top, :] .= real.(residual_rhs[1:top, :])
    @views system_rhs[(top + 1):(2top), :] .= imag.(residual_rhs[1:top, :])
    @views system_rhs[(2top + 1):end, :] .= w
    lmul!(system.Q', system_rhs)
    t = border_solution
    t .= @view system_rhs[1:nb, :]
    ldiv!(UpperTriangular(@view system.factors[1:nb, 1:nb]), t)
    @views mul!(w[1:kept, :], border_reduced, t, -1, 1)
    @views ldiv!(gram.U, w[1:kept, :])
    for p in 1:N
        y = @view ring_solution[:, :, p]
        @views y .= ring_rhs[1:c, :, p]
        @views mul!(y, coupling[:, ((p - 1) * c + 1):(p * c)]', w[1:kept, :], true, true)
        @views mul!(y, ring_border[1:c, :, p], t, -1, 1)
        ldiv!(UpperTriangular(@view ring_factors[1:c, :, p]), y)
    end
    return _shunt_coefficients!(x, plan, buffers)
end

"""
$(TYPEDSIGNATURES)

Write ``(A^T A)^{-1} w`` of the factored, equilibrated system ``A`` into `x`, which may be
`w`, for the implicit derivative. The blocks give ``c_p=R_p^{-H}\\hat w_p`` from the
transformed ring part ``\\hat w`` of `w`. The border coefficients ``q_t`` solve
``R_T^T R_T q_t=w_t-\\Re\\sum_p S_p^H c_p-\\tilde W^T g`` with the border factor ``R_T``,
``g=L^{-1}\\Re\\sum_p M_p c_p`` and ``\\tilde W=L^{-1}\\tilde H``. Each block's coefficients
are then ``R_p^{-1}(c_p-M_p^H L^{-H}(g+\\tilde W q_t)-S_p q_t)``, transformed back.
"""
function _shunt_normal!(x, factored, plan, buffers, w)
    N, c = plan.ring, plan.block_columns
    nb = plan.columns - N*c
    (; gram, system) = factored
    (; phases, ring_factors, ring_border, coupling, border_reduced, ring_solution,
        border_product, border_rhs, border_solution, border_normal) = buffers
    kept = size(coupling, 1)
    fill!(ring_solution, 0)
    for l in 0:(N - 1), p in 0:(N - 1)
        @views ring_solution[:, :, p + 1] .+= conj(phases[p + 1, l + 1]) .*
                                              w[(l * c + 1):((l + 1) * c), :]
    end
    for p in 1:N
        ldiv!(UpperTriangular(@view ring_factors[1:c, :, p])', @view ring_solution[:, :, p])
    end
    t = border_solution
    t .= @view w[(N * c + 1):end, :]
    fill!(border_normal, 0)
    fill!(border_product, 0)
    for p in 1:N
        @views mul!(border_normal, ring_border[1:c, :, p]', ring_solution[:, :, p], true, true)
        @views mul!(border_product, coupling[:, ((p - 1) * c + 1):(p * c)],
            ring_solution[:, :, p], true, true)
    end
    t .-= real.(border_normal)
    g = @view border_rhs[1:kept, :]
    g .= real.(border_product)
    ldiv!(gram.L, g)
    mul!(t, border_reduced', g, -1, 1)
    R = UpperTriangular(@view system.factors[1:nb, 1:nb])
    ldiv!(R', t)
    ldiv!(R, t)
    mul!(g, border_reduced, t, 1, 1)
    ldiv!(gram.U, g)
    for p in 1:N
        y = @view ring_solution[:, :, p]
        @views mul!(y, coupling[:, ((p - 1) * c + 1):(p * c)]', g, -1, 1)
        @views mul!(y, ring_border[1:c, :, p], t, -1, 1)
        ldiv!(UpperTriangular(@view ring_factors[1:c, :, p]), y)
    end
    return _shunt_coefficients!(x, plan, buffers)
end

# The equilibrated coefficients from the ring's Fourier coefficients and the border's: wire
# l's are the real part of the inverse transform, sum over p of exp(2pi i p l/N)/sqrt(N) y_p.
function _shunt_coefficients!(x, plan, buffers)
    N, c = plan.ring, plan.block_columns
    (; phases, ring_solution, border_solution) = buffers
    for l in 0:(N - 1)
        block = @view x[(l * c + 1):((l + 1) * c), :]
        fill!(block, 0)
        for p in 0:(N - 1)
            @views block .+= real.(phases[p + 1, l + 1] .* ring_solution[:, :, p + 1])
        end
    end
    @views x[(N * c + 1):end, :] .= border_solution
    return x
end

# Write the residual b - A y of the equilibrated system for the coefficients `y` and the
# right-hand sides `rhs` into `residual`. Block (q, l) of the ring's circulant is block
# (0, l - q mod N) of its first block row.
function _shunt_residual!(residual, plan, buffers, y, rhs)
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    (; ring_kernel, coupling_kernel, border) = buffers
    kept = size(coupling_kernel, 1)
    copyto!(residual, rhs)
    for q in 0:(N - 1), l in 0:(N - 1)
        d = mod(l - q, N)
        @views mul!(residual[(q * r + 1):((q + 1) * r), :],
            ring_kernel[:, (d * c + 1):((d + 1) * c)], y[(l * c + 1):((l + 1) * c), :], -1, 1)
    end
    @views mul!(residual[(N * r + 1):(N * r + kept), :], coupling_kernel, y[1:(N * c), :],
        -1, 1)
    @views mul!(residual, border, y[(N * c + 1):end, :], -1, 1)
    return residual
end

"""
$(TYPEDSIGNATURES)

Differentiate the fixed-resolution local least-squares problem in one physical
parameter direction. Centered differences evaluate only kernel and charge-map
derivatives. The nominal factors are reused. Include the nonzero collocation
residual in the differentiated stationarity equation:

```math
A^T A\\,dx=A^T(db-dA\\,x)+dA^T(b-Ax).
```

The factored solves evaluate this equation without forming normal equations. The perturbed
groups keep the nominal layout. Their points correspond to the nominal points.
Geometric boundary derivatives are streamed in bounded row blocks, not stored as a dense
matrix of uncertain scalars. Input directions retain their original physical
units. The direction parameter and `step` are dimensionless.

# Returns

- Directional derivative of terminal capacitance \\[F/m\\].
"""
function _shunt_tangent(selected::Formula{:boundary}, group, values, direction, solved, step)
    (; plan, buffers) = solved.workspace
    plus = DielectricGroup(group, values .+ step .* direction)
    minus = DielectricGroup(group, values .- step .* direction)
    dp = _shunt_plan(selected, plus, plan.layout)
    dm = _shunt_plan(selected, minus, plan.layout)
    m = plan.rows
    count = size(buffers.plus_collocation, 1)
    stationarity = fill!(buffers.stationarity, 0)
    for start in 1:count:m
        stop = min(start+count-1, m)
        rows = 1:(stop - start + 1)
        ap = @view buffers.plus_collocation[rows, :]
        am = @view buffers.minus_collocation[rows, :]
        xp, xm = @view(dp.points[start:stop]), @view(dm.points[start:stop])
        layers = @view plan.point_layers[start:stop]
        _shunt_matrix!(ap, xp, layers, dp.sources, plan.source_layers, dp.faces, dp.k;
            plan.integration, buffers)
        _shunt_matrix!(am, xm, layers, dm.sources, plan.source_layers, dm.faces, dm.k;
            plan.integration, buffers)
        ap .-= am
        ap ./= 2step
        rp = @view buffers.plus_excitation[rows, :]
        rm = @view buffers.minus_excitation[rows, :]
        _shunt_rhs!(rp, xp, layers, @view(dp.terminals[start:stop]), dp.k)
        _shunt_rhs!(rm, xm, layers, @view(dm.terminals[start:stop]), dm.k)
        target = @view buffers.drive[start:stop, :]
        target .= (rp .- rm) ./ (2step)
        mul!(target, ap, solved.coefficients, -1.0, 1.0)
        mul!(stationarity, transpose(ap), @view(buffers.residual[start:stop, :]), 1.0, 1.0)
    end
    stationarity ./= buffers.column_norms
    delta = _shunt_solve!(
        buffers.coefficient_derivative, solved.factored, plan, buffers, buffers.drive)
    delta .+= _shunt_normal!(stationarity, solved.factored, plan, buffers, stationarity)
    delta ./= buffers.column_norms
    charge_delta = (_shunt_charge_map(plus, dp) - _shunt_charge_map(minus, dm)) ./ (2step)
    result = (2pi*vacuum_permittivity(Float64)) .*
             (charge_delta*solved.coefficients + solved.charge*delta)
    result[1, 1] += 2pi*vacuum_permittivity(Float64)*(inv(dp.k.Rtotal)-inv(dm.k.Rtotal))/(2step)
    return result
end

"""
$(TYPEDSIGNATURES)

Compute the internal capacitance of one dielectric group: the terminal capacitance of its
inner anchor and open conductors relative to its reference \\[F/m\\]. This is the
calculation of the `:boundary` formula. The Measurements extension adds a method for
uncertain groups that preserves their correlations without introducing uncertain scalars
into the dense arrays.

# Arguments

- `selected`: the boundary shunt formula.
- `group`: the dielectric group, built once by `dielectric_groups`.
- `directions=nothing`: optional columns of perturbations of the group's physical
  descriptors, in the order of `Vector(group)`, per unit dimensionless direction parameter.

# Returns

- Capacitance `C` \\[F/m\\] and the numerical record. With `directions`, also the
  implicit capacitance derivatives, one row per direction.
"""
function internal_capacitance(selected::Formula{:boundary}, group::DielectricGroup,
        directions = nothing)
    values = Float64.(Vector(group))
    numerical = DielectricGroup(group, values)
    solved = try
        _shunt_capacitance(selected, numerical, directions)
    catch exception
        exception isa BoundarySolveError || rethrow()
        throw(BoundarySolveError(exception.category,
            merge(exception.context, (; design = group.design, terminals = group.terminals)),
            exception.message))
    end
    directions === nothing && return (; solved.C, solved.diagnostic)
    # The residual of the nominal least-squares solve enters the derivative.
    (; plan, buffers) = solved.workspace
    _shunt_residual!(buffers.residual, plan, buffers, buffers.coefficients, buffers.excitation)
    tangents = Matrix{Float64}(undef, size(directions, 2), length(solved.C))
    for column in axes(directions, 2)
        direction = @view directions[:, column]
        relative = maximum(eachindex(direction)) do i
            abs(direction[i])/(iszero(values[i]) ? abs(values[2]) : abs(values[i]))
        end
        iszero(relative) && (tangents[column, :] .= 0; continue)
        step = cbrt(eps(Float64))/relative
        derivative = _shunt_tangent(selected, numerical, values, direction, solved, step/2)
        all(isfinite, derivative) || throw(BoundarySolveError(:derivative,
            (; design = group.design, terminals = group.terminals, column), "nonfinite sensitivity"))
        tangents[column, :] .= vec(derivative)
    end
    return (; solved.C, solved.diagnostic, tangents)
end
