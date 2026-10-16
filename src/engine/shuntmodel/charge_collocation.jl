# Integrated finite-face charge element for the explicit boundary shunt model.
# Coordinates are in m. Green functions use charge/(2pi*epsilon0).
# The buffers of the boundary shunt formula hold the dense arrays of one solve.

"""
$(TYPEDSIGNATURES)

Propagate the radial Dirichlet-to-Neumann load of one positive angular mode
through concentric dielectric layers. In log-radius coordinates, the load
transformation is the admittance form of the uniform-line input-impedance
relation [Sunde1968](@cite), Section 1.6, Eqs. (1.39)-(1.42), p. 15:

```math
Y_{\\mathrm{in}}=Y_c\\frac{Y_L+Y_c\\tanh(m\\ell)}
{Y_c+Y_L\\tanh(m\\ell)},\\qquad
Y_c=\\varepsilon_r m,\\qquad\\ell=\\log(r_o/r_i).
```

Here the loads ``Y_L,Y_c,Y_{\\mathrm{in}}`` are dimensionless modal loads,
not electrical admittances. The recursion starts with ``Y_L=\\infty`` at
grounded metal. Applying the line relation to cylindrical Laplace modes is
the log-radius construction used here, not a cable propagation calculation.

# Arguments

- `layers`: radially ordered layers with inner and outer radii \\[m\\] and
  positive relative permittivities \\[dimensionless\\].
- `m`: positive angular mode order.

# Keywords

- `reverse_layers=false`: reverse the layer order to propagate from the
  outer grounded shield toward the host instead of from the inner conductor.

# Returns

- Modal load seen at the host boundary \\[dimensionless\\]. `Inf` when the
  layer sequence is empty, representing an immediately adjacent metal boundary.
"""
function _shunt_load(layers, m; reverse_layers = false)
    load = Inf
    for layer in (reverse_layers ? Iterators.reverse(layers) : layers)
        t = tanh(m * log(layer.ro/layer.ri))
        characteristic = layer.epsilon * m
        load = isinf(load) ? characteristic/t :
               characteristic * (load + characteristic*t)/(characteristic + load*t)
    end
    load
end

"""
$(TYPEDSIGNATURES)

Evaluate the reflected Fourier modes of the local electrostatic Green function
in a concentric layered annulus. The radial basis ``r^{\\pm m}``, with
``1,\\log r`` for the zero mode, is the classical cylindrical Laplace basis.
See [Schelkunoff1934](@cite), Eqs. (122)-(124), p. 573, in the small-radius-to-
wavelength limit of cylindrical fields. Here it is used for the electrostatic
potential. The layered-kernel strategy is also supported by [Bernal1997](@cite),
Sections II-III, whose planar Prony/matrix-pencil kernel is not used here.

With ``t=\\log(r/a)``, each source-free dielectric layer satisfies
``d^2u_m/dt^2-m^2u_m=0``. Matching potential and normal displacement gives
the radial load recursion in `_shunt_load`, terminated at grounded metal.
The zero Fourier mode is evaluated separately in `_shunt_kernel`:

```math
G_0(r,s)=\\frac{R_< (R_{\\mathrm{total}}-R_>)}{R_{\\mathrm{total}}},
\\qquad R(r)=R_{\\mathrm{left}}+\\frac{\\log(r/a)}{\\varepsilon_h}.
```

Here ``R_<`` and ``R_>`` are the smaller and larger of ``R(r),R(s)``.
The dielectric log-radius sums ``R`` and host relative permittivity
``\\varepsilon_h`` are dimensionless. The log-radius normalization follows
the coaxial capacitance ``C=2\\pi\\varepsilon/\\log(r_o/r_i)`` \\[F/m\\], with
absolute permittivity ``\\varepsilon`` \\[F/m\\]. See [Sunde1968](@cite),
Section 1.5, Eq. (1.20), p. 11. Direct and nearest-interface logarithms
are extracted analytically. Only the remaining reflections are truncated.
For an adjacent dielectric, the extracted contrast is
``(\\varepsilon_h-\\varepsilon_j)/(\\varepsilon_h+\\varepsilon_j)``, also used
in [Campione2018](@cite), Section 2, Eq. (3). That paper uses a local planar
approximation to the braid geometry. It supports these dielectric image
factors, not this concentric annular Green function or its radial load recursion.

# Arguments

- `g`: host radii and layer boundaries \\[m\\], with positive relative
  permittivities \\[dimensionless\\].
- `modes`: number of retained positive Fourier modes.

# Returns

- `A`, `B`, `D`: smooth Fourier reflection coefficients \\[dimensionless\\].
- `ra0`, `rb0`: reflection coefficients at the nearest inner and outer interfaces
  \\[dimensionless\\]. The resulting kernel multiplies charge per unit length
  divided by ``2\\pi\\varepsilon_0`` \\[V\\] to give potential \\[V\\].
- `Rleft`, `Rtotal`: the zero-mode dielectric log-radius sums from the inner anchor to the
  host and to the reference \\[dimensionless\\].
"""
function _shunt_kernel_coefficients(g, modes)
    Rleft = sum(l->log(l.ro/l.ri)/l.epsilon, g.left; init = zero(g.a))
    Rright = sum(l->log(l.ro/l.ri)/l.epsilon, g.right; init = zero(g.a))
    Rtotal = Rleft+log(g.b/g.a)/g.epsilon+Rright
    A, B, D = zeros(modes), zeros(modes), zeros(modes)
    ra0 = isempty(g.left) ? -1.0 :
          (g.epsilon-last(g.left).epsilon)/(g.epsilon+last(g.left).epsilon)
    rb0 = isempty(g.right) ? -1.0 :
          (g.epsilon-first(g.right).epsilon)/(g.epsilon+first(g.right).epsilon)
    for m in 1:modes
        l = _shunt_load(g.left, m)/g.epsilon
        r = _shunt_load(g.right, m; reverse_layers = true)/g.epsilon
        ra = isinf(l) ? -1.0 : (m-l)/(m+l)
        rb = isinf(r) ? -1.0 : (m-r)/(m+r)
        denominator = g.epsilon*m*(1-ra*rb*(g.a/g.b)^(2m))
        A[m] = ra/denominator-ra0/(g.epsilon*m)
        B[m] = rb/denominator-rb0/(g.epsilon*m)
        D[m] = ra*rb/denominator
    end
    (; A, B, D, ra0, rb0, Rleft, Rtotal)
end

# The direct, nearest-image and zero-mode terms of the Green function. The reflected
# Fourier modes enter by the batched products of `_shunt_kernel_matrix!`.
function _shunt_kernel(z, source, g, k; split_images = false)
    r, s = abs(z), abs(source)
    Rr = k.Rleft + log(r/g.a)/g.epsilon
    Rs = k.Rleft + log(s/g.a)/g.epsilon
    result = min(Rr, Rs)*(k.Rtotal-max(Rr, Rs))/k.Rtotal
    if split_images
        result += (log(max(r, s))+k.ra0*log(s)+k.rb0*(2log(g.b)-log(r)))/g.epsilon
    else
        result += log(max(r, s)/abs(z-source))/g.epsilon
        result -= k.ra0/g.epsilon*log(abs(1-g.a^2/(conj(z)*source)))
        result -= k.rb0/g.epsilon*log(abs(1-z*conj(source)/g.b^2))
    end
    result
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the arrays that solving `group` uses, sized by `plan`. The `N` ring
wires have `r` targets and `c` sources each. The border consists of the other `mb` targets
and `nb` charge coefficients:

- `ring_kernel` holds the potentials of the ring's sources at the first ring wire's
  targets, `coupling_kernel` at the border targets, and `border` the potentials of the
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
- `target_inner`, `target_outer`, `source_inner` and `source_outer` hold the image ratios
  ``a/\\bar z`` and ``z/b`` of the targets and of the sources or face nodes. The `*_power`
  arrays hold their running powers.
- `target_modes` and `source_modes` hold the Fourier factors of at most 64 modes at a time.
- `face_kernel` holds the Green function from the targets to the nodes of one face.
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

function _shunt_kernel_matrix!(matrix, targets, sources, g, k; split_images = false, buffers)
    # Same Green function, batched into small Fourier blocks for BLAS. No dense
    # N-by-modes array: the products use only 4*64 columns per geometric boundary set.
    @inbounds for j in eachindex(sources), i in eachindex(targets)

        matrix[i, j] = _shunt_kernel(targets[i], sources[j], g, k; split_images)
    end
    nr, nc = length(targets), length(sources)
    ta, tb = @view(buffers.target_inner[1:nr]), @view(buffers.target_outer[1:nr])
    sa, sb = @view(buffers.source_inner[1:nc]), @view(buffers.source_outer[1:nc])
    tpa = @view buffers.target_inner_power[1:nr]
    tpb = @view buffers.target_outer_power[1:nr]
    spa = @view buffers.source_inner_power[1:nc]
    spb = @view buffers.source_outer_power[1:nc]
    ta .= g.a ./ conj.(targets);
    tb .= targets ./ g.b
    sa .= g.a ./ conj.(sources);
    sb .= sources ./ g.b
    fill!(tpa, 1);
    fill!(tpb, 1);
    fill!(spa, 1);
    fill!(spb, 1)
    U, V = @view(buffers.target_modes[1:nr, :]), @view(buffers.source_modes[1:nc, :])
    for first_mode in 1:64:length(k.A)
        modes = first_mode:min(first_mode + 63, length(k.A))
        for (column, m) in enumerate(modes)
            tpa .*= ta
            tpb .*= tb
            spa .*= sa
            spb .*= sb
            cross = k.D[m]*(g.a/g.b)^m
            @inbounds for i in eachindex(targets)
                U[i, 4column - 3]=real(tpa[i])
                U[i, 4column - 2]=imag(tpa[i])
                U[i, 4column - 1]=real(tpb[i])
                U[i, 4column]=imag(tpb[i])
            end
            @inbounds for j in eachindex(sources)
                a=k.A[m]*spa[j]+cross*spb[j]
                b=k.B[m]*spb[j]+cross*spa[j]
                V[j, 4column - 3]=real(a)
                V[j, 4column - 2]=imag(a)
                V[j, 4column - 1]=real(b)
                V[j, 4column]=imag(b)
            end
        end
        @views mul!(
            matrix, U[:, 1:4length(modes)], transpose(V[:, 1:4length(modes)]), 1.0, 1.0)
    end
    matrix
end

_shunt_core_voltage(z, g, k) = 1-(k.Rleft+log(abs(z)/g.a)/g.epsilon)/k.Rtotal

# Point `j` of `count` on the circle at `fraction` of a wire's radius, in the wire's own
# polar frame. A rotation that maps one ring wire onto the next maps its points too.
# Round wires keep auxiliary sources inside. The tape has no source contour or inset: its
# physical faces are integrated below.
function _shunt_wire_point(wire, fraction, j, count)
    complex(wire.x, wire.y) + fraction*wire.r*cis(atan(wire.y, wire.x) + 2pi*j/count)
end

"""
$(TYPEDSIGNATURES)

Order the wires of `g` for collocation and find its ring. The wires of a ring share one
terminal, one radius and one pitch circle, and they are equally spaced. Each exposes the same
targets, at least as many as its sources. The ring's wires come first, in angular order.

# Returns

- `wires`: the wire indices in collocation order.
- `ring`: the number of ring wires, zero without a ring.
- `exposed`: for each wire in that order, the indices of its `2level.wire` surface points
  that lie outside every other wire of its terminal.
"""
function _shunt_layout(g, level)
    count = 2level.wire
    exposed = map(g.wires) do wire
        filter(0:(count - 1)) do j
            point = _shunt_wire_point(wire, 1, j, count)
            # Contact points and buried arcs do not expose dielectric to this terminal.
            all(other -> other === wire || other.terminal != wire.terminal ||
                abs(point-complex(other.x, other.y))-other.r > _shunt_tol(g.b), g.wires)
        end
    end
    angles = [mod(atan(wire.y, wire.x), 2pi) for wire in g.wires]
    order = sortperm(angles)
    n = length(order)
    ring = n > 0 && length(exposed[first(order)]) >= level.wire && all(1:n) do i
        wire, lead = g.wires[order[i]], g.wires[first(order)]
        pitch = hypot(wire.x, wire.y)
        spacing = angles[order[i]] - angles[first(order)] - 2pi*(i - 1)/n
        wire.terminal == lead.terminal && abs(wire.r - lead.r) <= _shunt_tol(g.b) &&
            abs(pitch - hypot(lead.x, lead.y)) <= _shunt_tol(g.b) &&
            abs(pitch*spacing) <= _shunt_tol(g.b) && exposed[order[i]] == exposed[first(order)]
    end
    wires = ring ? order : collect(eachindex(g.wires))
    return (; wires, ring = ring ? n : 0, exposed = exposed[wires])
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
dielectric interface. The host occupies angle ``\\pi/2`` and the adjacent
dielectric angle ``\\pi``. Continuity of potential and normal displacement,
with constant potential on the two metal faces, gives this specialization:

```math
\\varepsilon_h\\cos(\\nu\\pi/2)\\sin(\\nu\\pi)
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

- `epsilon_host`: host relative permittivity \\[dimensionless\\].
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
function _shunt_junction_exponent(epsilon_host, epsilon_outer)
    epsilon_host > 0 && epsilon_outer > 0 || throw(BoundarySolveError(:unsupported,
        (; epsilon_host, epsilon_outer), "boundary shunt requires positive permittivities"))
    # Host sector pi/2, outer dielectric sector pi. Determinant has no cotangent
    # poles. Exclude the trivial nu=0 root and bracket the first positive root.
    determinant(nu) = epsilon_host*cospi(nu/2)*sinpi(nu) +
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
# holds its geometry, its endpoint exponents, quadrature, weighted Jacobi values and the
# Jacobi recurrence coefficients.
const ShuntArc = NamedTuple{(:radius, :phi, :span), NTuple{3, Float64}}
const ShuntEnd = NamedTuple{(:phi, :mid, :half), NTuple{3, Float64}}
const ShuntFace{S} = NamedTuple{
    (:shape, :alpha, :beta, :p, :nodes, :weights, :weighted, :beta_norm, :terminal, :points,
        :recurrence),
    Tuple{S, Float64, Float64, Int, Vector{Float64}, Vector{Float64}, Matrix{Float64},
        Float64, Int, Vector{ComplexF64}, Vector{NTuple{3, Float64}}}}
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

- `g`: host and layer radii \\[m\\] and relative permittivities \\[dimensionless\\].
- `s`: tape geometry and terminal index. The inner and outer radii are in
  \\[m\\], and the center angle and angular span are in \\[rad\\].
- `p`: highest Jacobi degree on each face.
- `quadrature`: number of normalized Gauss–Jacobi nodes per face.

# Returns

- The 4 face records of the tape, with physical coordinates \\[m\\], endpoint exponents,
  quadrature weights and weighted polynomial values.
"""
function _shunt_tape_faces(g, s, p, quadrature)
    # A metal face on a dielectric interface meets that dielectric at its corners.
    at_interface = abs(s.ro - g.b) <= _shunt_tol(g.b)
    at_inner = abs(s.ri - g.a) <= _shunt_tol(g.b)
    nu_inner = at_inner ? _shunt_junction_exponent(g.epsilon, last(g.left).epsilon) : 2/3
    nu_outer = at_interface ? _shunt_junction_exponent(g.epsilon, first(g.right).epsilon) :
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
                recurrence = _shunt_jacobi_recurrence(face.alpha, face.beta, p)))
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

# Add the integrated charge columns of every tape face to `columns`. Count the moments
# whose QuadGK estimate exceeds the requested target.
function _shunt_tape_columns!(columns, targets, faces, g, k; integration, buffers)
    estimated_error = 0.0
    unmet = 0
    offset = 0
    for tape in faces, face in tape
        block = @view columns[:, (offset + 1):(offset + face.p + 1)]
        kernel = @view buffers.face_kernel[1:length(targets), 1:length(face.points)]
        _shunt_kernel_matrix!(
            kernel, targets, face.points, g, k; split_images = true, buffers)
        mul!(block, kernel, face.weighted)
        for (i, z) in pairs(targets)
            inner, outer = g.a^2/conj(z), g.b^2/conj(z)
            direct_coefficient = 1.0
            inner_same = abs(inner-z) <= _shunt_tol(g.b)
            outer_same = abs(outer-z) <= _shunt_tol(g.b)
            direct_coefficient += (inner_same ? k.ra0 : 0.0) + (outer_same ? k.rb0 : 0.0)
            for (point, coefficient) in ((inner, inner_same ? 0.0 : k.ra0),
                (outer, outer_same ? 0.0 : k.rb0), (z, direct_coefficient))
                coefficient == 0 && continue
                moments,
                err = _shunt_log_moments(point, face; integration...,
                    result = buffers.moments, segments = buffers.segments)
                @inbounds for n in eachindex(moments)
                    block[i, n] -= (coefficient/g.epsilon)*moments[n]
                end
                estimated_error = max(estimated_error, abs(coefficient/g.epsilon)*err)
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
together and the border. `ring` counts the wires that the Fourier reduction separates, zero
without a ring. `unknowns` and `equations` count the charge coefficients and collocation
equations, and `level` is the resolution.
"""
const InternalShuntDiagnostic = NamedTuple{
    (:log_moment_error, :condition, :ring, :unknowns, :equations, :level),
    Tuple{Float64, Float64, Int, Int, Int, typeof(DEFAULT_RESOLUTION)}}

# The fixed geometry of one solve in the wire order of `layout`: tape faces, collocation
# targets and their terminals, auxiliary sources and Fourier load coefficients, with the
# array sizes they give. The ring's targets and sources come first, wire by wire, then the
# other wires' and the tape faces'.
function _shunt_plan(selected::Formula{:boundary}, g, layout, derivatives = 0)
    (; resolution, integration) = selected.options.data
    level = resolution
    faces = ShuntTapeFaces[_shunt_tape_faces(g, tape, level.order, level.quadrature)
                           for tape in g.tapes]
    points, terminals, sources = ComplexF64[], Int[], ComplexF64[]
    for (index, exposed) in zip(layout.wires, layout.exposed)
        wire = g.wires[index]
        append!(points, _shunt_wire_point(wire, 1, j, 2level.wire) for j in exposed)
        append!(terminals, fill(wire.terminal, length(exposed)))
        append!(sources,
            _shunt_wire_point(wire, 0.75, j, level.wire) for j in 0:(level.wire - 1))
    end
    nface = max(24, 4(level.order+1))
    parameters = -cos.(pi .* ((1:nface) .- 0.5) ./ nface)
    for tape in faces, face in tape
        append!(points, [_shunt_face_point(face.shape, t) for t in parameters])
        append!(terminals, fill(face.terminal, nface))
    end
    k = _shunt_kernel_coefficients(g, level.modes)
    columns = length(sources)+4(level.order+1)*length(faces)
    ports = max(maximum(wire -> wire.terminal, g.wires; init = 1),
        maximum(tape -> tape.terminal, g.tapes; init = 1))
    ring = layout.ring
    return (; faces, points, terminals, sources, k, level, integration, layout, ring,
        block_rows = iszero(ring) ? 0 : length(first(layout.exposed)),
        block_columns = iszero(ring) ? 0 : level.wire,
        rows = length(points), columns, ports, modes = level.modes,
        quadrature = isempty(faces) ? 0 : level.quadrature, order = level.order, derivatives)
end

function _shunt_rhs!(rhs, points, terminals, g, k)
    fill!(rhs, 0)
    @inbounds for i in eachindex(points)
        rhs[i, 1] = -_shunt_core_voltage(points[i], g, k)
        rhs[i, terminals[i]] = 1
    end
    return rhs
end

function _shunt_matrix!(matrix, points, sources, faces, g, k; integration, buffers)
    count = length(sources)
    _shunt_kernel_matrix!(@view(matrix[:, 1:count]), points, sources, g, k; buffers)
    return _shunt_tape_columns!(
        @view(matrix[:, (count + 1):end]), points, faces, g, k; integration, buffers)
end

# Total charge of each terminal and the charge moment of the core potential, per face.
function _shunt_charge!(charge, offset, face::ShuntFace{ShuntArc}, g, k)
    charge[face.terminal, offset + 1] = 1
    # The core potential is constant on an arc: only the uniform mode has a moment.
    charge[1, offset + 1] = -_shunt_core_voltage(first(face.points), g, k)
    return charge
end

function _shunt_charge!(charge, offset, face::ShuntFace{ShuntEnd}, g, k)
    charge[face.terminal, offset + 1] = 1
    @views mul!(charge[1:1, (offset + 1):(offset + face.p + 1)],
        transpose(_shunt_core_voltage.(face.points, Ref(g), Ref(k))), face.weighted, -1.0, 0.0)
    return charge
end

function _shunt_charge_map(g, plan)
    (; sources, faces, level, k) = plan
    charge = zeros(plan.ports, plan.columns)
    @inbounds for (position, index) in pairs(plan.layout.wires)
        for j in ((position - 1) * level.wire + 1):(position * level.wire)
            charge[1, j] = -_shunt_core_voltage(sources[j], g, k)
            charge[g.wires[index].terminal, j] = 1
        end
    end
    offset = length(sources)
    for tape in faces, face in tape
        _shunt_charge!(charge, offset, face, g, k)
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
    (; faces, points, terminals, sources, k, integration) = plan
    N, r, c = plan.ring, plan.block_rows, plan.block_columns
    m, n = plan.rows, plan.columns
    m >= n || throw(BoundarySolveError(:rank, (; rows = m, columns = n),
        "fewer boundary equations than charge coefficients"))
    (; ring_kernel, coupling_kernel, border, column_norms) = buffers
    kept = size(coupling_kernel, 1)
    # The ring needs one block row of its kernel, and its columns at the border rows.
    ring = @view sources[1:(N * c)]
    _shunt_kernel_matrix!(ring_kernel, @view(points[1:r]), ring, g, k; buffers)
    _shunt_kernel_matrix!(
        coupling_kernel, @view(points[(N * r + 1):(N * r + kept)]), ring, g, k; buffers)
    tape = _shunt_matrix!(border, points, @view(sources[(N * c + 1):end]), faces, g, k;
        integration, buffers)
    tape.unmet == 0 ||
        @warn "Boundary logarithmic quadrature target was not met" integrals=tape.unmet log_moment_error=tape.log_moment_error rtol=integration.rtol atol=integration.atol maxevals=integration.maxevals
    # A ring column's norm over the ring rows is the same on every wire, which keeps the
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
    rhs = _shunt_rhs!(buffers.excitation, points, terminals, g, k)
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
        _shunt_matrix!(ap, xp, dp.sources, dp.faces, plus, dp.k; plan.integration, buffers)
        _shunt_matrix!(am, xm, dm.sources, dm.faces, minus, dm.k; plan.integration, buffers)
        ap .-= am
        ap ./= 2step
        rp = @view buffers.plus_excitation[rows, :]
        rm = @view buffers.minus_excitation[rows, :]
        _shunt_rhs!(rp, xp, @view(dp.terminals[start:stop]), plus, dp.k)
        _shunt_rhs!(rm, xm, @view(dm.terminals[start:stop]), minus, dm.k)
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
