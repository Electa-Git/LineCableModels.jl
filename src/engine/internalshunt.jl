# Integrated finite-face charge element, using the accepted dev prototype.
# Coordinates are in m; Green functions use charge/(2pi*epsilon0).
# Numerical arrays are local to a preparation, never a process-global cache.

"""
$(TYPEDSIGNATURES)

Propagate the radial Dirichlet-to-Neumann load of one positive angular mode
through concentric dielectric layers. In log-radius coordinates, the load
transformation is the admittance form of the uniform-line input-impedance
relation [Sunde1968](@cite), Section 1.6, Eqs. (1.39)–(1.42), p. 15:

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

- `layers`: Radially ordered layers with inner/outer radii \\[m\\] and
  positive relative permittivities \\[dimensionless\\].
- `m`: Positive angular mode order.

# Keywords

- `reverse_layers=false`: Reverse the layer order to propagate from the
  outer grounded shield toward the host instead of from the inner conductor.

# Returns

- Modal load seen at the host boundary \\[dimensionless\\]; `Inf` when the
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

Prepare the reflected Fourier modes of the local electrostatic Green function
in a concentric layered annulus. The radial basis ``r^{\\pm m}``, with
``1,\\log r`` for the zero mode, is the classical cylindrical Laplace basis;
see [Schelkunoff1934](@cite), Eqs. (122)–(124), p. 573, in the small-radius-to-
wavelength limit of cylindrical fields. Here it is used for the electrostatic
potential. The layered-kernel strategy is also supported by [Bernal1997](@cite),
Sections II–III, whose planar Prony/matrix-pencil kernel is not used here.

With ``t=\\log(r/a)``, each source-free dielectric layer satisfies
``d^2u_m/dt^2-m^2u_m=0``. Matching potential and normal displacement gives
the radial load recursion in `_shunt_load`, terminated at grounded metal.
The zero Fourier mode is evaluated separately in `_shunt_kernel`:

```math
G_0(r,s)=\\frac{R_< (R_{\\mathrm{total}}-R_>)}{R_{\\mathrm{total}}},
\\qquad R(r)=R_{\\mathrm{left}}+\\frac{\\log(r/a)}{\\varepsilon_h}.
```

Here ``R_<`` and ``R_>`` are the smaller and larger of ``R(r),R(s)``;
the dielectric log-radius sums ``R`` and host relative permittivity
``\\varepsilon_h`` are dimensionless. The log-radius normalization follows
the coaxial capacitance ``C=2\\pi\\varepsilon/\\log(r_o/r_i)`` \\[F/m\\], with
absolute permittivity ``\\varepsilon`` \\[F/m\\]; see [Sunde1968](@cite),
Section 1.5, Eq. (1.20), p. 11. Direct and nearest-interface logarithms
are extracted analytically; only the remaining reflections are truncated.
For an adjacent dielectric, the extracted contrast is
``(\\varepsilon_h-\\varepsilon_j)/(\\varepsilon_h+\\varepsilon_j)``, also used
in [Campione2018](@cite), Section 2, Eq. (3). That paper uses a local planar
approximation to the braid geometry; it supports these dielectric image
factors, not this concentric annular Green function or its radial load recursion.

# Arguments

- `g`: Host radii and layer boundaries \\[m\\], positive relative
  permittivities \\[dimensionless\\], and dielectric log-radius sums.
- `modes`: Number of retained positive Fourier modes.

# Returns

- `A`, `B`, `D`: Smooth Fourier reflection coefficients \\[dimensionless\\].
- `ra0`, `rb0`: Nearest inner/outer interface reflection coefficients
  \\[dimensionless\\]. The resulting kernel multiplies charge per unit length
  divided by ``2\\pi\\varepsilon_0`` \\[V\\] to give potential \\[V\\].
"""
function _shunt_kernel_coefficients(g, modes)
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
    (; A, B, D, ra0, rb0)
end

function _shunt_kernel(z, source, g, k; regular = false, split_images = false)
    r, s = abs(z), abs(source)
    Rr = g.Rleft + log(r/g.a)/g.epsilon
    Rs = g.Rleft + log(s/g.a)/g.epsilon
    result = min(Rr, Rs)*(g.Rtotal-max(Rr, Rs))/g.Rtotal
    if split_images
        result += (log(max(r,s))+k.ra0*log(s)+k.rb0*(2log(g.b)-log(r)))/g.epsilon
    else
        result += (regular ? log(max(r, s)) : log(max(r, s)/abs(z-source)))/g.epsilon
        result -= k.ra0/g.epsilon*log(abs(1-g.a^2/(conj(z)*source)))
        result -= k.rb0/g.epsilon*log(abs(1-z*conj(source)/g.b^2))
    end
    qa, qb = g.a^2/(r*s), r*s/g.b^2
    qd = (g.a/g.b)^2
    qc, qe = qd*r/s, qd*s/r
    pa, pb, pc, pe = qa, qb, qc, qe
    cosine = real(z*conj(source))/(r*s)
    previous, current = 1.0, cosine
    @inbounds for m in eachindex(k.A)
        result += current*(k.A[m]*pa + k.B[m]*pb + k.D[m]*(pc+pe))
        previous, current = current, 2cosine*current-previous
        pa *= qa
        pb *= qb
        pc *= qc
        pe *= qe
    end
    result
end

function _shunt_kernel_matrix!(matrix, targets, sources, g, k; regular = false, split_images = false)
    # Same Green function, batched into small Fourier blocks for BLAS. No dense
    # N-by-modes cache: working storage is only 4*64 columns per boundary set.
    zero_modes = merge(k, (; A = Float64[], B = Float64[], D = Float64[]))
    @inbounds for j in eachindex(sources), i in eachindex(targets)
        matrix[i,j] = _shunt_kernel(targets[i], sources[j], g, zero_modes; regular, split_images)
    end
    ta, tb = g.a ./ conj.(targets), targets ./ g.b
    sa, sb = g.a ./ conj.(sources), sources ./ g.b
    tpa, tpb = ones(ComplexF64, length(targets)), ones(ComplexF64, length(targets))
    spa, spb = ones(ComplexF64, length(sources)), ones(ComplexF64, length(sources))
    U = Matrix{Float64}(undef, length(targets), 4min(64,length(k.A)))
    V = Matrix{Float64}(undef, length(sources), size(U,2))
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
        @views mul!(matrix, U[:,1:4length(modes)], transpose(V[:,1:4length(modes)]), 1.0, 1.0)
    end
    matrix
end

function _shunt_kernel_matrix(targets, sources, g, k; kwargs...)
    return _shunt_kernel_matrix!(Matrix{Float64}(undef,length(targets),length(sources)),
        targets, sources, g, k; kwargs...)
end

_shunt_core_voltage(z, g) = 1-(g.Rleft+log(abs(z)/g.a)/g.epsilon)/g.Rtotal

# Round wires retain their existing auxiliary-source treatment. The tape has
# no source contour or inset: its physical faces are integrated below.
function _shunt_points(g, nw, fraction; shift = 0.0)
    targets, sources = ComplexF64[], ComplexF64[]
    for wire in g.wires
        centre = complex(wire.x, wire.y)
        for j in 0:(nw - 1)
            push!(targets, centre + wire.r*cis(2pi*(j+shift)/nw))
            push!(sources, centre + fraction*wire.r*cis(2pi*j/nw))
        end
    end
    (; targets, sources)
end

function _shunt_jacobi!(values,t,alpha,beta)
    p = length(values)-1
    values[1] = 1
    p == 0 && return values
    values[2] = ((alpha+beta+2)*t+alpha-beta)/2
    s = alpha+beta
    for n in 1:p-1
        A = (2n+s+1)*(2n+s+2)/(2(n+1)*(n+s+1))
        B = (alpha^2-beta^2)*(2n+s+1)/(2(n+1)*(n+s+1)*(2n+s))
        C = (n+alpha)*(n+beta)*(2n+s+2)/((n+1)*(n+s+1)*(2n+s))
        values[n+2] = (A*t+B)*values[n+1]-C*values[n]
    end
    values
end

_shunt_jacobi(t,alpha,beta,p) = _shunt_jacobi!(Vector{Float64}(undef,p+1),t,alpha,beta)

function _shunt_gauss_jacobi(n,alpha,beta)
    alpha > -1 && beta > -1 && n > 1 || error("Invalid Jacobi rule.")
    s = alpha+beta
    diagonal = [(beta-alpha)/(s+2);
        [(beta^2-alpha^2)/((2k+s)*(2k+s+2)) for k in 1:n-1]]
    off = [2/(s+2)*sqrt((1+alpha)*(1+beta)/(s+3));
        [2/(2k+s)*sqrt(k*(k+alpha)*(k+beta)*(k+s)/
            ((2k+s-1)*(2k+s+1))) for k in 2:n-1]]
    decomposition = eigen(SymTridiagonal(diagonal,off))
    # Normalized measure w(t)dt / integral(w); no arclength factor belongs here.
    decomposition.values,vec(decomposition.vectors[1,:].^2)
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
Material-dependent metal–dielectric wedge singularities, rather than a
universal homogeneous exponent, are treated by [VanBladel1985](@cite).
For equal permittivities this equation recovers ``\\nu=2/3``.

# Arguments

- `epsilon_host`: Host relative permittivity \\[dimensionless\\].
- `epsilon_outer`: Adjacent dielectric relative permittivity \\[dimensionless\\].

# Returns

- Leading exponent ``\\nu`` \\[dimensionless\\], found by bisection excluding
  the trivial zero root. This is the specified two-dielectric corner, not a
  general wedge solver.

# Errors

Nonpositive permittivities or an unbracketed root raise `ErrorException`.
"""
function _shunt_junction_exponent(epsilon_host,epsilon_outer)
    epsilon_host > 0 && epsilon_outer > 0 || error("Positive dielectric permittivities required.")
    # Host sector pi/2, outer dielectric sector pi. Determinant has no cotangent
    # poles. Exclude the trivial nu=0 root and bracket the first positive root.
    determinant(nu) = epsilon_host*cospi(nu/2)*sinpi(nu) +
        epsilon_outer*sinpi(nu/2)*cospi(nu)
    left,right = 1e-10,1.0
    determinant(left) > 0 && determinant(right) < 0 ||
        error("internal shunt: dielectric corner root is not bracketed")
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

_shunt_face_point(face,t) = face.kind == :arc ?
    face.radius*cis(face.phi+face.span*t/2) :
    (face.mid+face.half*t)*cis(face.phi)

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
\\[V\\]; ``dq`` is charge per cable length \\[C/m\\]. Each face's total charge
is ``2\\pi\\varepsilon_0 a_0``. Do not multiply this measure by a second
arclength Jacobian. Encoding known singularities in the approximation space
is also supported by [Classen2011](@cite), Sections 2 and 4; their FIT/DG
discretizations are not the boundary-integral element implemented here.

# Arguments

- `g`: Host/layer radii \\[m\\] and relative permittivities \\[dimensionless\\].
- `s`: Tape inner/outer radii \\[m\\], centre angle and angular span \\[rad\\],
  and terminal index.
- `p`: Highest Jacobi degree on each face.
- `quadrature`: Number of normalized Gauss–Jacobi nodes per face.

# Returns

- Four face records with physical coordinates \\[m\\], endpoint exponents,
  quadrature weights and weighted polynomial values.

# Errors

Insufficient quadrature order or tape contact with a bounding conductor
raises `ErrorException`.
"""
function _shunt_tape_faces(g,s,p,quadrature)
    p >= 0 && quadrature >= max(2,p+1) || error("Quadrature needs at least order+1 nodes.")
    at_interface = isapprox(s.ro,g.b; atol = 64eps(Float64)*g.b,rtol = 0)
    at_inner = isapprox(s.ri,g.a; atol = 64eps(Float64)*g.b,rtol = 0)
    at_inner && isempty(g.left) && error("internal shunt: tape contacts the inner conductor")
    nu_inner = at_inner ? _shunt_junction_exponent(g.epsilon,last(g.left).epsilon) : 2/3
    # A metal face coincident with the reference shield is not an open tape.
    at_interface && isempty(g.right) && error("internal shunt: tape contacts the reference shield")
    nu_outer = at_interface ? _shunt_junction_exponent(g.epsilon,first(g.right).epsilon) : 2/3
    descriptions = [
        (;kind = :arc,radius = s.ri,phi = s.phi,span = s.span,mid = 0.0,half = 0.0,
            alpha = nu_inner-1,beta = nu_inner-1),
        (;kind = :arc,radius = s.ro,phi = s.phi,span = s.span,mid = 0.0,half = 0.0,
            alpha = nu_outer-1,beta = nu_outer-1),
        [(;kind = :end,radius = 0.0,phi = s.phi+sign*s.span/2,span = 0.0,
            mid = (s.ri+s.ro)/2,half = (s.ro-s.ri)/2,
            alpha = nu_outer-1,beta = nu_inner-1) for sign in (-1,1)]...]
    map(descriptions) do face
        nodes,weights = _shunt_gauss_jacobi(quadrature,face.alpha,face.beta)
        polys = transpose(reduce(hcat,[_shunt_jacobi(t,face.alpha,face.beta,p) for t in nodes]))
        weighted = weights.*polys
        loggamma = SpecialFunctions.loggamma
        beta_norm = exp(loggamma(face.alpha+1)+loggamma(face.beta+1)-
            loggamma(face.alpha+face.beta+2))
        merge(face,(;p,nodes,weights,weighted,beta_norm,terminal=s.terminal,
            points = _shunt_face_point.(Ref(face),nodes)))
    end
end

function _shunt_face_projection(z,face)
    if face.kind == :arc
        radius = abs(z)
        angle_difference = atan(sin(angle(z)-face.phi),cos(angle(z)-face.phi))
        u = 2angle_difference/face.span
        delta = abs(radius-face.radius)/(face.radius*face.span/2)
    else
        rotated = z*cis(-face.phi)
        u = (real(rotated)-face.mid)/face.half
        delta = abs(imag(rotated))/face.half
    end
    u,delta
end

"""
$(TYPEDSIGNATURES)

Integrate the logarithmic kernel against the normalized face charge modes:

```math
M_n(z)=\\int_{-1}^{1} P_n^{(\\alpha,\\beta)}(t)
\\log\\!\\left(\\frac{|z-\\zeta(t)|}{1\\,\\mathrm{m}}\\right)d\\mu(t).
```

The face position is ``\\zeta(t)`` \\[m\\]; ``d\\mu`` is the normalized measure
in `_shunt_tape_faces`. Ordinary quadrature can lose accuracy close to a
source boundary [HelsingOjala2008](@cite), Section 1. Here direct and image
logarithms are integrated separately from the smooth Green remainder.
Nearby targets use ``t=\\cos\\theta`` and adaptive Gauss–Kronrod integration,
split at the projected singularity or near peak; distant targets use the
face's Gauss–Jacobi rule. This is not Helsing–Ojala's rational-quadrature
algorithm, and that reference does not prescribe the corner basis.

# Arguments

- `z`: Target coordinate in the complex cross-sectional plane \\[m\\].
- `face`: Physical face geometry, normalized weights and Jacobi degree.

# Keywords

- `rtol=1e-10`: Relative adaptive quadrature tolerance; the absolute tolerance
  is `0.01rtol` for these dimensionless moments.
- `result=zeros(face.p+1)`: In-place moment workspace, also returned.
- `segments=nothing`: Optional reusable QuadGK segment buffer.

# Returns

- Logarithmic moments \\[dimensionless\\].
- Adaptive error estimate; the far-field fixed rule returns zero because it
  supplies no estimate, not because its error is known to vanish.

# Errors

An unresolved zero-distance evaluation or failed adaptive quadrature raises
`ErrorException`; the logarithm is not replaced by an arbitrary finite floor.
"""
function _shunt_log_moments(z,face; rtol = 1e-10,
        result = zeros(face.p+1), segments = nothing)
    u,delta = _shunt_face_projection(z,face)
    if hypot(max(abs(u)-1,0),delta) > 0.2
        fill!(result,0)
        @inbounds for j in eachindex(face.points)
            value = log(abs(face.points[j]-z))
            for n in eachindex(result)
                result[n] += face.weighted[j,n]*value
            end
        end
        return result,0.0
    end
    # Weighted logarithmic product integration: compute the log moments of the
    # Jacobi modes, independently of the smooth Green remainder. t=cos(theta)
    # absorbs endpoint weights; split at the logarithmic singularity/near peak.
    theta0 = acos(clamp(u,-1,1))
    splits = [0.0,theta0,pi]
    if delta > 0
        width = min(sqrt(delta),delta/max(sin(theta0),eps(Float64)))
        for sign in (-1,1), factor in (1,4)
            push!(splits,clamp(theta0+sign*factor*width,0,pi))
        end
    end
    sort!(unique!(splits))
    function integrand!(output,theta)
        t = cos(theta)
        difference = abs(u) <= 1 ?
            -2sin((theta+theta0)/2)*sin((theta-theta0)/2) : t-u
        distance = if face.kind == :arc
            # Exact opposed-arc distance, also valid in the self limit.
            hypot(abs(z)-face.radius,
                2sqrt(abs(z)*face.radius)*sin(face.span*difference/4))
        else
            hypot(face.half*difference,imag(z*cis(-face.phi)))
        end
        # Roundoff may identify the endpoint only after its contribution is
        # below integration accuracy; never manufacture a finite log(0) floor.
        distance > 0 || error("Unresolved logarithmic integration point.")
        weight = sin(theta/2)^(2face.alpha+1)*cos(theta/2)^(2face.beta+1)/face.beta_norm
        _shunt_jacobi!(output,t,face.alpha,face.beta)
        output .*= log(distance)*weight
    end
    # Reuse quadrature work vectors; modal integration must not allocate a new
    # polynomial vector at each of its thousands of function evaluations.
    integral, estimate = quadgk!(integrand!,result,splits;
        segbuf=segments,rtol,atol = rtol*0.01,order = max(7,cld(face.p+1,2)),maxevals = 100000)
    estimate <= max(rtol*norm(integral),rtol*0.01) ||
        error("internal shunt: logarithmic quadrature did not converge")
    return integral, estimate
end

function _shunt_tape_columns!(columns,targets,faces,g,k; log_rtol = 1e-10)
    error_bound = 0.0
    offset = 0
    for face in faces
        block = @view columns[:,offset+1:offset+face.p+1]
        mul!(block,_shunt_kernel_matrix(targets,face.points,g,k; split_images = true),face.weighted)
        result = zeros(face.p+1)
        segments = alloc_segbuf(Float64, Vector{Float64}, Float64; size=32)
        for (i,z) in pairs(targets)
            inner,outer = g.a^2/conj(z),g.b^2/conj(z)
            direct_coefficient = 1.0
            images = Tuple{ComplexF64,Float64}[]
            for (point,coefficient) in ((inner,k.ra0),(outer,k.rb0))
                if abs(point-z) <= 8eps(Float64)*g.b
                    direct_coefficient += coefficient
                else
                    push!(images,(point,coefficient))
                end
            end
            push!(images,(z,direct_coefficient))
            for (point,coefficient) in images
                coefficient == 0 && continue
                moments,err = _shunt_log_moments(point,face; rtol = log_rtol, result, segments)
                @inbounds for n in eachindex(moments)
                    block[i,n] -= (coefficient/g.epsilon)*moments[n]
                end
                error_bound = max(error_bound,abs(coefficient/g.epsilon)*err)
            end
        end
        offset += face.p+1
    end
    (;columns,log_moment_error = error_bound)
end

function _shunt_tape_columns(targets,faces,g,k; kwargs...)
    columns = Matrix{Float64}(undef,length(targets),sum(f.p+1 for f in faces;init=0))
    return _shunt_tape_columns!(columns,targets,faces,g,k; kwargs...)
end

function _shunt_face_targets(faces,n; validation = false)
    parameters = -cos.(pi.*((1:n).-0.5)./n)
    if validation
        # Independent uniform targets plus progressively closer corner targets.
        parameters = sort!(unique!([parameters;collect(range(-1,1;length = n+1));
            [-1+10.0^-k for k in 2:10];[1-10.0^-k for k in 2:10]]))
    end
    [_shunt_face_point(face,t) for face in faces for t in parameters]
end

"Accepted finite-strip discretization; controls change numerical resolution, not geometry."
const INTERNAL_SHUNT_RESOLUTION = (wire=64, order=32, quadrature=256, modes=1024)
"Maximum dense collocation storage per local solve [bytes]."
const INTERNAL_SHUNT_MATRIX_BYTES = 384 * 1024^2

const ShuntTapeFace = NamedTuple{
    (:kind,:radius,:phi,:span,:mid,:half,:alpha,:beta,:p,:nodes,:weights,
        :weighted,:beta_norm,:terminal,:points),
    Tuple{Symbol,Float64,Float64,Float64,Float64,Float64,Float64,Float64,
        Int,Vector{Float64},Vector{Float64},Matrix{Float64},Float64,Int,Vector{ComplexF64}}}

const InternalShuntDiagnostic = NamedTuple{
    (:boundary_residual,:wire_residual,:tape_residual,:common_residual,
        :penetration_indicator,:reciprocity,:log_moment_error,:unknowns,
        :equations,:matrix_bytes,:level),
    Tuple{Float64,Float64,Float64,Float64,Float64,Float64,Float64,
        Int,Int,Int,typeof(INTERNAL_SHUNT_RESOLUTION)}}

"""
$(TYPEDEF)

Store one prepared lossless local terminal operator. `C` is shield-referenced
capacitance \\[F/m\\]; `P` is its charge-potential inverse \\[m/F\\]. No
boundary matrices survive preparation.

$(TYPEDFIELDS)
"""
struct InternalShuntBlock{T <: Real}
    "Global assembly conductor range."
    assembly::UnitRange{Int}
    "Global local-domain terminals including the reference shield."
    terminals::UnitRange{Int}
    "Local capacitance \\[F/m\\]."
    C::Matrix{T}
    "Local charge-potential coefficients \\[m/F\\]."
    P::Matrix{T}
end

"""
$(TYPEDEF)

Own prepared internal shunt blocks and their bounded diagnostic records for
one calculation. Read-only blocks can be shared across identical lossless
formulations; no mutable numerical workspace is shared.

$(TYPEDFIELDS)
"""
struct PreparedInternalShunt{T <: Real, D}
    "Local blocks, in physical terminal order."
    blocks::Vector{InternalShuntBlock{T}}
    "Conductor-owned radial intervals replaced by the local blocks."
    covered::BitVector
    "Per-domain numerical diagnostics without matrices or geometry."
    diagnostics::D
    "Number of distinct boundary solves performed."
    solves::Int
end

function _shunt_values(domain::InternalShuntDomain{T}, methods, frequency, temperature) where {T}
    epsilon(material) = begin
        law = material.kind === :semicon ? methods.semicon_admittance : methods.insulation_admittance
        kappa = constitutive(law, material, frequency, temperature;
            temperature_dependence=methods.temperature_dependence)
        imag(kappa)/(2pi*frequency*(one(T)*8.8541878128e-12))
    end
    values = T[domain.a,domain.b,epsilon(domain.material)]
    for layers in (domain.left,domain.right), layer in layers
        append!(values,(layer.ri,layer.ro,epsilon(layer.material)))
    end
    for wire in domain.wires
        append!(values,(wire.x,wire.y,wire.r))
    end
    for tape in domain.tapes
        append!(values,(tape.ri,tape.ro,tape.phi,tape.span))
    end
    return values
end

# The flat values vector is also the narrow differentiation boundary. It keeps
# correlated scalar graphs out of dense factorizations without discarding them.
function _shunt_data(values::AbstractVector{T}, domain) where {T}
    a,b,epsilon = values[1:3]
    Layer = NamedTuple{(:ri,:ro,:epsilon),Tuple{T,T,T}}
    left,right = Layer[],Layer[]
    cursor = 4
    for (destination, layers) in ((left,domain.left),(right,domain.right))
        for _ in layers
            push!(destination,(ri=values[cursor],ro=values[cursor+1],epsilon=values[cursor+2]))
            cursor += 3
        end
    end
    wires,tapes = ShuntWire{T}[],ShuntTape{T}[]
    for wire in domain.wires
        push!(wires,(x=values[cursor],y=values[cursor+1],r=values[cursor+2],terminal=wire.terminal))
        cursor += 3
    end
    for tape in domain.tapes
        push!(tapes,(ri=values[cursor],ro=values[cursor+1],phi=values[cursor+2],
            span=values[cursor+3],terminal=tape.terminal))
        cursor += 4
    end
    Rleft = sum(l->log(l.ro/l.ri)/l.epsilon,left;init=zero(T))
    Rright = sum(l->log(l.ro/l.ri)/l.epsilon,right;init=zero(T))
    return (;a,b,epsilon,left,right,wires,tapes,Rleft,Rright,
        Rtotal=Rleft+log(b/a)/epsilon+Rright, ports=length(domain.terminals)-1)
end

function _shunt_discretization(g,level; validation=false)
    faces = ShuntTapeFace[]
    for tape in g.tapes
        append!(faces,_shunt_tape_faces(g,tape,level.order,level.quadrature))
    end
    count = validation ? 3level.wire : 2level.wire
    points = _shunt_points(g,count,0.75;shift=validation ? 0.37 : 0.0).targets
    terminals = Int[w.terminal for w in g.wires for _ in 1:count]
    nface = validation ? max(41,12(level.order+1)) : max(24,4(level.order+1))
    for face in faces
        targets = _shunt_face_targets([face],nface;validation)
        append!(points,targets)
        append!(terminals,fill(face.terminal,length(targets)))
    end
    return (;faces,points,terminals,wire_targets=count*length(g.wires))
end

function _shunt_rhs!(rhs, points, terminals, g)
    fill!(rhs,0)
    @inbounds for i in eachindex(points)
        rhs[i,1] = -_shunt_core_voltage(points[i],g)
        rhs[i,terminals[i]] = 1
    end
    return rhs
end

function _shunt_matrix!(matrix, points, sources, faces, g, k)
    count = length(sources)
    _shunt_kernel_matrix!(@view(matrix[:,1:count]),points,sources,g,k)
    tape = _shunt_tape_columns!(@view(matrix[:,count+1:end]),points,faces,g,k)
    return tape.log_moment_error
end

function _shunt_charge_map(g, sources, faces, level)
    count = length(sources)+sum(f->f.p+1,faces;init=0)
    charge = zeros(g.ports,count)
    @inbounds for (index,wire) in pairs(g.wires)
        for j in (index-1)*level.wire+1:index*level.wire
            charge[1,j] = -_shunt_core_voltage(sources[j],g)
            charge[wire.terminal,j] = 1
        end
    end
    offset = length(sources)
    for face in faces
        charge[face.terminal,offset+1] = 1
        if face.kind === :arc
            charge[1,offset+1] = -_shunt_core_voltage(first(face.points),g)
        else
            @views mul!(charge[1:1,offset+1:offset+face.p+1],
                transpose(_shunt_core_voltage.(face.points,Ref(g))),face.weighted,-1.0,0.0)
        end
        offset += face.p+1
    end
    return charge
end

"""
$(TYPEDSIGNATURES)

Compute the lossless terminal capacitance of exposed wires and finite tapes
inside a layered circular reference shield. Round-wire auxiliary sources and
corner-weighted integrated tape charges share a column-equilibrated least-
squares solve. Charge is extracted by exact total-charge and core-potential
moments, not by symmetrizing the answer.

# Arguments

- `g`: Local radii/positions \\[m\\] and relative dielectric permittivities.

# Keywords

- `level`: Fixed numerical resolution, defaulting to the accepted prototype.
- `retain=false`: Retain the factorization only for local differentiation.

# Returns

- Local capacitance `C` \\[F/m\\], sampled numerical diagnostics, and optional
  differentiation state. Sampled residuals are not certified error bounds.

# Notes

The layered Green function and whole-face charge approach has precedent in
[Bernal1997](@cite), but this solve uses oversampled boundary collocation and
pivoted QR, not that paper's Galerkin system. The circular kernel, wire-source
placement at `0.75` times each wire radius, Jacobi adaptation, resolution and
acceptance checks are implementation choices. The cited methods do not
establish an accuracy bound for this combined discretization.

# Errors

Memory-budget, rank and boundary-resolution failures are `ErrorException`, not
geometry `DomainError`; Monte Carlo must not condition its samples on them.
"""
Base.@constprop :aggressive function _shunt_capacitance(g;level=INTERNAL_SHUNT_RESOLUTION,retain=false)
    return _shunt_capacitance(g,Val(retain);level)
end

function _shunt_capacitance(g,::Val{retain};level=INTERNAL_SHUNT_RESOLUTION) where {retain}
    estimated_unknowns = length(g.wires)*level.wire + 4length(g.tapes)*(level.order+1)
    estimated_rows = 2length(g.wires)*level.wire +
        4length(g.tapes)*max(24,4(level.order+1))
    estimated_unknowns > 0 && estimated_rows >= estimated_unknowns ||
        error("internal shunt: insufficient boundary equations")
    estimated_unknowns <= div(INTERNAL_SHUNT_MATRIX_BYTES,8estimated_rows) || error(
        "internal shunt: requested resolution exceeds the dense storage budget")
    discretization = _shunt_discretization(g,level)
    (;faces,points,terminals) = discretization
    sources = _shunt_points(g,level.wire,0.75).sources
    n = length(sources)+sum(f->f.p+1,faces;init=0)
    m = length(points)
    n > 0 && m >= n || error("internal shunt: insufficient boundary equations")
    n <= div(INTERNAL_SHUNT_MATRIX_BYTES,8m) || error(
        "internal shunt: $m by $n boundary matrix exceeds the $(INTERNAL_SHUNT_MATRIX_BYTES÷1024^2) MiB storage budget")
    k = _shunt_kernel_coefficients(g,level.modes)
    matrix = Matrix{Float64}(undef,m,n)
    moment_error = _shunt_matrix!(matrix,points,sources,faces,g,k)
    rhs = _shunt_rhs!(zeros(m,g.ports),points,terminals,g)
    scales = Vector{Float64}(undef,n)
    @inbounds for j in 1:n
        scales[j] = norm(@view matrix[:,j])
        isfinite(scales[j]) && scales[j] > 0 || error("internal shunt: degenerate charge column")
        @views matrix[:,j] ./= scales[j]
    end
    # In-place pivoted QR avoids retaining K, scaled K, U and V simultaneously.
    factor = qr!(matrix,ColumnNorm())
    diagonal = [abs(factor.factors[i,i]) for i in 1:n]
    cutoff = max(m,n)*eps(Float64)*maximum(diagonal)
    count(>(cutoff),diagonal) == n || error("internal shunt: unresolved rank-deficient charge basis")
    scaled_coefficients = factor\rhs
    coefficients = scaled_coefficients./scales
    charge = _shunt_charge_map(g,sources,faces,level)
    C = (2pi*8.8541878128e-12).*(charge*coefficients)
    C[1,1] += 2pi*8.8541878128e-12/g.Rtotal
    all(isfinite,C) || error("internal shunt: nonfinite terminal capacitance")
    checks = _shunt_discretization(g,level;validation=true)
    scratch = Matrix{Float64}(undef,min(192,length(checks.points)),n)
    errors = zeros(size(scratch,1),g.ports)
    boundary = common = wire_residual = tape_residual = 0.0
    for start in 1:192:length(checks.points)
        stop = min(start+191,length(checks.points))
        rows = 1:(stop-start+1)
        block = @view scratch[rows,:]
        selected = @view checks.points[start:stop]
        moment_error = max(moment_error,_shunt_matrix!(block,selected,sources,faces,g,k))
        residual = @view errors[rows,:]
        _shunt_rhs!(residual,selected,@view(checks.terminals[start:stop]),g)
        mul!(residual,block,coefficients,1.0,-1.0)
        for i in rows
            value = maximum(abs,@view residual[i,:])
            boundary = max(boundary,value)
            common = max(common,abs(sum(@view residual[i,:])))
            if start+i-1 <= checks.wire_targets
                wire_residual = max(wire_residual,value)
            else
                tape_residual = max(tape_residual,value)
            end
        end
    end
    boundary <= 0.06 || error("internal shunt: sampled boundary residual $boundary exceeds 0.06 V per unit excitation")
    reciprocity = norm(C-transpose(C))/norm(C)
    reciprocity <= 0.01 || error("internal shunt: terminal reciprocity did not resolve ($reciprocity)")
    minimum(eigvals!(copy((C+transpose(C))/2))) > 0 || error("internal shunt: nonpassive local capacitance")
    leakage = abs(sum(@view C[1,:]))
    coupling = sum(abs,@view C[1,2:end])
    indicator = iszero(leakage) ? Inf : common*coupling/leakage
    diagnostic = (;boundary_residual=boundary,wire_residual,tape_residual,
        common_residual=common,penetration_indicator=indicator,reciprocity,
        log_moment_error=moment_error,unknowns=n,equations=m,
        matrix_bytes=8m*n,level)
    state = if retain
        residual = factor.Q' * rhs
        residual[1:n,:] .= 0
        residual = factor.Q * residual
        (;factor,scales,scaled_coefficients,coefficients,charge,residual,
            points,terminals,sources,faces,k,level)
    else
        nothing
    end
    return (;C,diagnostic,state)
end

"""
$(TYPEDSIGNATURES)

Differentiate the fixed-resolution local least-squares problem in one physical
parameter direction. Centered differences evaluate only kernel/charge-map
derivatives; the nominal pivoted QR is reused. Include the nonzero collocation
residual in the differentiated stationarity equation:

```math
A^T A\\,dx=A^T(db-dA\\,x)+dA^T(b-Ax).
```

Triangular QR solves evaluate this equation without forming normal equations.
Boundary derivatives are streamed in bounded row blocks, not stored as a dense
matrix of uncertain scalars. Input directions retain their original physical
units; the direction parameter and `step` are dimensionless.

# Returns

- Directional derivative of terminal capacitance \\[F/m\\].
"""
function _shunt_tangent(values, direction, domain, state, step)
    plus = _shunt_data(values .+ step.*direction,domain)
    minus = _shunt_data(values .- step.*direction,domain)
    dp,dm = _shunt_discretization(plus,state.level),_shunt_discretization(minus,state.level)
    sp = _shunt_points(plus,state.level.wire,0.75).sources
    sm = _shunt_points(minus,state.level.wire,0.75).sources
    kp,km = _shunt_kernel_coefficients(plus,state.level.modes),
        _shunt_kernel_coefficients(minus,state.level.modes)
    m,n = size(state.factor)
    count = min(192,m)
    matrix_p,matrix_m = Matrix{Float64}(undef,count,n),Matrix{Float64}(undef,count,n)
    bp,bm = zeros(count,plus.ports),zeros(count,plus.ports)
    drive,stationarity = zeros(m,plus.ports),zeros(n,plus.ports)
    for start in 1:192:m
        stop = min(start+191,m)
        rows = 1:(stop-start+1)
        ap,am = @view(matrix_p[rows,:]),@view(matrix_m[rows,:])
        xp,xm = @view(dp.points[start:stop]),@view(dm.points[start:stop])
        _shunt_matrix!(ap,xp,sp,dp.faces,plus,kp)
        _shunt_matrix!(am,xm,sm,dm.faces,minus,km)
        ap .-= am
        ap ./= 2step
        rp,rm = @view(bp[rows,:]),@view(bm[rows,:])
        _shunt_rhs!(rp,xp,@view(dp.terminals[start:stop]),plus)
        _shunt_rhs!(rm,xm,@view(dm.terminals[start:stop]),minus)
        target = @view drive[start:stop,:]
        target .= (rp .- rm)./(2step)
        mul!(target,ap,state.coefficients,-1.0,1.0)
        mul!(stationarity,transpose(ap),@view(state.residual[start:stop,:]),1.0,1.0)
    end
    stationarity ./= state.scales
    delta = state.factor\drive
    R = UpperTriangular(@view state.factor.factors[1:n,1:n])
    correction = R\(transpose(R)\stationarity[state.factor.p,:])
    @inbounds for i in 1:n, j in axes(delta,2)
        delta[state.factor.p[i],j] += correction[i,j]
    end
    delta ./= state.scales
    charge_delta = (_shunt_charge_map(plus,sp,dp.faces,state.level) -
        _shunt_charge_map(minus,sm,dm.faces,state.level))./(2step)
    result = (2pi*8.8541878128e-12).*(charge_delta*state.coefficients + state.charge*delta)
    result[1,1] += 2pi*8.8541878128e-12*(inv(plus.Rtotal)-inv(minus.Rtotal))/(2step)
    return result
end

"""
$(TYPEDSIGNATURES)

Evaluate a prepared local shunt domain from its physical scalar descriptors.
This non-exported extension protocol lets Measurements preserve correlations
without introducing uncertain scalars into dense numerical workspaces.

# Arguments

- `values`: Host/layer radii, wire coordinates/radii and tape dimensions
  \\[m\\], relative permittivities and tape angles \\[rad\\], in the order
  supplied by Engine preparation.
- `domain`: Engine-prepared geometry and terminal ownership.

# Keywords

- `level`: Numerical resolution of the accepted finite-strip calculation.
- `retain=false`: Retain the nominal factorization for differentiation tests.
- `directions=nothing`: Optional columns of physical input perturbations per
  unit dimensionless direction parameter.

# Returns

- Capacitance `C` \\[F/m\\], numerical diagnostics and optional retained state.
  With `directions`, also return checked implicit capacitance derivatives.
"""
function internal_shunt_response(values::AbstractVector{<:Real}, domain;
        level=INTERNAL_SHUNT_RESOLUTION,retain=false,directions=nothing)
    numerical = Float64.(values)
    result = try
        _shunt_capacitance(_shunt_data(numerical,domain);level,
            retain=retain || directions !== nothing)
    catch exception
        exception isa ErrorException || rethrow()
        error("internal shunt, design $(domain.design), terminals $(domain.terminals): " * exception.msg)
    end
    directions === nothing && return result
    size(directions,1) == length(numerical) || throw(DimensionMismatch(
        "internal shunt directions must align with physical descriptors"))
    derivatives = Matrix{Float64}(undef,size(directions,2),length(result.C))
    for column in axes(directions,2)
        direction = @view directions[:,column]
        relative = maximum(eachindex(direction)) do i
            abs(direction[i])/(iszero(numerical[i]) ? abs(numerical[2]) : abs(numerical[i]))
        end
        iszero(relative) && (derivatives[column,:] .= 0; continue)
        step = cbrt(eps(Float64))/relative
        coarse = _shunt_tangent(numerical,direction,domain,result.state,step)
        fine = _shunt_tangent(numerical,direction,domain,result.state,step/2)
        all(isfinite,fine) && norm(fine-coarse) <=
            0.02max(norm(fine),norm(coarse)) + 256eps(Float64)*norm(result.C)/step ||
            error("internal shunt, design $(domain.design), terminals $(domain.terminals): " *
                "uncertainty sensitivity did not resolve on step refinement")
        derivatives[column,:] .= vec(fine)
    end
    return (;C=result.C,diagnostic=result.diagnostic,
        state=retain ? result.state : nothing,tangents=derivatives)
end

function prepare_internal_shunt(domains::Vector{InternalShuntDomain{T}}, count,
        methods, frequency, temperature;level=INTERNAL_SHUNT_RESOLUTION) where {T}
    isempty(domains) && return nothing
    _shunt_lossless(methods) || return nothing
    blocks = InternalShuntBlock{T}[]
    covered = falses(count)
    # A concrete vector is inferred from the first result, independent of cable
    # count or shape. It never retains the expensive collocation matrix.
    first_result = internal_shunt_response(_shunt_values(first(domains),methods,frequency,temperature),
        first(domains);level)
    diagnostics = typeof(first_result.diagnostic)[]
    solves = 0
    for (index,domain) in pairs(domains)
        previous = findfirst(i->_shunt_domain_equal(domain,domains[i]),1:(index-1))
        if previous === nothing
            result = index == 1 ? first_result :
                internal_shunt_response(_shunt_values(domain,methods,frequency,temperature),domain;level)
            C = Matrix{T}(result.C)
            P = lu(C)\Matrix{T}(I,size(C,1),size(C,1))
            push!(diagnostics,result.diagnostic)
            solves += 1
        else
            C,P = blocks[previous].C,blocks[previous].P
            push!(diagnostics,diagnostics[previous])
        end
        push!(blocks,InternalShuntBlock(domain.assembly,domain.terminals,C,P))
        covered[first(domain.terminals):last(domain.terminals)-1] .= true
    end
    return PreparedInternalShunt(blocks,covered,diagnostics,solves)
end

_shunt_covered(::Nothing, index) = false
_shunt_covered(shunt::PreparedInternalShunt, index) = shunt.covered[index]
_shunt_potential!(destination, ::Nothing) = destination
_shunt_admittance!(destination, ::Nothing, s) = destination

function _shunt_potential!(destination, shunt::PreparedInternalShunt)
    for block in shunt.blocks
        inner,reference = first(block.terminals),last(block.terminals)
        # Charge on the inner anchor includes all conductors shielded inside it.
        @inbounds for j in first(block.assembly):reference-1, i in first(block.assembly):reference-1
            destination[i,j] += block.P[max(1,i-inner+1),max(1,j-inner+1)]
        end
    end
    return destination
end

function _shunt_admittance!(destination, shunt::PreparedInternalShunt, s)
    for block in shunt.blocks
        first_index,reference = first(block.terminals),last(block.terminals)
        @inbounds for j in axes(block.C,2), i in axes(block.C,1)
            value = s*block.C[i,j]
            row,column = first_index+i-1,first_index+j-1
            destination[row,column] += value
            destination[row,reference] -= value
            destination[reference,column] -= value
            destination[reference,reference] += value
        end
    end
    return destination
end

function _shunt_details(domains, prepared)
    resolved = prepared !== nothing
    return (treatment=resolved ? :resolved_local : :equivalent_coaxial,
        domains=[(design=d.design,terminals=collect(d.terminals)) for d in domains],
        solves=resolved ? prepared.solves : 0,
        diagnostics=resolved ? prepared.diagnostics : InternalShuntDiagnostic[],
        reason=resolved ? :qualified_lossless_domain :
            isempty(domains) ? :equivalent_geometry : :equivalent_material_law)
end
