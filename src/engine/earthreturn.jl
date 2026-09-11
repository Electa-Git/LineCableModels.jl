function earth!(
        destination::AbstractMatrix, bindings::NamedTuple{(
            :selection, :cases)},
        earth, jω, formula, Γ, workspace,
        thickness
)
    bindings.selection === formula ||
        throw(ArgumentError("workspace is bound to a different earth formula selection"))
    prepared=workspace===nothing ? () : map(workspace.systems) do system
        prepare_earth_system!(system, jω, Γ, workspace)
    end
    foreach(bindings.cases) do binding
        resources=workspace
        if system_earth(binding.selection)&&haskey(binding.declaration.options, :integration)
            index=findfirst(system->system.selection===binding.selection, workspace.systems)
            index===nothing && throw(ArgumentError("missing prepared full earth system"))
            resources=merge(workspace, (
                unified = workspace.systems[index].response, state = prepared[index]))
        end
        earth!(destination, binding, earth, jω,
            binding.selection, Γ, resources, thickness)
    end
    return destination
end

function earth!(
        destination::AbstractMatrix{Complex{T}}, binding::NamedTuple{(
            :selection, :declaration, :interactions, :reductions)},
        earth, jω, formula, Γ, workspace,
        thickness
) where {T <: Real}
    thickness = media(formula) === Val(:stratified) ? thickness : nothing
    for interaction in binding.interactions
        index, pair=interaction.index, interaction.pair
        rho=@view earth.rho[:, index]
        epsilon=@view earth.epsilon[:, index]
        mu=@view earth.mu[:, index]
        functor=if system_earth(formula)&&haskey(binding.declaration.options, :integration)
            prepared_earth_functor(formula, workspace.state, pair,
                binding.declaration, interaction.physical_pair)
        else
            formula(rho, epsilon, mu, jω, pair, binding.declaration; Γ,
                thickness, physical_pair = interaction.physical_pair)
        end
        resources = haskey(binding.declaration.options, :integration) ? workspace : nothing
        destination[pair.row, pair.column]=functor(resources)
    end
    return destination
end

@inline _gamma(::Nothing, frequency::Int) = nothing
@inline _gamma(values::AbstractVector, frequency::Int) = values[frequency]

# Numerical errors have the same units as the kernel value. They are distinct
# from the correlated physical uncertainty carried by Measurements values.
struct EarthKernelEstimate{V, R}
    value::V
    error::R
end
function EarthKernelEstimate(estimate::SpectralEstimate)
    EarthKernelEstimate(estimate.value, estimate.error)
end

function earth_kernel_estimate(estimate::SpectralEstimate, numerical)
    if numerical!==nothing&&haskey(numerical, :report)
        report=numerical.report
        report.integrals[]+=1
        report.evaluations[]+=estimate.evaluations
        report.samples[]+=estimate.samples
        report.cutoff[]=max(report.cutoff[], estimate.cutoff)
    end
    return EarthKernelEstimate(estimate)
end
Base.zero(a::EarthKernelEstimate) = EarthKernelEstimate(zero(a.value), zero(a.error))
function Base.:*(a::Number, b::EarthKernelEstimate)
    return EarthKernelEstimate(a*b.value, spectral_magnitude(a)*b.error)
end
Base.:*(a::EarthKernelEstimate, b::Number) = b*a
Base.:/(a::EarthKernelEstimate, b::Number) = inv(b)*a
Base.:-(a::EarthKernelEstimate) = EarthKernelEstimate(-a.value, a.error)
function Base.:+(a::EarthKernelEstimate, b::EarthKernelEstimate)
    EarthKernelEstimate(a.value+b.value, a.error+b.error)
end
Base.:-(a::EarthKernelEstimate, b::EarthKernelEstimate) = a+(-b)
function Base.:+(a::EarthKernelEstimate, b::Number)
    R=typeof(a.error)
    return EarthKernelEstimate(a.value+b, a.error+8eps(R)*spectral_magnitude(b))
end
Base.:+(a::Number, b::EarthKernelEstimate) = b+a
Base.:-(a::EarthKernelEstimate, b::Number) = a+(-b)

"""
$(TYPEDEF)

Own the complete exterior K/H/L matrices and factorization work arrays.
The kernel matrices use a common source-column scaling; final Ze/Pe/Ye are
in physical current and generalized-charge coordinates.

$(TYPEDFIELDS)
"""
struct EarthReturnWorkspace{T <: Real, G, R, Key, S}
    "Exterior circles; lengths \\[m\\]."
    geometry::G
    "Scaled negative axial-field kernel \\[Ω/m\\]."
    K::Matrix{Complex{T}}
    "Scaled line-voltage kernel \\[m/F\\]."
    H::Matrix{Complex{T}}
    "Scaled total-current map \\[dimensionless\\]."
    L::Matrix{Complex{T}}
    "Physical exterior series impedance \\[Ω/m\\]."
    Ze::Matrix{Complex{T}}
    "Physical exterior potential coefficient \\[m/F\\]."
    Pe::Matrix{Complex{T}}
    "Physical exterior admittance \\[S/m\\]."
    Ye::Matrix{Complex{T}}
    "Factorization buffer for transposed right solves."
    factor::Matrix{Complex{T}}
    "Transposed right-hand side buffer."
    rhs::Matrix{Complex{T}}
    "Estimated absolute errors of K and H, in their respective units."
    errors::NamedTuple{(:K, :H), Tuple{Matrix{R}, Matrix{R}}}
    "Per-interaction relative integration controls."
    tolerances::Matrix{R}
    "Estimated absolute final Ze/Pe/Ye errors, in their respective units."
    output_errors::NamedTuple{(:Ze, :Pe, :Ye), Tuple{Matrix{R}, Matrix{R}, Matrix{R}}}
    "Resolved physical state identifying a reusable complete response."
    cache_key::Key
    "Whether the complete response belongs to the current frequency evaluation."
    prepared::typeof(Ref(false))
    "Number of complete response preparations, for workspace diagnostics."
    assemblies::typeof(Ref(0))
    "Reusable state, numerical sampling and solve-sensitivity buffers."
    scratch::S
end

function EarthReturnWorkspace(geometry::EarthReturnGeometry{T}) where {T}
    n=length(geometry.radius)
    arrays=ntuple(_->zeros(Complex{T}, n, n), 8)
    R=typeof(float(nominal(one(T))))
    errors=(K = zeros(R, n, n), H = zeros(R, n, n))
    output_errors=(Ze = zeros(R, n, n), Pe = zeros(R, n, n), Ye = zeros(R, n, n))
    key=(s = zero(Complex{T}), Γ = zero(Complex{T}), sigma = (zero(T), zero(T)),
        epsilon = (zero(T), zero(T)), mu = (zero(T), zero(T)), gamma2 = (
            zero(Complex{T}), zero(Complex{T})))
    scratch=(
        current = (x = zeros(Complex{T}, n), scaling = zeros(T, n),
            A = zeros(Complex{T}, n), F = zeros(Complex{T}, n)),
        complex = ntuple(_->zeros(Complex{R}, n, n), 5), real = ntuple(_->zeros(R, n, n), 14),
        F = zeros(R, n), dirty = falses(n, n),
        cim_order = collect(CartesianIndices((n, n)))[:],
        report = (integrals = Ref(0), evaluations = Ref(0), samples = Ref(0),
            cutoff = Ref(zero(R)), refinements = Ref(0)),
        numerical = (
            segments = alloc_segbuf(R, Complex{T}, R; size = 128), images = Complex{T}[],
            seeds = alloc_segbuf(R, Complex{T}, R; size = 128), seed = alloc_segbuf(R, Complex{T}, R; size = 1),
            exponents = Complex{T}[], rules = (), spectral = spectral_scratch(R),
            cim = CIMWorkspace(),
            resolution = (phase = Ref(zero(R)), envelope = Ref(zero(R)), panels = Ref(0)),
            statistics = (evaluations = Ref(0), cutoff = Ref(zero(R)))))
    return EarthReturnWorkspace{T, typeof(geometry), R, typeof(Ref(key)), typeof(scratch)}(
        geometry, arrays..., errors,
        zeros(R, n, n), output_errors, Ref(key), Ref(false), Ref(0), scratch)
end

function earth_spectral_term(kind, pmedium, qmedium, state, hp, hq, y, radius, logscale,
        method, controls, numerical)
    R=typeof(float(nominal(real(state.s))))
    h=hp+hq
    angle=min(R(π)/6, atan(R(nominal(h))/(2max(R(nominal(y+radius)), eps(R)))))
    if method===Val(:cim)&&pmedium===Val(2)&&qmedium===Val(2)&&iszero(radius)
        cache=cim_workspace(numerical)
        cache===nothing || (angle=min(angle, R(cache.angle_limit[])))
    end
    angle=earth_contour_angle(state, angle)
    features=earth_spectral_features(state, h, y, radius, angle, numerical)
    scale=max(R(abs(nominal(state.k[1]))), R(abs(nominal(state.k[2]))), inv(R(nominal(h))))
    kernel=earth_spectrum(kind, pmedium, qmedium, state, hp, hq, logscale)
    integral=if kind===Val(:phi)&&method===Val(:cim)&&iszero(radius)&&!iszero(state.k[1]) &&
                abs(nominal(state.k[1]))<abs(nominal(state.k[2]))
        # The scalar kernel has a strong small-air-root feature. Its radial
        # coordinate must remove that root even when both wires are buried;
        # using ag leaves the feature compressed beside a much larger shift.
        radial=RadializedEarthSpectrum(kernel, state.k[1], h)
        SpectralIntegral(
            Val(:radial), radial, (height = h, separation = y, q = state.k[1]),
            scale; angle, features)
    elseif pmedium===Val(2)&&qmedium===Val(2)&&iszero(radius) &&
           kind in (Val(:Z), Val(:phi), Val(:voltage))&&method===Val(:cim)
        radial=EarthRadialSpectrum{
            typeof(kind).parameters[1], typeof(state), typeof(logscale)}(state, logscale)
        SpectralIntegral(
            Val(:radial), radial, (height = h, separation = y, q = state.k[2]),
            scale; angle, features)
    elseif iszero(radius)
        SpectralIntegral(
            Val(:cosine), kernel, (height = h, separation = y), scale; angle, features)
    else
        SpectralIntegral(Val(:besselcosine), kernel,
            (height = h, separation = y, radius), scale; angle, features)
    end
    return earth_kernel_estimate(spectral_estimate(method, integral, controls, numerical), numerical)
end

function unified_earth_pair(::Val{P}, ::Val{Q}, state, geometry, p, q, reference,
        method, controls, numerical) where {P, Q}
    u=state
    hp=abs(geometry.height[p])
    hq=abs(geometry.height[q])
    r=geometry.radius[p]
    y=abs(geometry.horizontal[p]-geometry.horizontal[q])
    zero_r=zero(r)
    sp=u.scaling[p]
    sq=u.scaling[q]
    πT=one(u.s)*π
    term(kind,
        hp,
        hq,
        radius,
        logscale) = earth_spectral_term(kind, Val(P), Val(Q), u,
        hp, hq, y, radius, logscale, method, controls, numerical)
    direct=P==Q ? earth_direct(u, geometry, p, q, P) : zero(u.s)
    z=u.s/πT*u.A[p]*term(Val(:Z), hp, hq, zero_r, sp+sq)
    P==Q && (z+=u.s*u.mu[P]/(2πT)*direct)
    phi=zero(z)
    if !iszero(u.Γ) || reference === :scalar
        phi=u.s/πT*u.A[p]*term(Val(:phi), hp, hq, zero_r, sp+sq)
        P==Q && (phi+=u.s/(2πT*u.sh[P])*direct)
    end
    k=z-u.Γ^2/u.s*phi
    if reference === :scalar
        return k, phi
    end
    air_limit=P==1&&abs(nominal(u.k[1]*hp))<0.01
    interface_path=reference===:interface&&abs(real(nominal(u.x[p])))<300
    combined_path=air_limit||interface_path
    h=if combined_path
        R=typeof(float(nominal(real(u.s))))
        padding=hq/2
        angle=min(R(π)/6, atan(R(nominal(hq))/(4max(R(nominal(y+r)), eps(R)))))
        angle=earth_contour_angle(u, angle)
        features=earth_spectral_features(u, hq-padding, y, r, angle, numerical)
        g=(hp, hq, radius = r, padding, logscale = sq, i0minus = bessel_i0m1(u.x[p]))
        kernel=reference===:interface ?
               EarthPathVoltageSpectrum{P, Q, :interface, typeof(u), typeof(g)}(u, g) :
               EarthPathVoltageSpectrum{P, Q, :deep, typeof(u), typeof(g)}(u, g)
        scale=max(R(abs(nominal(u.k[2]))), inv(R(nominal(hq))))
        integral=if method===Val(:cim)&&!iszero(u.k[P])
            height=hq-padding
            radial=RadializedEarthSpectrum(kernel, u.k[P], height)
            SpectralIntegral(
                Val(:radial), radial, (; height, separation = y, q = u.k[P]), scale;
                angle, features)
        else
            SpectralIntegral(
                Val(:cosine), kernel, (height = hq-padding, separation = y), scale;
                angle, features)
        end
        u.s/πT*earth_kernel_estimate(spectral_estimate(method, integral, controls, numerical), numerical)
    else
        u.s/πT*u.A[p]*term(Val(:voltage), hp, hq, zero_r, sp+sq)
    end
    P==Q && (h+=u.s/(2πT*u.sh[P])*direct)
    P==1&&!combined_path && (h+=u.s/πT*term(Val(:endpoint), zero_r, hq, r, sq))
    if reference === :interface&&!interface_path
        h-=u.s/πT*term(Val(:surface), zero_r, hq, r, sq)
    elseif reference isa Real
        if Q==1
            h-=u.s/πT*earth_spectral_term(Val(:surface), Val(2), Val(1), u,
                reference, hq, y, r, sq, method, controls, numerical)
        else
            for (kind, depth) in ((Val(:finite_direct), reference-hq), (
                Val(:finite_image), reference+hq))
                h-=u.s/πT*earth_spectral_term(kind, Val(2), Val(2), u,
                    depth, zero_r, y, r, sq, method, controls, numerical)
            end
        end
    end
    return k, h
end

# Solve X*A=B in reusable transposed buffers, preserving nonsymmetric entries.
function earth_right_solve!(destination, B, A, workspace)
    copyto!(workspace.factor, transpose(A))
    copyto!(workspace.rhs, transpose(B))
    factor=lu!(workspace.factor)
    ldiv!(factor, workspace.rhs)
    copyto!(destination, transpose(workspace.rhs))
    return destination
end

"""
$(TYPEDSIGNATURES)

Assemble the complete manuscript earth matrices with a common voltage reference:

```math
P_e L=H,\\qquad Z_e L=K+\\Gamma^2 H/s,\\qquad Y_e H=sL.
```

# Arguments

- `workspace`: Complete exterior geometry and reusable matrix storage.
- `state`: Evaluated medium data and prescribed Γ, prepared by the formula owner.
- `integration`: Normalized spectral method and controls.

# Keywords

- `reference`: `:deep` (default), `:interface`, `:scalar`, or a positive earth
  reference depth \\[m\\] below every circumference.
- `numerical`: Reusable integration workspace, or `nothing`.

# Returns

- `workspace`, with exterior Ze \\[Ω/m\\], Pe \\[m/F\\] and Ye \\[S/m\\].
  Internal impedance and insulation potential coefficients are not included.
"""
function unified_earth!(workspace::EarthReturnWorkspace, state, integration;
        reference = :deep, numerical = nothing)
    geometry=workspace.geometry
    workspace.assemblies[]+=1
    for counter in values(workspace.scratch.report)
        counter[]=zero(counter[])
    end
    numerical===nothing && (numerical=workspace.scratch.numerical)
    if integration.method===Val(:cim)
        cache=cim_workspace(numerical)
        if cache!==nothing
            for counter in values(cache.statistics)
                counter[]=0
            end
            buried=filter(i->geometry.height[i]<0, eachindex(geometry.height))
            if !isempty(buried)
                height=2minimum(i->abs(nominal(geometry.height[i])), buried)
                extent=maximum(nominal, geometry.horizontal)-minimum(nominal, geometry.horizontal)
                cache.angle_limit[]=min(pi/6, atan(height/(2max(extent, eps(Float64)))))
            end
        end
        # Difficult distant weights are prepared first, making their image
        # representations available to the less demanding nearby interactions.
        sort!(workspace.scratch.cim_order; by = index->(
            abs(nominal(geometry.height[index[1]]))+abs(nominal(geometry.height[index[2]])),
            -abs(nominal(geometry.horizontal[index[1]]-geometry.horizontal[index[2]]))))
    end
    numerical=merge(numerical, (report = workspace.scratch.report,))
    if isempty(numerical.rules)&&integration.method===Val(:trapz)
        T=eltype(geometry.radius)
        R=typeof(float(nominal(one(T))))
        controls=integration.options
        numerical=merge(numerical,
            (rules = ((controls, precision = precision(R),
                rule = spectral_de_rule(R, controls, nothing)),),))
    end
    if reference isa Real
        isfinite(reference)&&reference>maximum(-geometry.height .+ geometry.radius) ||
            throw(DomainError(reference, "finite reference must be below every circumference"))
    else
        reference in (:deep, :interface, :scalar) ||
            throw(ArgumentError("unknown earth voltage reference"))
    end
    u=unified_earth_state(state, geometry, workspace.scratch.current)
    controls=integration.options
    fill!(workspace.tolerances, controls.rtol)
    dirty=workspace.scratch.dirty
    fill!(dirty, true)
    for refinement in 0:7
        workspace.scratch.report.refinements[]=refinement
        for linear in eachindex(dirty)
            index=integration.method===Val(:cim) ? workspace.scratch.cim_order[linear] : CartesianIndices(dirty)[linear]
            p, q=Tuple(index)

            dirty[p, q] || continue
            P=geometry.height[p]>0 ? Val(1) : Val(2)
            Q=geometry.height[q]>0 ? Val(1) : Val(2)
            local_controls=merge(controls, (rtol = workspace.tolerances[p, q],))
            k,
            h=unified_earth_pair(P, Q, u, geometry, p, q, reference,
                integration.method, local_controls, numerical)
            workspace.K[p, q]=k.value
            workspace.H[p, q]=h.value
            workspace.errors.K[p, q]=k.error
            workspace.errors.H[p, q]=h.error
            workspace.L[p, q]=(p==q ? inv(u.A[p]) : zero(k.value))-u.F[p]*k.value
        end
        earth_right_solve!(workspace.Pe, workspace.H, workspace.L, workspace)
        @. workspace.Ze=workspace.K+u.Γ^2/u.s*workspace.H
        earth_right_solve!(workspace.Ze, workspace.Ze, workspace.L, workspace)
        @. workspace.Ye=u.s*workspace.L
        earth_right_solve!(workspace.Ye, workspace.Ye, workspace.H, workspace)
        earth_refine_errors!(dirty, workspace, u, controls) || return workspace
    end
    throw(ErrorException("earth matrix entries did not meet their propagated numerical error budgets"))
end

# Use the exact perturbed right-solve inequality, with estimated input errors.
# This is an estimated envelope, not a proof that backend estimates are bounds.
function earth_inverse_magnitude!(destination, A, scratch)
    factor, rhs=scratch.complex[4:5]
    @. factor=nominal(A)
    fill!(rhs, zero(eltype(rhs)))
    for i in axes(rhs, 1)
        rhs[i, i]=one(eltype(rhs))
    end
    ldiv!(lu!(factor), rhs)
    @. destination=abs(rhs)
    return destination
end

function earth_solve_error!(destination, X, A, B, EA, EB, C, scratch)
    R=eltype(EA)
    Xn, An, product=scratch.complex[1:3]
    absX, absA, tmp, residual, feedback, factor, rhs=scratch.real[1:7]
    @. Xn=nominal(X)
    @. An=nominal(A)
    @. absX=abs(Xn)
    @. absA=abs(An)
    mul!(product, Xn, An)
    mul!(tmp, absX, absA)
    @. residual=abs(product-nominal(B))+8eps(R)*(tmp+abs(nominal(B)))
    mul!(feedback, EA, C)
    if norm(feedback, Inf)>=1
        fill!(destination, R(Inf))
        return destination
    end
    mul!(tmp, absX, EA)
    @. residual+=EB+tmp
    mul!(destination, residual, C)
    # Right solve by I-EA*abs(inv(A)), in reusable transposed storage.
    copyto!(factor, transpose(feedback))
    factor .*= -one(R)
    for i in axes(factor, 1)
        factor[i, i]+=one(R)
    end
    copyto!(rhs, transpose(destination))
    ldiv!(lu!(factor), rhs)
    copyto!(destination, transpose(rhs))
    return destination
end

function earth_refine_errors!(dirty, workspace, state, controls)
    w=workspace
    R=eltype(w.tolerances)
    n=size(w.K, 1)
    scratch=w.scratch
    error_rhs, unused, EL, CL, CH, budgetK, budgetH=scratch.real[8:14]
    # These solves are for sensitivity estimates only. Physical matrices above
    # are obtained directly by right solves and never from a stored inverse.
    earth_inverse_magnitude!(CL, w.L, scratch)
    earth_inverse_magnitude!(CH, w.H, scratch)
    F=scratch.F
    @. F=abs(nominal(state.F))
    s=nominal(state.s)
    g=nominal(state.Γ^2/state.s)
    @. EL=F*w.errors.K+8eps(R)*abs(nominal(w.L))
    earth_solve_error!(w.output_errors.Pe, w.Pe, w.L, w.H, EL, w.errors.H, CL, scratch)
    @. w.rhs=w.K+g*w.H
    @. error_rhs=w.errors.K+abs(g)*w.errors.H
    earth_solve_error!(w.output_errors.Ze, w.Ze, w.L, w.rhs, EL, error_rhs, CL, scratch)
    @. w.rhs=s*w.L
    @. error_rhs=abs(s)*EL
    earth_solve_error!(
        w.output_errors.Ye, w.Ye, w.H, w.rhs, w.errors.H, error_rhs, CH, scratch)
    fill!(budgetK, R(Inf))
    fill!(budgetH, R(Inf))
    diagnostic=(kind = :Ze, row = 0, column = 0, estimate = zero(R), target = one(R))
    diagnostic=earth_entry_budgets!(
        Val(:Ze), 1e-14, w, state, controls, CL, CH, F, s, g, diagnostic)
    diagnostic=earth_entry_budgets!(
        Val(:Pe), 1e-14, w, state, controls, CL, CH, F, s, g, diagnostic)
    diagnostic=earth_entry_budgets!(
        Val(:Ye), 1e-10, w, state, controls, CL, CH, F, s, g, diagnostic)
    fill!(dirty, false)
    diagnostic.row==0 && return false
    for i in eachindex(dirty)
        ek=w.errors.K[i]
        eh=w.errors.H[i]
        ek<=budgetK[i]&&eh<=budgetH[i] && continue
        factor=min(R(0.5), budgetK[i]/max(2ek, floatmin(R)), budgetH[i]/max(2eh, floatmin(R)))
        next=max(8eps(R), w.tolerances[i]*factor)
        if next<w.tolerances[i]
            w.tolerances[i]=next
            dirty[i]=true
        end
    end
    any(dirty) || throw(ErrorException(
        "earth matrix error budget is limited by arithmetic precision or the absolute integral tolerance: $diagnostic"))
    return true
end

function earth_entry_budgets!(
        ::Val{Kind}, floor, w, state, controls, CL, CH, F, s, g, diagnostic) where {Kind}
    R=eltype(w.tolerances)
    n=size(w.K, 1)
    budgetK, budgetH=w.scratch.real[13:14]
    values=getproperty(w, Kind)
    errors=getproperty(w.output_errors, Kind)
    for q in 1:n, p in 1:n

        value=nominal(values[p, q])
        target=R(floor) + controls.rtol*min(abs(real(value)), abs(imag(value))) +
               32eps(R)*abs(value)
        errors[p, q]<=target && continue
        if errors[p, q]/target>diagnostic.estimate/diagnostic.target
            diagnostic=(kind = Kind, row = p, column = q, estimate = errors[p, q], target)
        end
        for b in 1:n, a in 1:n

            wk,
            wh=if Kind===:Pe
                (abs(nominal(w.Pe[p, a]*state.F[a]))*CL[b, q], p==a ? CL[b, q] : zero(R))
            elseif Kind===:Ze
                (abs((p==a ? one(R) : zero(R))+nominal(w.Ze[p, a]*state.F[a]))*CL[b, q],
                    p==a ? abs(g)*CL[b, q] : zero(R))
            else
                (p==a ? abs(s)*F[a]*CH[b, q] : zero(R), abs(nominal(w.Ye[p, a]))*CH[b, q])
            end
            wk>0 && (budgetK[a, b]=min(budgetK[a, b], target/(8n^2*wk)))
            wh>0 && (budgetH[a, b]=min(budgetH[a, b], target/(8n^2*wh)))
        end
    end
    return diagnostic
end

earth_state_equal(a, b) = isequal(a, b)
function earth_state_equal(a::Number, b::Number)
    a===b||(isequal(a, b)&&iszero(spectral_magnitude(a-b)))
end
function earth_state_equal(a::Tuple, b::Tuple)
    length(a)==length(b)&&all(pair->earth_state_equal(pair...), zip(a, b))
end
function earth_state_equal(a::NamedTuple, b::NamedTuple)
    keys(a)==keys(b)&&earth_state_equal(values(a), values(b))
end

function prepare_earth_system!(system, jω, Γ, numerical)
    data=system.materials
    for values in (data.rho, data.epsilon, data.mu)
        for column in axes(values, 2), row in axes(values, 1)

            earth_state_equal(values[row, column], values[row, 1]) || throw(ArgumentError(
                "the full-current default requires one globally consistent equivalent earth; pair-dependent reductions remain available through author formulas"))
        end
    end
    interaction=first(system.binding.interactions)
    declaration=first(system.declarations)
    functor=system.selection(
        @view(data.rho[:, 1]), @view(data.epsilon[:, 1]), @view(data.mu[:, 1]),
        jω, interaction.pair, declaration; Γ, physical_pair = interaction.physical_pair)
    if Γ===nothing
        for pair in system.binding.interactions
            value=declaration.hooks.Γ(jω, functor.state, pair.pair.layers)
            earth_state_equal(value, functor.state.Γ) || throw(ArgumentError(
                "the full-current default requires one common prescribed Γ for the complete system"))
        end
    end
    T=eltype(system.response.geometry.radius)
    state=functor.state
    key=(s = Complex{T}(state.jω), Γ = Complex{T}(state.Γ),
        sigma = ntuple(i->T(state.sigma[i]), 2), epsilon = ntuple(i->T(state.epsilon[i]), 2),
        mu = ntuple(i->T(state.mu[i]), 2), gamma2 = ntuple(i->Complex{T}(state.gamma_medium_squared[i]), 2))
    if !system.response.prepared[]||!earth_state_equal(system.response.cache_key[], key)
        unified_earth!(system.response, state, declaration.options.integration;
            reference = get(system.selection.parameters, :reference, :deep), numerical)
        system.response.cache_key[]=key
        system.response.prepared[]=true
    end
    return functor.state
end

function prepared_earth_functor(
        ::EarthImpedance.Formula{ID}, state, pair, selected, physical_pair) where {ID}
    binding=(; pair, physical_pair, kind = selected.kind, equation = selected.equation)
    return EarthImpedance.Functor{ID, typeof(binding), typeof(selected.hooks),
        typeof(state), typeof(selected.options)}(
        binding, selected.hooks, state, selected.options)
end
function prepared_earth_functor(
        ::EarthAdmittance.Formula{ID}, state, pair, selected, physical_pair) where {ID}
    binding=(; pair, physical_pair, kind = selected.kind, equation = selected.equation)
    return EarthAdmittance.Functor{ID, typeof(binding), typeof(selected.hooks),
        typeof(state), typeof(selected.options)}(
        binding, selected.hooks, state, selected.options)
end

function unified_entry(::Val{:impedance}, functor, pair, workspace)
    workspace!==nothing&&haskey(workspace, :unified) || throw(ArgumentError(
        "the default earth formula requires a prepared full-system context; use compute for physical matrices or select an author formula for an isolated pair"))
    return workspace.unified.Ze[pair.row, pair.column]
end
function unified_entry(::Val{:potential}, functor, pair, workspace)
    workspace!==nothing&&haskey(workspace, :unified) || throw(ArgumentError(
        "the default earth formula requires a prepared full-system context; use compute for physical matrices or select an author formula for an isolated pair"))
    return workspace.unified.Pe[pair.row, pair.column]
end
