function earth_spectral_term(
        kind::Val{Kind}, pmedium, qmedium, state, hp, hq, y, radius, logscale,
        method, controls, numerical; context = (;)) where {Kind}
    R=typeof(float(nominal(real(state.s))))
    h=hp+hq
    angle=min(R(π)/6, atan(R(nominal(h))/(2max(R(nominal(y+radius)), eps(R)))))
    angle=earth_contour_angle(state, angle)
    features=_unified_points!(numerical.unified, state, h, y, radius, angle)
    scale=max(R(abs(nominal(state.k[1]))), R(abs(nominal(state.k[2]))), inv(R(nominal(h))))
    kernel=earth_spectrum(kind, pmedium, qmedium, state, hp, hq, logscale)
    contour=scale*cis(angle)
    integral=SpectralIntegral(t->contour*earth_spectral_value(
        kernel, contour*t, h, y, radius))
    features ./= scale
    push!(features, one(scale))
    value,
    _ = integrate(method, integral, controls,
        numerical.quadrature;
        points = features, coordinate_type = R, context = merge(context, (term = Kind,)),
        observations = numerical.observations)
    return value
end

function unified_earth_pair(::Val{P}, ::Val{Q}, state, geometry, p, q,
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
    context=(formula = :unified, family = :earth_current, frequency = imag(u.s)/(2π),
        receiver = p, source = q)
    term(kind,
        hp,
        hq,
        radius,
        logscale) = earth_spectral_term(kind, Val(P), Val(Q), u,
        hp, hq, y, radius, logscale, method, controls, numerical; context)
    direct=P==Q ? earth_direct(u, geometry, p, q, P) : zero(u.s)
    z=u.s/πT*u.A[p]*term(Val(:Z), hp, hq, zero_r, sp+sq)
    P==Q && (z+=u.s*u.mu[P]/(2πT)*direct)
    phi=zero(z)
    if !iszero(u.Γ)
        phi=u.s/πT*u.A[p]*term(Val(:phi), hp, hq, zero_r, sp+sq)
        P==Q && (phi+=u.s/(2πT*u.sh[P])*direct)
    end
    k=z-u.Γ^2/u.s*phi
    combined_path=P==1&&abs(real(nominal(u.x[p])))<300
    h=if combined_path
        R=typeof(float(nominal(real(u.s))))
        padding=hq/2
        angle=min(R(π)/6, atan(R(nominal(hq))/(4max(R(nominal(y+r)), eps(R)))))
        angle=earth_contour_angle(u, angle)
        features=_unified_points!(numerical.unified, u, hq-padding, y, r, angle)
        g=(hp, hq, radius = r, padding, logscale = sq, i0minus = bessel_i0m1(u.x[p]))
        kernel=AirVoltageSpectrum{Q, typeof(u), typeof(g)}(u, g)
        scale=max(R(abs(nominal(u.k[2]))), inv(R(nominal(hq))))
        contour=scale*cis(angle)
        integral=SpectralIntegral(t->contour*earth_spectral_value(
            kernel, contour*t, hq-padding, y, zero_r))
        features ./= scale
        push!(features, one(scale))
        value,
        _ = integrate(method, integral, controls,
            numerical.quadrature;
            points = features, coordinate_type = R, context = merge(context, (term = :air_voltage,)),
            observations = numerical.observations)
        u.s/πT*value
    else
        u.s/πT*u.A[p]*term(Val(:voltage), hp, hq, zero_r, sp+sq)
    end
    P==Q && (h+=u.s/(2πT*u.sh[P])*direct)
    P==1&&!combined_path && (h+=u.s/πT*term(Val(:air_reference), zero_r, hq, r, sq))
    return k, h
end

"""
$(TYPEDSIGNATURES)

Assemble the physical exterior matrices by factorizing L once:

```math
P_e L=H,\\qquad Z_e L=K+\\Gamma^2 H/s.
```

Return `buffers` with Ze \\[Ω/m\\] and Pe \\[m/F\\]. Quadrature estimates
are reported by the numerical integrator; they do not reject a finite solve.
"""
function _unified_current!(numerical::NamedTuple, geometry, state, integration)
    buffers=numerical.unified
    u=_unified_state!(buffers.current, state, geometry)
    for q in axes(buffers.K, 2), p in axes(buffers.K, 1)

        P=geometry.height[p]>0 ? Val(1) : Val(2)
        Q=geometry.height[q]>0 ? Val(1) : Val(2)
        k,
        h=unified_earth_pair(P, Q, u, geometry, p, q,
            integration.method, integration.options, numerical)
        buffers.K[p, q]=k
        buffers.H[p, q]=h
        buffers.L[p, q]=(p==q ? inv(u.A[p]) : zero(k))-u.F[p]*k
    end
    copyto!(buffers.factor, transpose(buffers.L))
    factor=lu!(buffers.factor)
    copyto!(buffers.rhs, transpose(buffers.H))
    ldiv!(factor, buffers.rhs)
    copyto!(buffers.Pe, transpose(buffers.rhs))
    @. buffers.Ze=buffers.K+u.Γ^2/u.s*buffers.H
    copyto!(buffers.rhs, transpose(buffers.Ze))
    ldiv!(factor, buffers.rhs)
    copyto!(buffers.Ze, transpose(buffers.rhs))
    return buffers
end

function _unified_current!(workspace, materials, binding, frequency)
    selected=binding.selection
    s=workspace.input.jω[frequency]
    isfinite(s) && !iszero(s) || throw(DomainError(s, "jω must be finite and nonzero"))
    prescribed = first(binding.equations).declaration.options.data.Γ
    longitudinal = prescribed isa Number ? prescribed : prescribed[frequency]
    longitudinal isa Number && isfinite(longitudinal) ||
        throw(ArgumentError("Γ must be one finite scalar [1/m]"))
    for column in axes(materials.rho, 2)
        validate(selected, @view(materials.rho[:, column]),
            @view(materials.epsilon[:, column]), @view(materials.mu[:, column]), materials.thickness)
    end
    for values in (materials.rho, materials.epsilon, materials.mu)
        for column in axes(values, 2), row in axes(values, 1)

            same_physical_state(values[row, column], values[row, 1]) ||
                throw(ArgumentError("unified current closure requires globally consistent equivalent media"))
        end
    end
    state = (jω = s, Γ = oftype(s, longitudinal),
        sigma = ntuple(m -> conductivity(materials.rho[m, 1]), 2),
        epsilon = ntuple(m -> materials.epsilon[m, 1], 2),
        mu = ntuple(m -> materials.mu[m, 1], 2))
    return _unified_current!(workspace.buffers, workspace.invariants.geometry,
        state, first(binding.equations).declaration.options.data.integration)
end

# Full-current closure needs every exterior circumference, including those whose
# response entries are supplied by another selected formula.
function _unified_geometry(pairs::AbstractVector{<:EarthPair})
    radii = Dict(pair.row => something(pair.radius)
    for pair in pairs if pair.row == pair.column)
    for pair in pairs
        validate(pair)
        if pair.row == pair.column
            radii[pair.row] < abs(pair.heights[1]) || throw(DomainError(pair.radius,
                "each exterior circumference must lie wholly in one half-space"))
        else
            hypot(pair.separation, pair.heights[1] - pair.heights[2]) >
            radii[pair.row] + radii[pair.column] || throw(DomainError(
                (pair.row, pair.column), "exterior circumferences must not overlap"))
        end
    end
    return pairs
end
