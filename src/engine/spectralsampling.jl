"""
$(TYPEDEF)

Declare mandatory positive spectral breakpoints before numerical refinement.
Locations and widths come from the kernel's mathematical structure. The optional
`tail` callable returns an absolute remainder bound after a physical cutoff;
`nothing` retains estimated error semantics. A declared contour rotation must
be admissible for the complete kernel and analytic weight.

$(TYPEDFIELDS)
"""
struct SpectralFeatures{T, F}
    "Feature breakpoints in the real contour parameter \\[1/m\\]."
    points::Vector{T}
    "Absolute integrated tail envelope, or nothing."
    tail::F
end

function SpectralFeatures(points::AbstractVector{T}; tail = nothing) where {T <:
                                                                            AbstractFloat}
    all(x->isfinite(x)&&x>=0, points) ||
        throw(DomainError(points, "spectral feature points must be finite and nonnegative"))
    return SpectralFeatures{T, typeof(tail)}(sort!(unique!(collect(points))), tail)
end

"""
$(TYPEDEF)

Return a spectral value and numerical diagnostics in the integral's units.
`error` is an estimate unless all contributing errors have analytic bounds.

$(TYPEDFIELDS)
"""
struct SpectralEstimate{V, T}
    "Integral value."
    value::V
    "Estimated total absolute numerical error."
    error::T
    "Separate cutoff remainder; full-range integrations/certificates include tails in error."
    tail::T
    "Number of kernel evaluations."
    evaluations::Int
    "Largest quad/trapz evaluation, or CIM physical validation cutoff \\[1/m\\]."
    cutoff::T
end

# Supply the backend's error norm without changing physical values. In
# particular abs(Complex{Measurement}) is nondifferentiable at zero nominal
# magnitude; a numerical error norm must not introduce NaN derivatives there.
struct SpectralSample{V <: Number}<:Number
    value::V
    SpectralSample(value::V) where {V <: Union{Real, Complex}} = new{V}(value)
end
Base.zero(::Type{SpectralSample{V}}) where {V} = SpectralSample(zero(V))
Base.zero(x::SpectralSample) = SpectralSample(zero(x.value))
Base.:+(a::SpectralSample, b::SpectralSample) = SpectralSample(a.value+b.value)
Base.:-(a::SpectralSample, b::SpectralSample) = SpectralSample(a.value-b.value)
Base.:*(a::Real, b::SpectralSample) = SpectralSample(a*b.value)
Base.:*(a::SpectralSample, b::Real) = SpectralSample(a.value*b)
Base.:/(a::SpectralSample, b::Real) = SpectralSample(a.value/b)
norm(x::SpectralSample) = spectral_magnitude(x.value)

@inline spectral_rotation(integral) = integral.features === nothing ?
                                      one(integral.scale) : cis(integral.angle)

function spectral_scratch(::Type{R}) where {R}
    return (; features = sizehint!(R[], 128), seeds = sizehint!(R[], 128),
        scales = sizehint!(R[], 16), points = sizehint!(R[], 128), mapped = sizehint!(R[], 128))
end

function spectral_buffer(workspace, ::Type{R}, name) where {R}
    values=workspace===nothing||!haskey(workspace, :spectral) ? R[] :
           getproperty(workspace.spectral, name)
    empty!(values)
    return values
end

function spectral_sort_unique!(points)
    sort!(points; alg = Base.Sort.QuickSort)
    isempty(points) && return points
    count=1
    for index in 2:length(points)
        if !isequal(points[index], points[count])
            count+=1
            points[count]=points[index]
        end
    end
    return resize!(points, count)
end

function spectral_breakpoints(integral, workspace = nothing)
    R=typeof(integral.scale)
    points=spectral_buffer(workspace, R, :points)
    push!(points, zero(R), integral.scale)
    if integral.pole !== nothing
        push!(points, R(abs(integral.pole.location)))
    end
    if integral.features !== nothing
        append!(points, integral.features.points)
    end
    spectral_sort_unique!(points)
    return points
end

function spectral_estimate(::Val{:quad}, integral::SpectralIntegral, controls, workspace)
    R=typeof(integral.scale)
    scale=integral.scale
    rotation=spectral_rotation(integral)
    physical_points=spectral_breakpoints(integral, workspace)
    # Map explicitly: QuadGK's vector-domain API reuses segment coordinates,
    # so physical finite breakpoints must not accompany an implicit Inf map.
    points=spectral_buffer(workspace, R, :mapped)
    for x in physical_points
        push!(points, x/(scale+x))
    end
    push!(points, one(R))
    spectral_sort_unique!(points)
    segments=workspace===nothing ? nothing : workspace.segments
    statistics=workspace!==nothing&&haskey(workspace, :statistics) ? workspace.statistics :
               (evaluations = Ref(0), cutoff = Ref(zero(R)))
    evaluations=statistics.evaluations
    cutoff=statistics.cutoff
    evaluations[]=0
    cutoff[]=zero(R)
    f=t->begin
        evaluations[]+=1
        den=inv(one(R)-t)
        cutoff[]=max(cutoff[], t*den*scale)
        integral(rotation*(t*den*scale))*(rotation*scale*den^2)
    end
    value,
    error=if workspace!==nothing&&haskey(workspace, :seeds)
        # QuadGK's vector-domain convenience method allocates seed segments.
        # Populate reusable segments through its public evaluation API instead
        # of depending on its private Segment constructor or field layout.
        seeds=workspace.seeds
        empty!(seeds)
        sample=zero(eltype(workspace.images))
        for i in 1:(length(points) - 1)
            quadgk(_->sample, points[i], points[i + 1];
                segbuf = workspace.seed, norm = spectral_magnitude)
            append!(seeds, workspace.seed)
        end
        quadgk(f, zero(R), one(R); rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, segbuf = segments, eval_segbuf = seeds, norm = spectral_magnitude)
    else
        quadgk(f, points; rtol = controls.rtol, atol = controls.atol,
            maxevals = controls.maxevals, segbuf = segments, norm = spectral_magnitude)
    end
    isfinite(value)&&isfinite(error) &&
    error<=max(controls.atol, controls.rtol*spectral_magnitude(value)) ||
        throw(ErrorException("spectral :quad did not converge (estimated error=$error)"))
    return SpectralEstimate(value, R(error), zero(R), evaluations[], cutoff[])
end

# Reusable tables for the standard scalar types and default work limits.
# Only physical integrands carry uncertainty; sampling coordinates stay nominal.
const EARTH_DE_FLOAT64=QuadDE(Float64; maxlevel = 12)
const EARTH_DE_FLOAT32=QuadDE(Float32; maxlevel = 12)
function spectral_de_rule(::Type{R}, controls, workspace) where {R}
    if workspace!==nothing&&haskey(workspace, :rules)
        for table in workspace.rules
            table.controls.samples==controls.samples &&
                table.controls.max_refinements==controls.max_refinements &&
                table.precision==precision(R) && return table.rule
        end
    end
    controls.max_refinements==12&&controls.samples==128&&R===Float64 &&
        return EARTH_DE_FLOAT64
    controls.max_refinements==12&&controls.samples==128&&R===Float32 &&
        return EARTH_DE_FLOAT32
    return QuadDE(R; maxlevel = controls.max_refinements, h0 = min(one(R), R(128)/controls.samples))
end

function spectral_rule_bindings(bindings, ::Type{R}) where {R}
    controls=unique([case.declaration.options.integration.options
                     for case in bindings.cases
                     if haskey(case.declaration.options,
        :integration) &&
        case.declaration.options.integration.method===Val(:trapz)])
    return Tuple(map(controls) do options
        (controls = options, precision = precision(R),
            rule = spectral_de_rule(R, options, nothing))
    end)
end

# Geometry changes guide initial panel resolution; physical features and residual
# checks remain necessary because phase and decay do not locate narrow kernels.
function spectral_resolution(integral::SpectralIntegral{Kind}, left, right) where {Kind}
    w=integral.weight
    rotation=spectral_rotation(integral)
    delta=rotation*(right-left)
    root_delta=if Kind===:radial
        q=complex(nominal(real(w.q)), nominal(imag(w.q)))
        sqrt((rotation*right)^2+q^2)-sqrt((rotation*left)^2+q^2)
    else
        delta
    end
    height=nominal(w.height)
    extent=nominal(w.separation+get(w, :radius, zero(w.height)))
    return (phase = abs(height*imag(root_delta))+extent*abs(real(delta)),
        envelope = abs(height*real(root_delta))+extent*abs(imag(delta)))
end

function spectral_trapz_points(integral, controls, workspace)
    R=typeof(integral.scale)
    points=spectral_breakpoints(integral, workspace)
    rotation=spectral_rotation(integral)
    w=integral.weight
    decay=nominal(w.height)*real(rotation)-nominal(w.separation+get(w, :radius, 0))*abs(imag(rotation))
    # This is a finite seeding horizon only: the remaining interval to infinity
    # is always integrated and verified, never dropped on this heuristic.
    horizon=decay>0 ? R(max(20, -log(max(controls.rtol, eps(R))))/decay) : integral.scale
    isfinite(horizon)&&horizon>0 || (horizon=integral.scale)
    # Existing material features often already extend beyond this horizon;
    # do not introduce a redundant tail panel in that case.
    last(points)<horizon && push!(points, horizon)
    spectral_sort_unique!(points)
    seeds=spectral_buffer(workspace, R, :seeds)
    append!(seeds, points)
    maxphase=zero(R); maxenvelope=zero(R)
    for i in 1:(length(seeds)-1)
        left, right=seeds[i], min(seeds[i+1], horizon)
        left<right || continue
        demand=spectral_resolution(integral, left, right)
        demand_ratio=max(demand.phase/max(R(pi), controls.samples*R(pi)/40),
            demand.envelope/max(one(R), R(controls.samples)/32))
        count=ceil(Int, clamp(demand_ratio, one(R), R(controls.max_tail_refinements)))
        for j in 1:(count-1)
            push!(points, left+(right-left)*j/count)
        end
        maxphase=max(maxphase, R(demand.phase/count))
        maxenvelope=max(maxenvelope, R(demand.envelope/count))
    end
    spectral_sort_unique!(points)
    if workspace!==nothing&&haskey(workspace, :resolution)
        report=workspace.resolution
        report.phase[]=maxphase
        report.envelope[]=maxenvelope
        report.panels[]=length(points)
    end
    mapped=spectral_buffer(workspace, R, :mapped)
    for x in points
        push!(mapped, x/(integral.scale+x))
    end
    push!(mapped, one(R))
    spectral_sort_unique!(mapped)
    return mapped
end

function spectral_estimate(::Val{:trapz}, integral::SpectralIntegral, controls, workspace)
    R=typeof(integral.scale)
    points=spectral_trapz_points(integral, controls, workspace)
    rotation=spectral_rotation(integral)
    evaluations=Ref(0)
    cutoff=Ref(zero(R))
    f=t->begin
        evaluations[]+=1
        coordinate=min(t, prevfloat(one(R)))
        den=inv(one(R)-coordinate)
        parameter=integral.scale*coordinate*den
        cutoff[]=max(cutoff[], parameter)
        integral(rotation*parameter)*(rotation*integral.scale*den^2)
    end
    rule=spectral_de_rule(R, controls, workspace)
    segments=workspace===nothing ? nothing : get(workspace, :segments, nothing)
    return spectral_de_estimate(rule, f, points, controls, evaluations, cutoff, segments)
end

function de_reference_panel(f::F, a, b, controls, count, segments; target = nothing) where {F}
    R=typeof(a)
    value, error=quadgk(f, a, b;
        rtol = target===nothing ? max(8eps(R), controls.rtol/32) : zero(R),
        atol = target===nothing ? controls.atol/(32count) : target/(32count),
        maxevals = get(controls, :maxevals, 10^7), segbuf = segments, norm = spectral_magnitude)
    return (left = a, right = b, value, error)
end

# Independent panel references serve both normalization and verification. They
# replace the redundant full-domain reference plus a second set of panel solves.
function spectral_de_estimate(rule, f::F, points, controls, evaluations, cutoff, segments) where {F}
    R=eltype(points)
    n=length(points)-1
    references=[de_reference_panel(f, points[i], points[i+1], controls, n, segments) for i in 1:n]
    reference=sum(panel->panel.value, references)
    target=max(R(controls.atol), R(controls.rtol)*R(spectral_magnitude(reference)))
    if sum(panel->panel.error, references)>target/16
        for i in eachindex(references)
            panel=references[i]
            if panel.error>target/(16n)
                references[i]=de_reference_panel(f, panel.left, panel.right, controls, n, segments; target)
            end
        end
    end
    amplitude=max(sum(panel->R(spectral_magnitude(panel.value)), references),
        R(controls.atol), floatmin(R))
    panels=references
    for i in eachindex(panels)
        panels[i]=verified_de_panel(rule, f, panels[i], amplitude, target, controls, n, 0)
    end
    for refinement in 0:(controls.max_tail_refinements-1)
        value=sum(panel->panel.value, panels)
        rounding=8eps(R)*sum(panel->R(spectral_magnitude(panel.value)), panels)
        error=sum(panel->panel.error, panels)+rounding
        target=max(R(controls.atol), R(controls.rtol)*R(spectral_magnitude(value)))
        if isfinite(value)&&isfinite(error)&&error<=target
            return SpectralEstimate(value, error, zero(R), evaluations[], cutoff[])
        end
        index=argmax(map(panel->panel.error, panels))
        panel=panels[index]
        mid=(panel.left+panel.right)/2
        panel.left<mid<panel.right || break
        count=length(panels)+1
        left=de_reference_panel(f, panel.left, mid, controls, count, segments; target)
        right=de_reference_panel(f, mid, panel.right, controls, count, segments; target)
        panels[index]=verified_de_panel(rule, f, left, amplitude, target,
            controls, count, refinement+1)
        push!(panels, verified_de_panel(rule, f, right, amplitude, target,
            controls, count, refinement+1))
    end
    throw(ErrorException("spectral :trapz failed local error verification; no quadrature value was substituted"))
end

function verified_de_panel(rule, f, reference, amplitude, target, controls, count, refinement)
    R=typeof(reference.left)
    tolerance=target/amplitude/(32count)
    value, estimate=rule(t->SpectralSample(f(t)/amplitude), reference.left, reference.right;
        rtol = zero(R), atol = max(eps(R), 2tolerance/4^refinement))
    physical=value.value*amplitude
    error=max(R(nominal(estimate))*amplitude,
        R(spectral_magnitude(physical-reference.value))+reference.error)
    return (left = reference.left, right = reference.right, value = physical, error)
end

function cim_windows(integral::SpectralIntegral{Kind}, limit) where {Kind}
    R=typeof(integral.scale)
    points=spectral_breakpoints(integral)
    if Kind===:radial
        q=integral.weight.q
        points=R[abs(x^2/(sqrt(x^2+q^2)+q)) for x in points]
    end
    filter!(x->zero(R)<x<=limit, points)
    isempty(points)&&push!(points, min(integral.scale, limit))
    push!(points, limit)
    sort!(unique!(points))
    seeds=copy(points)
    for i in 1:(length(seeds) - 1)
        x=4seeds[i]
        while x<seeds[i + 1]
            push!(points, x)
            x*=4
        end
    end
    return sort!(unique!(points))
end

function cim_sample_weights(integral, nodes, rotation, shift)
    w=integral.weight
    radius=get(w, :radius, zero(w.height))
    rate=real(w.height*rotation)-abs(imag((w.separation+radius)*rotation))
    weights=similar(nodes)
    for i in eachindex(nodes)
        left=i==1 ? nodes[i] : nodes[i - 1]
        right=i==length(nodes) ? nodes[i] : nodes[i + 1]
        weights[i]=sqrt((right-left)/2)*exp(-rate*nodes[i])
    end
    return weights
end

function cim_fit_samples(integral::SpectralIntegral{Kind}, nodes, rotation, shift) where {Kind}
    if Kind !== :radial
        return nodes, integral.kernel.(rotation .* nodes .+ shift),
        cim_sample_weights(integral, nodes, rotation, shift)
    end
    # The radial fitting variable follows a curved path on the physical
    # integration contour. Uniform pencils supply exponents, but amplitudes
    # must be checked and fitted on that path, not extrapolated from real u-q.
    w=integral.weight
    λ=spectral_rotation(integral) .* nodes
    u=sqrt.(λ .^ 2 .+ w.q^2)
    coordinates=λ .^ 2 ./ (u .+ w.q)
    values=[radial_kernel_value(integral.kernel, u[i], λ[i]) for i in eachindex(u)]
    weights=similar(nodes)
    for i in eachindex(nodes)
        left=i==1 ? nodes[i] : nodes[i - 1]
        right=i==length(nodes) ? nodes[i] : nodes[i + 1]
        envelope=exp(-real(w.height*u[i])+abs(imag(w.separation*λ[i])))/abs(u[i])
        weights[i]=sqrt((right-left)/2)*envelope
    end
    return coordinates, values, weights
end

function cim_tail_penalty(integral::SpectralIntegral{Kind}, poles, rotation, limit) where {Kind}
    R=typeof(integral.scale)
    Kind===:radial && return zeros(R, length(poles))
    w=integral.weight
    rate=real(w.height*rotation)-abs(imag((w.separation+get(w, :radius, zero(R)))*rotation))
    return [exp(-(rate+real(b))*limit)/sqrt(2(rate+real(b))) for b in poles]
end

function cim_residual_error(integral::SpectralIntegral{Kind}, amplitudes, poles,
        rotation, shift, limit, target, controls) where {Kind}
    R=typeof(integral.scale)
    true_residual=SpectralIntegral(
        Val(Kind), integral.kernel, integral.weight, integral.scale;
        angle = integral.angle, features = integral.features)
    scale=integral.scale
    physical=sort!(unique!(vcat(spectral_breakpoints(integral), limit)))
    points=sort!(unique!(vcat(physical ./ (scale .+ physical), one(R))))
    f=t->begin
        den=inv(one(R)-t)
        λ=(Kind===:radial ? spectral_rotation(integral) : rotation)*scale*t*den
        fit=sum(i->cim_weighted_image(integral, amplitudes[i], poles[i], λ, rotation), eachindex(poles))
        abs(true_residual(λ)-fit)*scale*den^2
    end
    value,
    error=quadgk(f, points; rtol = 1e-3, atol = target/32, maxevals = controls.maxevals)
    return value+error
end

function cim_weighted_image(
        integral::SpectralIntegral{Kind}, a, b, λ, rotation) where {Kind}
    w=integral.weight
    if Kind===:radial
        u=sqrt(λ^2+w.q^2)
        coordinate=λ^2/(u+w.q)
        return a/u*(exp(-b*coordinate-w.height*u+im*w.separation*λ) +
                    exp(-b*coordinate-w.height*u-im*w.separation*λ))/2
    end
    growth=Kind===:besselcosine ? abs(imag(w.radius*λ)) : zero(real(λ))
    value=a*(exp(-b*λ/rotation-w.height*λ+im*w.separation*λ+growth) +
             exp(-b*λ/rotation-w.height*λ-im*w.separation*λ+growth))/2
    iszero(spectral_magnitude(value)) && return value
    return Kind===:besselcosine ? value*special_besseljx(0, w.radius*λ) : value
end

function cim_true_tail(integral, limit, rotation, target, controls)
    R=typeof(integral.scale)
    physical=filter(>(limit), spectral_breakpoints(integral))
    points=sort!(unique!(vcat(zero(R), (physical .- limit) ./ physical, one(R))))
    f=t->begin
        den=inv(one(R)-t)
        abs(integral(rotation*limit*den))*limit*den^2
    end
    value,
    error=quadgk(f, points; rtol = 1e-3, atol = target/64, maxevals = controls.maxevals)
    if integral.features!==nothing&&integral.features.tail!==nothing
        bound=integral.features.tail(limit)
        bound isa Real&&isfinite(bound)&&bound>=0 ||
            throw(DomainError(bound, "the spectral tail envelope must be finite and nonnegative"))
        value-error<=bound+target/64 ||
            throw(ErrorException("the declared spectral tail envelope is contradicted by independent integration"))
        return R(bound)
    end
    return value+error
end
