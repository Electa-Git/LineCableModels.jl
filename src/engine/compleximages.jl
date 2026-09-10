"""
$(TYPEDEF)

Identify an immutable earth spectral response for complex-image reuse. Physical
state and kernel geometry are copied as values; mutable current-map arrays are
excluded. Arbitrary user callables are not cached without an identity method.

$(TYPEDFIELDS)
"""
struct CIMIdentity
    "Spectral kernel family and radial transformation."
    family::Tuple{Symbol, Symbol, Bool}
    "Ordered receiver and source media."
    media::Tuple{Int, Int}
    "Frequency, propagation constant, material data and medium roots in SI units."
    state::NTuple{10, ComplexF64}
    "Parameters occurring inside the fitted kernel, in SI units."
    geometry::NTuple{8, ComplexF64}
end

"""
$(TYPEDEF)

Record an absolute numerical error certificate for a complex-image expansion.
An envelope certificate applies to greater heights and smaller separation plus
radius on the same contour; a point certificate applies only to its geometry.
Errors remain estimates because their residual integrals use numerical quadrature.

$(TYPEDFIELDS)
"""
struct CIMCertificate
    "Weight height \\[m\\]."
    height::Float64
    "Horizontal separation \\[m\\]."
    separation::Float64
    "Bessel weight radius \\[m\\]."
    radius::Float64
    "Physical contour angle \\[rad\\]."
    angle::Float64
    "Absolute integration error in the integral's units."
    error::Float64
    "Integrated tail estimate in the integral's units."
    tail::Float64
    "Physical validation cutoff \\[1/m\\]."
    cutoff::Float64
    "Whether the certificate integrates the geometry envelope."
    envelope::Bool
end

"""
$(TYPEDEF)

Own one reusable exponential expansion and its geometry certificates. Radial
exponents multiply u−q; other exponents multiply λ divided by `rotation`.

$(TYPEDFIELDS)
"""
struct CIMImageFit
    "Exact built-in kernel identity."
    identity::CIMIdentity
    "Analytic weight family."
    kind::Symbol
    "Physical scalar precision, before Float64 fitting."
    precision::Int
    "Source/receiver logarithmic normalization \\[dimensionless\\]."
    logscale::Float64
    "Exponential amplitudes in the residual kernel's units."
    images::Vector{ComplexF64}
    "Exponential decay lengths \\[m\\]."
    exponents::Vector{ComplexF64}
    "Affine fitting direction \\[dimensionless\\]."
    rotation::ComplexF64
    "Radial fitting origin q \\[1/m\\], or zero for a cosine fit."
    shift::ComplexF64
    "Verified point geometries or monotone geometry envelopes."
    certificates::Vector{CIMCertificate}
end

"""
$(TYPEDEF)

Own a bounded collection of complex-image fits, factorization buffers and work
diagnostics. The cache is local to a computation workspace; no physical values
or uncertainty correlations are identified by a hash or nominal approximation.

$(TYPEDFIELDS)
"""
struct CIMWorkspace
    "At most 64 recently constructed image expansions."
    fits::Vector{CIMImageFit}
    "Reusable matrix-pencil sample buffer."
    samples::Vector{ComplexF64}
    "Reusable Hankel matrix storage."
    hankel::Vector{ComplexF64}
    "Common admissible angle cap for same-earth weights \\[rad\\]."
    angle_limit::Base.RefValue{Float64}
    "Construction, reuse and certification counters."
    statistics::NamedTuple{(:fits, :pencils, :hits, :certifications), NTuple{4, Base.RefValue{Int}}}
end
CIMWorkspace() = CIMWorkspace(CIMImageFit[], ComplexF64[], ComplexF64[], Ref(pi/2),
    (fits = Ref(0), pencils = Ref(0), hits = Ref(0), certifications = Ref(0)))
cim_workspace(workspace) = workspace===nothing ? nothing : get(workspace, :cim, nothing)

"""
$(TYPEDEF)

Retain one sampled window's projected matrix pencil. Candidate image orders
reuse leading blocks of this projection and its singular values; changing
order requires no new kernel samples or Hankel factorization.

$(TYPEDFIELDS)
"""
struct CIMPencil
    "Projected shifted Hankel matrix, in the sampled kernel's units."
    shifted::Matrix{ComplexF64}
    "Retained Hankel singular values, in the sampled kernel's units."
    singular::Vector{Float64}
    "Uniform fitting-coordinate increment \\[1/m\\]."
    step::Float64
    "Right fitting-window endpoint \\[1/m\\]."
    width::Float64
end

function cim_pencil(data, step, width, tolerance, maxrank, cache)
    samples=length(data)
    rows=samples÷2
    columns=samples-rows
    storage=cache===nothing ? Vector{ComplexF64}(undef,2rows*columns) :
            resize!(cache.hankel,2rows*columns)
    H0=reshape(view(storage,1:(rows*columns)),rows,columns)
    H1=reshape(view(storage,(rows*columns+1):(2rows*columns)),rows,columns)
    for j in 1:columns, i in 1:rows
        H0[i,j]=data[i+j-1]
        H1[i,j]=data[i+j]
    end
    cache===nothing || (cache.statistics.pencils[]+=1)
    factor=svd!(H0)
    rank=min(maxrank,count(>(max(eps(Float64)*samples,tolerance*1e-3)*first(factor.S)),factor.S))
    U=@view factor.U[:,1:rank]
    V=@view factor.V[:,1:rank]
    return CIMPencil(U'*H1*V, factor.S[1:rank],step,width)
end

function cim_poles!(poles, pencils, terms, integral::SpectralIntegral{Kind}, rotation) where {Kind}
    w=integral.weight
    rate=real(w.height*rotation)-abs(imag((w.separation+get(w,:radius,0))*rotation))
    for pencil in pencils
        rank=min(terms,length(pencil.singular))
        rank==0 && continue
        reduced=pencil.shifted[1:rank,1:rank]
        for j in 1:rank, i in 1:rank
            reduced[i,j]/=pencil.singular[j]
        end
        for z in eigvals!(reduced)
            iszero(z) && continue
            b=-log(z)/pencil.step
            abs(b)*pencil.width<100eps(Float64) && (b=zero(b))
            physical_rate=Kind===:radial ? real((w.height+b)*spectral_rotation(integral))-
                abs(imag(w.separation*spectral_rotation(integral))) : real(b)+rate
            isfinite(b) && real(b)+rate>sqrt(eps(Float64))*rate && physical_rate>0 && push!(poles,b)
        end
    end
    return poles
end

cim_identity(kernel) = nothing
cim_identity(kernel::SpectralKernelCounter) = cim_identity(kernel.kernel)
cim_logscale(kernel) = 0.0
cim_logscale(kernel::SpectralKernelCounter) = cim_logscale(kernel.kernel)
cim_logscale(kernel::EarthRadialSpectrum) = Float64(kernel.logscale)

function cim_state_identity(u)
    return ComplexF64.((get(u, :s, 0), get(u, :Γ, 0), u.sh..., u.mu..., u.k2..., u.k...))
end
function cim_identity(kernel::EarthRadialSpectrum{Kind}) where {Kind}
    return CIMIdentity((:earth_radial, Kind, false), (2, 2),
        cim_state_identity(kernel.state), ntuple(_->zero(ComplexF64), 8))
end
function cim_identity(kernel::EarthSpectrum{Kind, P, Q}) where {Kind, P, Q}
    g=kernel.geometry
    return CIMIdentity((:earth, Kind, false), (P, Q), cim_state_identity(kernel.state),
        ComplexF64.((g.hp, g.hq, g.logscale, 0, 0, 0, 0, 0)))
end
function cim_identity(kernel::EarthPathVoltageSpectrum{P, Q, Reference}) where {P, Q, Reference}
    g=kernel.geometry
    return CIMIdentity((:path, Reference, false), (P, Q), cim_state_identity(kernel.state),
        ComplexF64.((g.hp, g.hq, g.radius, g.padding, g.logscale, g.i0minus, 0, 0)))
end
function cim_identity(kernel::RadializedEarthSpectrum)
    key=cim_identity(kernel.kernel)
    key===nothing && return nothing
    return CIMIdentity((key.family[1], key.family[2], true), key.media, key.state,
        (key.geometry[1:6]..., ComplexF64(kernel.q), ComplexF64(kernel.height)))
end

function cim_certificate(integral, error, tail, cutoff; envelope = false)
    w=integral.weight
    return CIMCertificate(Float64(w.height), Float64(w.separation),
        Float64(get(w, :radius, 0)), Float64(angle(spectral_rotation(integral))),
        Float64(error), Float64(tail), Float64(cutoff), envelope)
end
function cim_contains(certificate, integral)
    w=integral.weight
    certificate.angle==angle(spectral_rotation(integral)) || return false
    radius=get(w, :radius, 0)
    return certificate.envelope ?
           w.height>=certificate.height && w.separation+radius<=certificate.separation+certificate.radius :
           w.height==certificate.height && w.separation==certificate.separation && radius==certificate.radius
end

function cim_image_value(::Val{Kind}, w, images, poles, rotation) where {Kind}
    total=zero(ComplexF64)
    magnitude=0.0
    for i in eachindex(images, poles)
        a, b=images[i], poles[i]
        term=if Kind===:cosine
            p=w.height*rotation+b
            rotation*a*p/(p^2+(w.separation*rotation)^2)
        elseif Kind===:besselcosine
            p=w.height*rotation+b
            rotation*a*(inv(sqrt((p+im*w.separation*rotation)^2+(w.radius*rotation)^2))+
                        inv(sqrt((p-im*w.separation*rotation)^2+(w.radius*rotation)^2)))/2
        else
            d=sqrt((w.height+b)^2+w.separation^2)
            a*image_besselk0x(w.q*d)*exp(-w.q*(w.height+w.separation^2/(d+w.height+b)))
        end
        total+=term
        magnitude+=abs(term)
    end
    return total, 8eps(Float64)*magnitude
end

function cim_store_fit!(integral::SpectralIntegral{Kind}, controls, workspace,
        images, poles, rotation, error, tail, cutoff, ::Type{Result}) where {Kind, Result}
    cache=cim_workspace(workspace)
    cache===nothing && return
    key=cim_identity(integral.kernel)
    key===nothing && return
    integral.pole===nothing || return
    spectral_cacheable_tail(spectral_tail(integral)) || return
    length(cache.fits)>=64 && popfirst!(cache.fits)
    certificate=cim_certificate(integral, error, tail, cutoff)
    certificates=[certificate]
    target=max(controls.atol, controls.rtol*abs(first(cim_image_value(
        Val(Kind), integral.weight, images, poles, rotation))))
    cache.statistics.certifications[]+=1
    envelope=cim_envelope_error(integral, images, poles, rotation, 1.0, cutoff, target, controls)
    if isfinite(envelope.error)
        push!(certificates, cim_certificate(integral, envelope.error, envelope.tail,
            envelope.cutoff; envelope = true))
    end
    push!(cache.fits, CIMImageFit(key, Kind, precision(real(Result)), cim_logscale(integral.kernel),
        ComplexF64.(images), ComplexF64.(poles), ComplexF64(rotation),
        ComplexF64(Kind===:radial ? integral.weight.q : 0), certificates))
    cache.statistics.fits[]+=1
    return
end

function cim_reuse_estimate(integral::SpectralIntegral{Kind}, controls, workspace,
        evaluations, ::Type{Result}) where {Kind, Result}
    cache=cim_workspace(workspace)
    cache===nothing && return nothing
    integral.pole===nothing || return nothing
    spectral_cacheable_tail(spectral_tail(integral)) || return nothing
    key=cim_identity(integral.kernel)
    key===nothing && return nothing
    for fit in Iterators.reverse(cache.fits)
        fit.identity==key && fit.kind===Kind && fit.precision==precision(real(Result)) || continue
        length(fit.images)<=controls.max_terms || continue
        Kind===:radial && fit.shift!=integral.weight.q && continue
        Kind===:radial || fit.rotation==cis(integral.angle) || continue
        factor=exp(cim_logscale(integral.kernel)-fit.logscale)
        isfinite(factor)&&factor>0 || continue
        value, rounding=cim_image_value(Val(Kind), integral.weight, fit.images, fit.exponents, fit.rotation)
        value*=factor
        rounded=Result(value)
        rounding=rounding*factor+abs(value-rounded)
        isfinite(value) || continue
        target=max(controls.atol, controls.rtol*abs(value))
        for certificate in fit.certificates
            cim_contains(certificate, integral) || continue
            error=factor*certificate.error+rounding
            if error<=target
                cache.statistics.hits[]+=1
                return SpectralEstimate(rounded, error, factor*certificate.tail,
                    evaluations[], certificate.cutoff)
            end
        end
        # A changed geometry is certified against the full weighted residual,
        # not a finite collection of representative geometry samples.
        cutoff=maximum(c->c.cutoff, fit.certificates)
        cache.statistics.certifications[]+=1
        envelope=cim_envelope_error(integral, fit.images, fit.exponents,
            fit.rotation, factor, cutoff, target, controls)
        if isfinite(envelope.error)&&envelope.error+rounding<=target
            length(fit.certificates)>=32 && popfirst!(fit.certificates)
            push!(fit.certificates, cim_certificate(integral, envelope.error/factor,
                envelope.tail/factor, envelope.cutoff; envelope = true))
            cache.statistics.hits[]+=1
            return SpectralEstimate(rounded, envelope.error+rounding, envelope.tail,
                evaluations[], envelope.cutoff)
        end
    end
    return nothing
end

# Integrating the absolute residual times an analytic geometry envelope gives
# a certificate for every dominated weight on this same admissible contour.
function cim_envelope_error(integral::SpectralIntegral{Kind}, images, poles,
        rotation, factor, limit, target, controls) where {Kind}
    w=integral.weight
    contour=spectral_rotation(integral)
    for b in poles
        rate=Kind===:radial ? real((w.height+b)*contour)-w.separation*abs(imag(contour)) :
             real((w.height+b/rotation)*contour)-(w.separation+get(w, :radius, 0))*abs(imag(contour))
        rate>0 || return (error=Inf, tail=0.0, cutoff=Inf)
    end
    limit,remainder=cim_residual_limit(integral,images,poles,rotation,factor,limit,target,controls)
    physical=spectral_breakpoints(integral;limit)
    points=physical ./ (integral.scale .+ physical)
    isfinite(limit) || push!(points,1.0)
    f=t->begin
        den=inv(1-t)
        λ=contour*integral.scale*t*den
        if Kind===:radial
            u=sqrt(λ^2+w.q^2)
            coordinate=λ^2/(u+w.q)
            envelope=exp(-w.height*real(u)+w.separation*abs(imag(λ)))/abs(u)
            iszero(envelope) && return 0.0
            residual=radial_kernel_value(integral.kernel, u, λ)-factor*sum(
                i->images[i]*exp(-poles[i]*coordinate), eachindex(poles))
        else
            envelope=exp(-w.height*real(λ)+(w.separation+get(w, :radius, 0))*abs(imag(λ)))
            iszero(envelope) && return 0.0
            residual=integral.kernel(λ)-factor*sum(
                i->images[i]*exp(-poles[i]*λ/rotation), eachindex(poles))
        end
        return abs(residual)*envelope*integral.scale*den^2
    end
    value, error=try
        quadgk(f, points; rtol = 1e-3, atol = target/32, maxevals = controls.maxevals)
    catch error
        # An expansion can overflow outside its original geometry. Reject that
        # proposed reuse and let the normal construction resolve the new case.
        error isa DomainError || rethrow()
        return (error=Inf, tail=0.0, cutoff=Inf)
    end
    return (error=value+error+remainder, tail=remainder, cutoff=Float64(limit))
end

"""
$(TYPEDSIGNATURES)

Bound the absolute tail of the complete image expansion after `limit` \\[1/m\\].
Cosine and Bessel weights use the exponential envelope on the physical contour.
Radial images also bound the outgoing root and its reciprocal; an inadmissible
or insufficiently decaying image returns an infinite bound.
"""
function cim_image_tail(integral::SpectralIntegral{Kind}, images, poles, rotation, factor, limit) where {Kind}
    w=integral.weight
    contour=spectral_rotation(integral)
    result=0.0
    beta=0.0
    lower=1.0
    if Kind===:radial
        eta=abs(w.q^2)/limit^2
        eta<=1/4 || return Inf
        beta=abs(w.q^2)/(1+sqrt(1-eta))
        delta=beta/limit^2
        delta<real(contour)/2 || return Inf
        lower=1-delta
    end
    for (a,b) in zip(images,poles)
        iszero(a) && continue
        if Kind===:radial
            rate=real((w.height+b)*contour)-w.separation*abs(imag(contour))
            rate>0 || return Inf
            logvalue=log(abs(factor*a))+real(b*w.q)+abs(w.height+b)*beta/limit-
                     rate*limit-log(lower*rate*limit)
        else
            rate=real((w.height+b/rotation)*contour)-
                 (w.separation+get(w,:radius,0))*abs(imag(contour))
            rate>0 || return Inf
            logvalue=log(abs(factor*a))-rate*limit-log(rate)
        end
        result+=max(nextfloat(0.0),exp(logvalue))
    end
    return result*(1+32eps(Float64)*length(images))
end

function cim_residual_limit(integral,images,poles,rotation,factor,limit,target,controls)
    spectral_bounded(integral)&&isfinite(limit) || return (Inf,0.0)
    for _ in 1:controls.max_tail_refinements
        remainder=spectral_tail(integral)(limit)+cim_image_tail(integral,images,poles,rotation,factor,limit)
        remainder<=target/32 && return (limit,remainder)
        limit*=2
    end
    # A poor image continuation is rejected by the existing full-range
    # residual verifier; no truncated or quadrature value is substituted.
    return (Inf,0.0)
end
