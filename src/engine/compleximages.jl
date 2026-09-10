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
    angle_limit::typeof(Ref(0.0))
    "Construction, reuse and certification counters."
    statistics::NamedTuple{(:fits, :pencils, :hits, :certifications), NTuple{4, typeof(Ref(0))}}
end
CIMWorkspace() = CIMWorkspace(CIMImageFit[], ComplexF64[], ComplexF64[], Ref(pi/2),
    (fits = Ref(0), pencils = Ref(0), hits = Ref(0), certifications = Ref(0)))
cim_workspace(workspace) = workspace===nothing ? nothing : get(workspace, :cim, nothing)

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
    integral.features===nothing || integral.features.tail===nothing || return
    length(cache.fits)>=64 && popfirst!(cache.fits)
    certificate=cim_certificate(integral, error, tail, cutoff)
    certificates=[certificate]
    target=max(controls.atol, controls.rtol*abs(first(cim_image_value(
        Val(Kind), integral.weight, images, poles, rotation))))
    cache.statistics.certifications[]+=1
    envelope_error=cim_envelope_error(integral, images, poles, rotation, 1.0, cutoff, target, controls)
    if isfinite(envelope_error)
        push!(certificates, cim_certificate(integral, envelope_error, 0.0, cutoff; envelope = true))
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
    integral.features===nothing || integral.features.tail===nothing || return nothing
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
        error=cim_envelope_error(integral, fit.images, fit.exponents,
            fit.rotation, factor, cutoff, target, controls)
        if isfinite(error)&&error+rounding<=target
            length(fit.certificates)>=32 && popfirst!(fit.certificates)
            push!(fit.certificates, cim_certificate(integral, error/factor, 0.0,
                cutoff; envelope = true))
            cache.statistics.hits[]+=1
            return SpectralEstimate(rounded, error+rounding, 0.0, evaluations[], cutoff)
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
        rate>0 || return Inf
    end
    physical=sort!(unique!(vcat(spectral_breakpoints(integral), limit)))
    points=sort!(unique!(vcat(physical ./ (integral.scale .+ physical), 1.0)))
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
        return Inf
    end
    return value+error
end
