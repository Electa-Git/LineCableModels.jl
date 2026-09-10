"""
$(TYPEDEF)

Represent a real-axis spectral integral with a separately declared analytic
weight. `:cosine` integrates `kernel(λ) exp(-hλ) cos(yλ)`; `:radial` integrates
`kernel(u) exp(-hu) cos(yλ)/u`, where `u=sqrt(λ²+q²)` on the decaying branch.
The spectral coordinate and `scale` have units [1/m]; `h`, `y` have units [m].
CIM fits along the spectral coordinate at fixed physical frequency.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{Kind, K, W, T, P, F}
    "Residual spectral kernel, excluding the declared analytic weight."
    kernel::K
    "Analytic weight parameters: height, separation, and radial q where applicable."
    weight::W
    "Positive real spectral scale [1/m]; the path is λ ∈ [0,∞)."
    scale::T
    "Admissible upper-half-plane rotation for cosine-kernel CIM fitting [rad]."
    angle::T
    "Optional exact simple pole, residue/(λ+location), separated from the kernel."
    pole::P
    "Optional typed spectral feature and tail metadata."
    features::F
end

# Kernel owners may supply a bound for the complete analytic weight. Arbitrary
# callbacks keep the independently verified, infinite-range numerical path.
spectral_kernel_tail(kernel, weight, angle) = nothing

function SpectralIntegral(::Val{Kind}, kernel::K, weight::NamedTuple, scale::Real;
        angle = zero(scale), pole = nothing, features = nothing) where {Kind, K}
    isfinite(scale) && scale > 0 ||
        throw(DomainError(scale, "spectral scale must be positive and finite"))
    all(x -> isfinite(x) && x >= 0, (weight.height, weight.separation)) ||
        throw(ArgumentError("spectral weight lengths must be finite and nonnegative"))
    if Kind===:besselcosine
        haskey(weight, :radius)&&isfinite(weight.radius)&&weight.radius>=0 ||
            throw(ArgumentError("Bessel weight radius must be finite and nonnegative"))
    end
    0 <= angle < pi/2 ||
        throw(ArgumentError("spectral fitting rotation must lie in [0,π/2)"))
    Kind === :radial && pole !== nothing &&
        throw(ArgumentError("simple-pole extraction requires a cosine weight"))
    pole === nothing ||
        (isfinite(pole.residue) && isfinite(pole.location) && !iszero(pole.location)) ||
        throw(ArgumentError("the extracted pole must have finite residue and nonzero finite location"))
    if features !== nothing && features.tail === nothing && pole === nothing
        tail = spectral_kernel_tail(kernel, weight, angle)
        features = SpectralFeatures{eltype(features.points), typeof(tail)}(features.points, tail)
    end
    return SpectralIntegral{
        Kind, typeof(kernel), typeof(weight), typeof(scale), typeof(pole), typeof(features)}(
        kernel, weight, scale, typeof(scale)(angle), pole, features)
end

@inline spectral_magnitude(z) = abs(complex(nominal(real(z)), nominal(imag(z))))

@inline function spectral_cosine(height, separation, λ)
    # Avoid 0*Inf on an admissible complex contour.
    return (exp(-(height+im*separation)*λ)+exp(-(height-im*separation)*λ))/2
end

function (integral::SpectralIntegral{:cosine})(λ)
    w = integral.weight
    if integral.kernel isa Union{EarthSpectrum, EarthPathVoltageSpectrum} &&
       earth_combined_weight(integral.kernel)
        return (earth_weighted_spectrum(integral.kernel, λ, -(w.height+im*w.separation)*λ) +
                earth_weighted_spectrum(integral.kernel, λ, -(w.height-im*w.separation)*λ))/2
    end
    kernel = integral.kernel(λ)
    integral.pole === nothing ||
        (kernel += integral.pole.residue / (λ + integral.pole.location))
    return kernel * spectral_cosine(w.height, w.separation, λ)
end

function (integral::SpectralIntegral{:besselcosine})(λ)
    w=integral.weight
    growth=abs(imag(w.radius*λ))
    weight=(exp(-(w.height+im*w.separation)*λ+growth) +
            exp(-(w.height-im*w.separation)*λ+growth))/2
    weighted=if integral.kernel isa EarthSpectrum&&earth_combined_weight(integral.kernel)
        (earth_weighted_spectrum(integral.kernel, λ, -(w.height+im*w.separation)*λ+growth) +
         earth_weighted_spectrum(integral.kernel, λ, -(w.height-im*w.separation)*λ+growth))/2
    else
        integral.kernel(λ)*weight
    end
    isfinite(weighted) || throw(DomainError((λ, weighted, w),
        "nonfinite weighted Bessel spectral kernel"))
    iszero(spectral_magnitude(weighted)) && return weighted
    return weighted*special_besseljx(0, w.radius*λ)
end

function (integral::SpectralIntegral{:radial})(λ)
    w = integral.weight
    u = sqrt(λ^2 + w.q^2)
    return radial_kernel_value(integral.kernel, u, λ)/u *
           (exp(-w.height*u+im*w.separation*λ)+exp(-w.height*u-im*w.separation*λ))/2
end

# Kernels may use the original physical coordinate to avoid subtracting two
# nearly equal squared roots when evaluating their radial representation.
@inline radial_kernel_value(kernel, u, λ) = kernel(u)

struct SpectralKernelCounter{K, R}
    kernel::K
    evaluations::R
end
@inline function (counter::SpectralKernelCounter)(x)
    counter.evaluations[]+=1
    return counter.kernel(x)
end
@inline function radial_kernel_value(counter::SpectralKernelCounter, u, λ)
    counter.evaluations[]+=1
    return radial_kernel_value(counter.kernel, u, λ)
end

"""
$(TYPEDSIGNATURES)

Evaluate one declared spectral integral with the selected numerical algorithm.
Controls come from `computation_options`; mutable samples and image coefficients
belong to the supplied computation workspace. Nonconvergence raises an error.
The reported algorithm always produces the returned value.
"""
function integrate(::Val{:quad}, integral::SpectralIntegral, controls, workspace)
    return spectral_estimate(Val(:quad), integral, controls, workspace).value
end

function integrate(::Val{:trapz}, integral::SpectralIntegral, controls, workspace)
    return spectral_estimate(Val(:trapz), integral, controls, workspace).value
end

"""
$(TYPEDSIGNATURES)

Approximate a spectral kernel by complex exponentials using matrix pencils on
adaptive spectral regions. Fit amplitudes against independent logarithmic and
uniform samples, then integrate the images analytically. For the cosine weight
an image contributes `a*(h+b)/((h+b)^2+y^2)`; the radial Sommerfeld weight gives
`a*K₀(q*sqrt((h+b)^2+y^2))`. These are different analytic identities.

Built-in earth kernels reuse workspace-owned fits when their complete material
identity matches and a geometry certificate meets the requested error budget.
New fits require independent quadrature and full weighted-residual verification;
new geometry envelopes require full residual certification. Cache hits evaluate
only the images. Arbitrary callables without a declared identity are not cached.
Matrix-pencil factorization admits Float32/Float64 kernels; uncertainty and
higher precision inputs fail explicitly. Quadrature never supplies a value
returned as `:cim`.

The matrix-pencil/GPOF construction follows the spectral exponential approach
used in discrete complex images; see Rallis, doctoral thesis (National Archive record 10442/34633), and
Y. P. Chen, W. C. Chew and L. Jiang, IEEE AWPL 10 (2011), 419–422,
DOI 10.1109/LAWP.2011.2152358.
"""
function integrate(::Val{:cim}, integral::SpectralIntegral, controls, workspace)
    estimate=spectral_estimate(Val(:cim), integral, controls, workspace)
    target=max(controls.atol, controls.rtol*spectral_magnitude(estimate.value))
    estimate.error<=target || throw(ErrorException(
        "spectral :cim did not converge (estimated error=$(estimate.error), target=$target); no quadrature fallback was used"))
    return estimate.value
end

function spectral_estimate(::Val{:cim}, integral::SpectralIntegral{Kind}, controls, workspace) where {Kind}
    original=integral.kernel
    prototype=original(Kind===:radial ? integral.weight.q : zero(integral.scale))
    physical_type=typeof(float(real(prototype)))
    physical_type in (Float32, Float64) || throw(ArgumentError(
        "CIM matrix pencils require Float32 or Float64 physical inputs; got $physical_type"))
    weight_types=map(value->typeof(float(real(value))), values(integral.weight))
    all(type->type in (Float32, Float64), weight_types) || throw(ArgumentError(
        "CIM matrix pencils require Float32 or Float64 physical weights; got $weight_types"))
    evaluations=Ref(1)
    counted=SpectralKernelCounter(original, evaluations)
    # Fit and verify in Float64 also for Float32 physical data. Float32
    # sampling noise otherwise overwhelms small singular values and makes the
    # absolute residual verifier refine roundoff. Account for output rounding.
    sampled=SpectralIntegral(Val(Kind), counted, integral.weight, Float64(integral.scale);
        angle = integral.angle, pole = integral.pole, features = integral.features)
    reused=cim_reuse_estimate(sampled, controls, workspace, evaluations, Complex{physical_type})
    # Keep the physical scalar contract explicit across the cache/construction
    # join; all fitting and numerical-error arithmetic uses Float64.
    reused===nothing || return reused::SpectralEstimate{Complex{physical_type}, Float64}
    return cim_estimate(sampled, controls, workspace, evaluations,
        Complex{physical_type})::SpectralEstimate{Complex{physical_type}, Float64}
end

function cim_estimate(integral::SpectralIntegral{Kind}, controls, workspace, evaluations,
        ::Type{Result}) where {Kind, Result}
    Kind in (:cosine, :besselcosine, :radial) ||
        throw(ArgumentError("CIM has no analytic identity for weight :$Kind"))
    if workspace !== nothing
        empty!(workspace.images)
        empty!(workspace.exponents)
    end
    first_value = zero(Result)
    T = Float64
    fit_tolerance = max(controls.rtol, 100eps(T))
    R = typeof(integral.scale)
    w = integral.weight
    magnitude(z) = abs(z)
    shift = Kind === :radial ? w.q : zero(first_value)
    rotation = Kind === :radial ? one(shift) : cis(integral.angle)
    kernel(x) = integral.kernel(rotation*x + shift)
    budget=SpectralSampleBudget(controls.samples)
    construct(x)=spectral_sample!(budget,kernel,x)
    verification=spectral_reference(integral,controls)
    reference=verification.value
    target=max(controls.atol,controls.rtol*abs(reference))
    bounded=spectral_bounded(integral)&&isfinite(verification.cutoff)
    limit=bounded ? verification.cutoff : integral.scale
    amplitude_nodes=spectral_breakpoints(integral;limit=bounded ? limit : Inf)
    if Kind===:radial
        amplitude_nodes=R[abs(x^2/(sqrt(x^2+w.q^2)+w.q)) for x in amplitude_nodes]
    end
    append!(amplitude_nodes, (zero(R), integral.scale/16, integral.scale))
    bounded && filter!(x->x<=limit, amplitude_nodes)
    bounded || push!(amplitude_nodes,16integral.scale)
    amplitude = max(
        maximum(x -> abs(construct(x)), amplitude_nodes),
        floatmin(T))
    # Establish a finite fitting interval where the weighted kernel has decayed.
    decayed = bounded
    for step in 1:controls.max_tail_refinements
        decayed && break
        limit *= 4
        if exp(log(abs(construct(limit))) - real(w.height * (rotation*limit + shift)) +
               abs(imag(w.separation * rotation)) * limit) <=
           fit_tolerance * amplitude / 100
            decayed = true
            break
        end
    end
    decayed ||
        throw(ErrorException("CIM spectral fitting interval did not reach a decaying tail"))
    # Independent integral validation also detects errors between sampled points
    # and use of a transformed spectral fit outside its verified domain.
    tail_error=bounded ? verification.tail : R(Inf)
    for tail in 1:controls.max_tail_refinements
        bounded && break
        tail_error=cim_true_tail(integral, limit,
            Kind===:radial ? spectral_rotation(integral) : rotation, target, controls)
        tail_error<=target/8 && break
        tail==controls.max_tail_refinements &&
            throw(ErrorException("CIM could not control the integrated spectral tail"))
        limit*=2
    end
    fit_tolerance = min(fit_tolerance,
        max(eps(T), target / max(amplitude*integral.scale, floatmin(T))))
    pole_value = zero(reference)
    if integral.pole !== nothing
        p = integral.pole
        for length in (complex(w.height, w.separation), complex(w.height, -w.separation))
            z = length * p.location
            term = SpecialFunctions.expintx(z)
            # Continue E₁ along the Laplace parameter, rather than resetting its
            # logarithm to the principal branch when the product crosses the cut.
            phase = angle(length) + angle(p.location)
            phase > π && (term -= 2π * im * exp(z))
            phase < -π && (term += 2π * im * exp(z))
            pole_value += p.residue * term / 2
        end
    end
    best_value=zero(reference)
    best_error=R(Inf)
    best_cutoff=R(Inf)
    best_tail=zero(R)
    diagnostic = nothing
    pencils=CIMPencil[]
    pencil_plan=(false,0)
    constant_pole=false
    sampled_nonzero=false
    fit_data=nothing
    for refinement in 0:(controls.max_refinements+1)
        # Start with a compact collection of regions. Expand coverage before
        # increasing the uniform pencil resolution on difficult kernels.
        compact=refinement<=1
        samples = controls.samples===nothing ? 128*2^max(0,refinement-2)+1 :
                  min(129,max(16,controls.samples÷32))
        widths = cim_windows(integral, limit)
        if compact && length(widths)>8
            widths=widths[unique(round.(Int, range(1, length(widths); length = 8)))]
        end
        poles = Complex{T}[]
        term_limit=refinement==0 ? min(96, controls.max_terms) : controls.max_terms
        terms = max(4, term_limit ÷ length(widths))
        origins = vcat(zero(R), widths[1:(end - 1)])
        windows=compact ? collect(zip(origins, widths)) :
                vcat(collect(zip(origins, widths)), [(zero(R), width) for width in widths[2:end]])
        cache=cim_workspace(workspace)
        if pencil_plan!=(compact,samples)
            empty!(pencils)
            fit_data=nothing
            constant_pole=false
            sampled_nonzero=false
            for (origin,width) in windows
                step=(width-origin)/(samples-1)
                data=cache===nothing ? Vector{Complex{T}}(undef,samples) : resize!(cache.samples,samples)
                for i in 1:samples
                    data[i]=construct(origin+step*(i-1))
                end
                sampled_nonzero |= any(!iszero,data)
                variation=maximum(z->abs(z-first(data)),data)
                if variation<=max(8eps(T),fit_tolerance/100)*amplitude
                    constant_pole=true
                else
                    push!(pencils,cim_pencil(data,step,width,fit_tolerance,
                        min(controls.max_terms,max(32,terms)),cache))
                end
            end
            pencil_plan=(compact,samples)
        end
        constant_pole && push!(poles,zero(Complex{T}))
        cim_poles!(poles,pencils,terms,integral,rotation)
        if !compact
            # Supplement an unresolved pencil basis with real decays covering
            # every declared scale. These independent candidates can resolve a
            # smooth residual when near-duplicate complex pencil poles make
            # its representation cancellation-sensitive.
            minimum_rate=inv(limit)
            maximum_rate=8/first(widths)
            rate_count=clamp(ceil(Int, log(maximum_rate/minimum_rate)/log(sqrt(2))), 16, 256)
            append!(poles, Complex{T}.(exp.(range(log(minimum_rate), log(maximum_rate); length = rate_count))))
        end
        if isempty(poles)
            if !sampled_nonzero && abs(pole_value-reference) <= target
                rounded=Result(pole_value)
                return SpectralEstimate(rounded,
                    R(abs(rounded-pole_value)+abs(pole_value-reference)+verification.error),
                    R(tail_error), evaluations[], R(limit), budget.used[])
            end
            throw(ErrorException("CIM matrix pencils produced no decaying images"))
        end
        # Logarithmic samples resolve material and geometric scales separately.
        if fit_data===nothing
            nodes = sort!(unique!(vcat(zero(R), spectral_breakpoints(integral;limit),
            exp.(range(
                log(first(widths) / samples), log(limit); length = 4samples)),
            [origin+(width-origin)*i/16 for (origin, width) in zip(origins, widths)
             for i in 0:16])))
            fit_data=cim_fit_samples(integral,nodes,rotation,shift,budget)
        end
        coordinates,values,weights=fit_data
        basis = [exp(-b * coordinates[i]+log(weights[i]))
                 for i in eachindex(coordinates), b in poles]
        column_scales=[maximum(abs, column) for column in eachcol(basis)]
        finite_columns=findall(x->isfinite(x)&&x>0, column_scales)
        isempty(finite_columns)&&throw(ErrorException(
            "CIM candidates overflow on the physical fitting contour"))
        if length(finite_columns)!=length(poles)
            # A pole can decay at infinity yet have unrepresentable growth on
            # the curved radial contour. Such candidates cannot enter LAPACK;
            # the complete residual still verifies the remaining image basis.
            poles=poles[finite_columns]
            basis=basis[:, finite_columns]
            column_scales=column_scales[finite_columns]
        end
        if maximum(column_scales)>1e30
            basis ./= transpose(column_scales)
        else
            fill!(column_scales, one(T))
        end
        if length(poles)>term_limit
            # A first-come cap discarded the large-scale windows when many
            # small material features were present. Select a diverse basis
            # from every window using the physical weighted fitting metric.
            selected=qr(basis, ColumnNorm()).p[1:term_limit]
            poles=poles[selected]
            basis=basis[:, selected]
            column_scales=column_scales[selected]
        end
        # Penalize each image's analytic continuation beyond the fitting
        # interval. A nearly growing image must not hide between the final
        # training point and its very long weighted tail.
        penalty=cim_tail_penalty(integral, poles, rotation, limit) ./ column_scales
        decomposition = svd(vcat(basis, Diagonal(penalty)))
        rank = count(
            >(eps(T) * max(sqrt(max(size(basis)...)),
                  max(size(basis)...)/4^refinement) * first(decomposition.S)),
            decomposition.S)
        amplitudes = decomposition.V[:, 1:rank] *
                     ((decomposition.U[:, 1:rank]' *
                       vcat(values .* weights, zeros(Complex{T}, length(poles)))) ./
                      decomposition.S[1:rank])
        amplitudes ./= column_scales
        held = exp.(range(log(first(widths) / (samples + 1)), log(limit);
            length = 3samples + 1))
        residual = maximum(
            x -> abs(sum(i->amplitudes[i]*exp(-poles[i]*x), eachindex(poles)) - kernel(x)) *
                 exp(-real(w.height * (rotation*x + shift)) +
                     abs(imag(w.separation*rotation))*x),
            held)
        images = Kind===:radial ? amplitudes : amplitudes .* exp.(poles .* shift)
        value, image_rounding=cim_image_value(Val(Kind), w, images, poles, rotation)
        value += pole_value
        residual_estimate = isfinite(value) &&
                              (abs(value-reference)<=target ||
                               refinement==controls.max_refinements+1) ?
                              cim_residual_error(
            integral, amplitudes, poles, rotation, shift, limit, target, controls) :
            (error=R(Inf), tail=zero(R), cutoff=R(Inf))
        integrated_residual=residual_estimate.error
        # The complete residual certificate already bounds the accepted fit.
        # Its distance from an earlier rejected fit adds no accuracy evidence.
        if isfinite(value) && integrated_residual+verification.error+image_rounding <= target &&
           abs(value - reference) <= target
            if workspace !== nothing
                resize!(workspace.images, length(images))
                copyto!(workspace.images, images)
                resize!(workspace.exponents, length(poles))
                copyto!(workspace.exponents, poles)
            end
            rounded=Result(value)
            cim_store_fit!(integral, controls, workspace, images, poles, rotation,
                integrated_residual+verification.error, residual_estimate.tail,
                residual_estimate.cutoff, Result)
            return SpectralEstimate(
                rounded, R(abs(rounded-value)+integrated_residual+verification.error+image_rounding),
                R(residual_estimate.tail), evaluations[], residual_estimate.cutoff, budget.used[])
        end
        diagnostic = (residual = residual/amplitude, integrated_residual,
            integral_error = abs(value-reference),
            target = target, images = length(poles))
        if integrated_residual+verification.error+image_rounding<best_error
            best_value=value
            best_error=integrated_residual+verification.error+image_rounding
            best_cutoff=residual_estimate.cutoff
            best_tail=residual_estimate.tail
        end
        if refinement==controls.max_refinements+1 && isfinite(integrated_residual)
            # A system consumer may allocate a larger absolute budget to an
            # insignificant correction. The value-only wrapper still enforces
            # its requested integral tolerance.
            rounded=Result(best_value)
            return SpectralEstimate(rounded, R(abs(rounded-best_value)+best_error),
                R(best_tail), evaluations[], R(best_cutoff), budget.used[])
        end
    end
    throw(ErrorException("spectral :cim did not converge: $diagnostic; no quadrature fallback was used"))
end

"""
$(TYPEDSIGNATURES)

Normalize a spectral integration algorithm and its controls. The selected
method is `:quad`, `:trapz`, or `:cim`. Controls belong to that method; changing
methods does not inherit controls from the previous method. CIM samples the
spectral variable at fixed physical frequency.
"""
function computation_options(::Type{SpectralIntegral}, values::NamedTuple)
    isempty(setdiff(keys(values), (:method, :options))) ||
        throw(ArgumentError("integration accepts only method and options"))
    method = get(values, :method, :quad)
    method in (:quad, :trapz, :cim) ||
        throw(ArgumentError("integration method must be :quad, :trapz, or :cim"))
    supplied = get(values, :options, (;))
    supplied isa NamedTuple ||
        throw(ArgumentError("integration options must be a named tuple"))
    defaults = if method === :quad
        (rtol = 1e-8, atol = 0.0, maxevals = 10^7)
    elseif method === :trapz
        (rtol = 1e-6, atol = 0.0, samples = nothing,
            max_refinements = 12, max_tail_refinements = 32)
    else
        (rtol = 1e-6, atol = 0.0, samples = nothing, max_terms = 192, max_refinements = 3,
            max_tail_refinements = 20, maxevals = 10^7)
    end
    unknown_controls = setdiff(keys(supplied), keys(defaults))
    isempty(unknown_controls) ||
        throw(ArgumentError("unknown :$method integration controls: $(collect(unknown_controls))"))
    controls = merge(defaults, supplied)
    for key in (:rtol, :atol)
        value = getproperty(controls, key)
        value isa Real && isfinite(value) && value >= 0 ||
            throw(ArgumentError("$key must be finite and nonnegative"))
    end
    controls.rtol > 0 || controls.atol > 0 ||
        throw(ArgumentError("at least one integration tolerance must be positive"))
    for key in setdiff(keys(controls), (:rtol, :atol))
        value = getproperty(controls, key)
        key===:samples && value===nothing && continue
        value isa Integer && !(value isa Bool) && value > 0 ||
            throw(ArgumentError("$key must be a positive integer"))
    end
    haskey(controls, :samples) && controls.samples!==nothing && controls.samples < 16 &&
        throw(ArgumentError("integration samples must be at least 16"))
    return (method = Val(method), options = controls)
end

function computation_options(::FormulaMethod, ::Val{:integration}, defaults::NamedTuple,
        supplied::NamedTuple)
    isempty(setdiff(keys(supplied), (:method, :options))) ||
        throw(ArgumentError("integration accepts only method and options"))
    get(supplied, :options, (;)) isa NamedTuple ||
        throw(ArgumentError("integration options must be a NamedTuple"))
    method = get(supplied, :method, defaults.method)
    controls = method === defaults.method ?
               merge(defaults.options, get(supplied, :options, (;))) :
               get(supplied, :options, (;))
    return computation_options(SpectralIntegral, (; method, options = controls))
end
