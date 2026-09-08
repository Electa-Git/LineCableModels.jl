"""
$(TYPEDEF)

Represent a real-axis spectral integral with a separately declared analytic
weight. `:cosine` integrates `kernel(λ) exp(-hλ) cos(yλ)`; `:radial` integrates
`kernel(u) exp(-hu) cos(yλ)/u`, where `u=sqrt(λ²+q²)` on the decaying branch.
The spectral coordinate and `scale` have units [1/m]; `h`, `y` have units [m].
CIM fits along the spectral coordinate at fixed physical frequency.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{Kind, K, W, T, P}
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
end

function SpectralIntegral(::Val{Kind}, kernel::K, weight::NamedTuple, scale::Real;
        angle = zero(scale), pole = nothing) where {Kind, K}
    isfinite(scale) && scale > 0 ||
        throw(DomainError(scale, "spectral scale must be positive and finite"))
    all(x -> isfinite(x) && x >= 0, (weight.height, weight.separation)) ||
        throw(ArgumentError("spectral weight lengths must be finite and nonnegative"))
    0 <= angle < pi/2 ||
        throw(ArgumentError("spectral fitting rotation must lie in [0,π/2)"))
    Kind === :radial && !iszero(angle) &&
        throw(ArgumentError("radial image fitting uses the declared u coordinate"))
    Kind === :radial && pole !== nothing &&
        throw(ArgumentError("simple-pole extraction requires a cosine weight"))
    pole === nothing ||
        (isfinite(pole.residue) && isfinite(pole.location) && !iszero(pole.location)) ||
        throw(ArgumentError("the extracted pole must have finite residue and nonzero finite location"))
    return SpectralIntegral{
        Kind, typeof(kernel), typeof(weight), typeof(scale), typeof(pole)}(
        kernel, weight, scale, angle, pole)
end

function (integral::SpectralIntegral{:cosine})(λ)
    w = integral.weight
    kernel = integral.kernel(λ)
    integral.pole === nothing ||
        (kernel += integral.pole.residue / (λ + integral.pole.location))
    return kernel * exp(-w.height * λ) * cos(w.separation * λ)
end

function (integral::SpectralIntegral{:radial})(λ)
    w = integral.weight
    u = sqrt(λ^2 + w.q^2)
    return integral.kernel(u) * exp(-w.height * u) / u * cos(w.separation * λ)
end

"""
$(TYPEDSIGNATURES)

Evaluate one declared spectral integral with the selected numerical algorithm.
Controls come from `computation_options`; mutable samples and image coefficients
belong to the supplied computation workspace. Nonconvergence raises an error.
The reported algorithm always produces the returned value.
"""
function integrate(::Val{:quad}, integral::SpectralIntegral, controls, workspace)
    R = typeof(integral.scale)
    magnitude(z) = abs(complex(nominal(real(z)), nominal(imag(z))))
    scale = integral.scale
    segments = workspace === nothing ? nothing : workspace.segments
    segments === nothing || empty!(segments)
    # Scaling resolves kernels whose transition is far from λ=1 [1/m].
    value,
    error = quadgk(x -> integral(x * scale) * scale, zero(R),
        (integral.pole === nothing ? one(R) : R(abs(integral.pole.location)) / scale),
        R(Inf);
        rtol = controls.rtol, atol = controls.atol, maxevals = controls.maxevals,
        segbuf = segments, norm = magnitude)
    isfinite(value) && isfinite(error) &&
    error <= max(controls.atol, controls.rtol * magnitude(value)) ||
        throw(ErrorException("spectral :quad did not converge (estimated error=$error)"))
    return value
end

function integrate(::Val{:trapz}, integral::SpectralIntegral, controls, workspace)
    R = typeof(integral.scale)
    magnitude(z) = abs(complex(nominal(real(z)), nominal(imag(z))))
    scale = integral.pole === nothing ? integral.scale :
            min(integral.scale, R(abs(integral.pole.location)))
    # Logarithmic mapping resolves widely separated material and geometric scales.
    rotation = iszero(integral.angle) ? one(R) : cis(integral.angle)
    transformed(t) = integral(rotation * scale * expm1(t)) * rotation * scale * exp(t)
    previous_tail = nothing
    for tail in 1:controls.max_tail_refinements
        upper = R(tail) * log(R(4))
        n = controls.samples
        previous_grid = nothing
        converged = false
        value = transformed(zero(R)) * zero(R)
        for refinement in 0:controls.max_refinements
            step = one(R) / n
            value = zero(value)
            for index in 1:(n - 1)
                v = step * index
                u = upper * sinpi(v / 2)^2
                jacobian = upper * (R(π) / 2) * sinpi(v)
                value += transformed(u) * jacobian
            end
            value *= step
            isfinite(value) ||
                throw(ErrorException("nonfinite transformed trapezoidal integral"))
            if previous_grid !== nothing &&
               magnitude(value - previous_grid) <=
               max(controls.atol, controls.rtol * magnitude(value)) / 4
                converged = true
                break
            end
            previous_grid = value
            n *= 2
        end
        converged ||
            throw(ErrorException("spectral :trapz grid refinement did not converge"))
        if previous_tail !== nothing &&
           magnitude(value - previous_tail) <=
           max(controls.atol, controls.rtol * magnitude(value)) / 2
            return value
        end
        previous_tail = value
    end
    throw(ErrorException("spectral :trapz tail refinement did not converge"))
end

"""
$(TYPEDSIGNATURES)

Approximate a spectral kernel by complex exponentials using matrix pencils on
nested spectral intervals. Fit amplitudes against independent logarithmic and
uniform samples, then integrate the images analytically. For the cosine weight
an image contributes `a*(h+b)/((h+b)^2+y^2)`; the radial Sommerfeld weight gives
`a*K₀(q*sqrt((h+b)^2+y^2))`. These are different analytic identities.

Fits are local to this integral and frequency. Matrix-pencil factorization
currently admits Float32/Float64 kernels; uncertainty and higher precision
inputs fail explicitly. A held-out residual and an independent quadrature
comparison must both pass. Quadrature is a validation calculation and never
supplies a value returned as `:cim`.

The matrix-pencil/GPOF construction follows the spectral exponential approach
used in discrete complex images; see Rallis, doctoral thesis (National Archive record 10442/34633), and
Y. P. Chen, W. C. Chew and L. Jiang, IEEE AWPL 10 (2011), 419–422,
DOI 10.1109/LAWP.2011.2152358.
"""
function integrate(::Val{:cim}, integral::SpectralIntegral{Kind}, controls, workspace) where {Kind}
    Kind in (:cosine, :radial) ||
        throw(ArgumentError("CIM has no analytic identity for weight :$Kind"))
    if workspace !== nothing
        empty!(workspace.images)
        empty!(workspace.exponents)
    end
    first_value = integral.kernel(Kind === :radial ? integral.weight.q :
                                  zero(integral.scale))
    T = typeof(float(real(first_value)))
    T in (Float32, Float64) || throw(ArgumentError(
        "CIM matrix pencils require Float32 or Float64 physical inputs; got $T"))
    fit_tolerance = max(controls.rtol, 100eps(T))
    R = typeof(integral.scale)
    w = integral.weight
    magnitude(z) = abs(z)
    shift = Kind === :radial ? w.q : zero(first_value)
    rotation = Kind === :cosine ? cis(integral.angle) : one(shift)
    kernel(x) = integral.kernel(rotation*x + shift)
    amplitude = max(
        maximum(x -> abs(kernel(x)),
            (zero(R), integral.scale / 16, integral.scale, 16integral.scale)),
        floatmin(T))
    limit = integral.scale
    # Establish a finite fitting interval where the weighted kernel has decayed.
    decayed = false
    for step in 1:controls.max_tail_refinements
        limit *= 4
        if exp(log(abs(kernel(limit))) - real(w.height * (rotation*limit + shift)) +
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
    reference = integrate(Val(:quad), integral,
        (rtol = controls.rtol / 10, atol = controls.atol / 10,
            maxevals = controls.maxevals), nothing)
    target = max(controls.atol, controls.rtol * abs(reference))
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
    previous = nothing
    diagnostic = nothing
    for refinement in 0:controls.max_refinements
        samples = controls.samples * 2^refinement
        intervals = max(1, ceil(Int, log(limit / integral.scale) / log(16)))
        widths = exp.(range(log(integral.scale), log(limit); length = intervals + 1))
        poles = Complex{T}[]
        sampled_nonzero = false
        terms = max(2, controls.max_terms ÷ length(widths))
        for width in widths
            origin = width > 2integral.scale ? width / 16 : zero(width)
            step = (width-origin) / (samples - 1)
            data = Complex{T}[kernel(origin + step * i) for i in 0:(samples - 1)]
            sampled_nonzero |= any(!iszero, data)
            rows = samples ÷ 2
            columns = samples - rows
            H0 = [data[i + j - 1] for i in 1:rows, j in 1:columns]
            H1 = [data[i + j] for i in 1:rows, j in 1:columns]
            pencil = svd(H0)
            rank = min(terms, count(
                >(max(eps(T) * samples, fit_tolerance * 1e-3) *
                  first(pencil.S)), pencil.S))
            rank == 0 && continue
            U = @view pencil.U[:, 1:rank]
            V = @view pencil.V[:, 1:rank]
            reduced = U' * H1 * V * Diagonal(inv.(pencil.S[1:rank]))
            for z in eigvals(reduced)
                iszero(z) && continue
                b = -log(z) / step
                isfinite(b) && real(b) > zero(T) && length(poles) < controls.max_terms &&
                    push!(poles, b)
            end
        end
        if isempty(poles)
            !sampled_nonzero && abs(pole_value-reference) <= target && return pole_value
            throw(ErrorException("CIM matrix pencils produced no decaying images"))
        end
        # Logarithmic samples resolve material and geometric scales separately.
        nodes = sort!(unique!(vcat(zero(R),
            exp.(range(
                log(integral.scale / samples), log(limit); length = 4samples)),
            collect(range(zero(R), limit; length = samples)))))
        values = Complex{T}[kernel(x) for x in nodes]
        basis = [exp(-b * x) for x in nodes, b in poles]
        weights = inv.(sqrt.(max.(abs.(values), amplitude * fit_tolerance * 1e-4)))
        decomposition = svd(basis .* weights)
        rank = count(>(eps(T) * max(size(basis)...) * first(decomposition.S)), decomposition.S)
        amplitudes = decomposition.V[:, 1:rank] *
                     ((decomposition.U[:, 1:rank]' * (values .* weights)) ./
                      decomposition.S[1:rank])
        held = exp.(range(log(integral.scale / (samples + 1)), log(limit);
            length = 3samples + 1))
        residual = maximum(
            x -> abs(sum(amplitudes .* exp.(-poles .* x)) - kernel(x)) *
                 exp(-real(w.height * (rotation*x + shift)) +
                     abs(imag(w.separation*rotation))*x),
            held)
        images = amplitudes .* exp.(poles .* shift)
        value = if Kind === :cosine
            sum(rotation .* images .* ((w.height*rotation .+ poles) ./
                 ((w.height*rotation .+ poles) .^ 2 .+ (w.separation*rotation)^2)))
        else
            sum(images .* map(
                b -> special_besselk(0,
                    w.q * sqrt((w.height + b)^2 + w.separation^2)), poles))
        end
        value += pole_value
        if isfinite(value) && residual <= sqrt(fit_tolerance) * amplitude &&
           abs(value - reference) <= target &&
           (previous === nothing || abs(value - previous) <= 2target)
            if workspace !== nothing
                resize!(workspace.images, length(images));
                copyto!(workspace.images, images)
                resize!(workspace.exponents, length(poles));
                copyto!(workspace.exponents, poles)
            end
            return value
        end
        diagnostic = (residual = residual/amplitude, integral_error = abs(value-reference),
            target = target, images = length(poles))
        previous = value
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
        (rtol = 1e-6, atol = 0.0, samples = 128,
            max_refinements = 12, max_tail_refinements = 32)
    else
        (rtol = 1e-6, atol = 0.0, samples = 128, max_terms = 192, max_refinements = 3,
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
        value isa Integer && !(value isa Bool) && value > 0 ||
            throw(ArgumentError("$key must be a positive integer"))
    end
    haskey(controls, :samples) && controls.samples < 16 &&
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
