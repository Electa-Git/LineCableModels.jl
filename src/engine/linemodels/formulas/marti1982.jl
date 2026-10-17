import ..ModalAnalysis: Tv, orient_modes!

"""
$(TYPEDSIGNATURES)

**Identification.** J. Marti constant real modal transformation (1982).

**Expression.** Modal characteristic impedance and propagation magnitudes have
Bode asymptotes with slopes in multiples of 20 dB per decade. Their corners
supply starting real poles and zeros. Levenberg–Marquardt least squares refines
positive pole magnitudes, zero magnitudes and gain in logarithmic coordinates:

```math
F(s)=k\\frac{\\prod_j(s+z_j)}{\\prod_j(s+p_j)},
\\qquad k>0,\\quad z_j>0,\\quad p_j>0.
```

The poles are ``-p_j`` and the zeros are ``-z_j`` \\[s⁻¹\\]. The product is
minimum phase by construction. Characteristic impedance has equal numbers of
poles and zeros. Lossy propagation is strictly proper and uses the full band
with a magnitude weight floor.

Factor products give real residues directly, without expanding polynomials.
Propagation retains a separate travel delay \\[s\\] consistent with the final
rational phase.

The options select the transformation frequency \\[Hz\\], Bode tolerance \\[dB\\]
and propagation magnitude weight floor \\[dimensionless\\].
"""
description(::Type{<:Formula{:marti}}; compact::Bool=false) =
    compact ? "Marti" : "Constant-transformation modal line model (1982)"

function formulation_options(::Expression{<:Formula{:marti}, typeof(fit)})
    return FormulationOptions((transformation=(frequency=nothing,),
        bode=(tolerance_db=0.5, threshold=1e-2)))
end

function formulation_options(expression::Expression{<:Formula{:marti}, typeof(fit)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :marti formulation options"))
    sections = map(keys(defaults)) do name
        explicit = get(supplied.data, name, (;))::NamedTuple
        section = getproperty(defaults, name)
        isempty(setdiff(keys(explicit), keys(section))) ||
            throw(ArgumentError("unknown :marti $name options"))
        return merge(section, explicit)
    end
    return FormulationOptions(NamedTuple{keys(defaults)}(sections))
end

function validate(options::FormulationOptions, ::Formula{:marti})
    tolerance = options.data.bode.tolerance_db::Real
    threshold = options.data.bode.threshold::Real
    frequency = options.data.transformation.frequency
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("Bode tolerance must be finite and positive"))
    0 < threshold < 1 ||
        throw(ArgumentError("Bode threshold must lie strictly between zero and one"))
    frequency === nothing || (isfinite(frequency::Real) && frequency > 0) ||
        throw(ArgumentError("transformation frequency must be finite and positive"))
    return options
end

function initialize_buffers(::Formula{:marti}, ::Type{T}, input, plan,
        buffers) where {T}
    R = eltype(input.frequencies)
    n = length(input.coordinates)
    count = length(input.frequencies)
    section = NamedTuple{(:constant, :poles, :residues, :delay),
        Tuple{R, Vector{R}, Vector{R}, R}}
    return merge(buffers, (omega=2R(π) .* input.frequencies,
        impedance=Matrix{Complex{R}}(undef, n, count),
        propagation=Matrix{Complex{R}}(undef, n, count),
        roots=Matrix{Complex{R}}(undef, n, count),
        modal_impedance=Matrix{Complex{R}}(undef, n, n),
        modal_admittance=Matrix{Complex{R}}(undef, n, n),
        impedance_sections=Vector{section}(undef, n),
        propagation_sections=Vector{section}(undef, n),
        propagation_count=zeros(Int, n),
        offdiagonal_ratio=zeros(R, 2), phase_error_degrees=zeros(R, n),
        transformation_imaginary_norm=zeros(R, n)))
end

function initialize_buffers(::Formula{:marti}, ::Val{:bode}, ::Type{T}, input, plan,
        buffers) where {T}
    return merge(buffers, (log_frequency=log10.(input.omega),
        decibels=20 .* log10.(abs.(input.data)), poles=T[], zeros=T[]))
end

function bode_asymptotes!(::Formula{:marti}, buffers, omega, tolerance)
    x, y = buffers.log_frequency, buffers.decibels
    slope = zero(eltype(x))
    intercept = first(y)
    origin = first(x)
    for k in 2:length(x)
        difference = y[k] - (intercept + slope * (x[k] - origin))
        if abs(difference) > tolerance
            corner = k - 1
            change = 20 * sign(difference)
            push!(change < 0 ? buffers.poles : buffers.zeros, omega[corner])
            intercept += slope * (x[corner] - origin)
            origin = x[corner]
            slope += change
        end
    end
    return buffers
end

function validate(samples::Tuple{<:AbstractVector, <:AbstractVector, <:Real}, ::Formula{:marti})
    omega, magnitude, tolerance = samples
    length(omega) == length(magnitude) && !isempty(omega) ||
        throw(ArgumentError("Bode samples must have equal nonzero lengths"))
    all(>(0), omega) && all(>(0), diff(omega)) && all(isfinite, omega) ||
        throw(ArgumentError("Bode frequencies must be finite, positive and increasing"))
    all(>(0), magnitude) && all(isfinite, magnitude) ||
        throw(ArgumentError("Bode magnitudes must be finite and positive"))
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("Bode tolerance must be finite and positive"))
    return samples
end

function bode_fit(formula::Formula{:marti}, omega, data, tolerance;
        delay=zero(eltype(omega)), threshold=zero(eltype(omega)), proper=false)
    validate((omega, max.(abs.(data), threshold), tolerance), formula)
    first_below = findfirst(a -> abs(a) < threshold, data)
    count = first_below === nothing ? length(omega) : first_below - 1
    support = (omega=omega[1:count], data=data[1:count])
    validate((support.omega, abs.(support.data), tolerance), formula)
    input = (; omega, data, threshold)
    asymptotes = initialize_buffers(formula, Val(:bode), eltype(omega), support, (;), (;))
    bode_asymptotes!(formula, asymptotes, support.omega, tolerance)
    # Complete the asymptotes above the band to give the required relative degree.
    while length(asymptotes.poles) < length(asymptotes.zeros) + proper
        push!(asymptotes.poles, isempty(asymptotes.poles) ? 10last(omega) :
            max(10last(omega), 1.25last(asymptotes.poles)))
    end
    if !proper
        while length(asymptotes.zeros) < length(asymptotes.poles)
            push!(asymptotes.zeros, isempty(asymptotes.zeros) ? 10last(omega) :
                max(10last(omega), 1.25last(asymptotes.zeros)))
        end
    end
    plan = (poles=asymptotes.poles, zeros=asymptotes.zeros, delay=delay, proper=proper)
    buffers = initialize_buffers(formula, Val(:factors), eltype(omega), input, plan, (;))
    return refine_poles!(formula, input, plan, buffers)
end

function initialize_buffers(::Formula{:marti}, ::Val{:factors}, ::Type{T}, input, plan,
        buffers) where {T}
    poles = copy(plan.poles)
    upper = 10maximum(poles; init=first(input.omega))
    # Distinct simple poles give nonzero denominators in the residue formula.
    for j in 2:length(poles)
        poles[j] = max(poles[j], 1.25poles[j - 1])
    end
    upper = max(upper, 1.25maximum(poles; init=first(input.omega)))
    zeros = T[j <= length(plan.zeros) ? plan.zeros[j] : (poles[j] + upper) / 2
        for j in 1:(length(poles) - plan.proper)]
    if plan.proper
        for j in eachindex(zeros)
            zeros[j] = clamp(zeros[j], 1.01poles[j], 0.99upper)
        end
    end
    gain = log(abs(first(input.data))) + sum(log, poles; init=zero(T)) -
        sum(log, zeros; init=zero(T))
    theta = T[gain; log.(zeros); log.(poles)]
    count, columns = length(input.omega), length(theta)
    evaluations = ntuple(2) do _
        (poles=similar(poles), zeros=similar(zeros),
            gain=Vector{T}(undef, count), residual=Vector{T}(undef, count),
            jacobian=Matrix{T}(undef, count, columns))
    end
    return merge(buffers, (theta=theta, trial_theta=similar(theta),
        current=first(evaluations), trial=last(evaluations),
        normal=Matrix{T}(undef, columns, columns), gradient=Vector{T}(undef, columns),
        weights=inv.(max.(abs.(input.data), input.threshold)),
        target=abs.(input.data), upper=log(upper), gap=log(T(1.2))))
end

function bode_magnitude!(::Formula{:marti}, theta, input, buffers, weights, target)
    nz = length(buffers.zeros)
    buffers.zeros .= exp.(@view theta[2:(nz + 1)])
    buffers.poles .= exp.(@view theta[(nz + 2):end])
    fill!(buffers.gain, first(theta) + sum(@view theta[2:(nz + 1)]) -
        sum(@view theta[(nz + 2):end]))
    buffers.jacobian[:, 1] .= 1
    for j in 2:length(theta), k in eachindex(input.omega)
        direction = j <= nz + 1 ? 1 : -1
        x = 2(log(input.omega[k]) - theta[j])
        buffers.gain[k] += direction * (max(x, 0) + log1p(exp(-abs(x)))) / 2
        # Derivatives use DC gain coordinates to avoid subtracting nearly equal terms.
        buffers.jacobian[k, j] = -direction *
            (x >= 0 ? inv(1 + exp(-x)) : exp(x) / (1 + exp(x)))
    end
    buffers.gain .= exp.(buffers.gain)
    @. buffers.residual = (buffers.gain - target) * weights
    for j in axes(buffers.jacobian, 2), k in eachindex(input.omega)
        buffers.jacobian[k, j] *= buffers.gain[k] * weights[k]
    end
    return nothing
end

function refine_poles!(formula::Formula{:marti}, input, plan, buffers)
    theta, trial_theta = buffers.theta, buffers.trial_theta
    current, trial = buffers.current, buffers.trial
    bode_magnitude!(formula, theta, input, current, buffers.weights, buffers.target)
    error = norm(current.residual)
    damping = eltype(theta)(1e-3)
    n, nz = length(current.poles), length(current.zeros)
    for iteration in 1:300
        mul!(buffers.normal, transpose(current.jacobian), current.jacobian)
        mul!(buffers.gradient, transpose(current.jacobian), current.residual)
        for j in eachindex(theta)
            buffers.normal[j, j] += damping * max(buffers.normal[j, j], eps(eltype(theta)))
        end
        step = buffers.normal \ buffers.gradient
        # Map the balanced DC gain step to the stored log k coordinate.
        step[1] += sum(@view step[(nz + 2):end]) - sum(@view step[2:(nz + 1)])
        norm(step) < 1e-10 * max(norm(theta), one(eltype(theta))) && break
        trial_theta .= theta .- step
        gain = trial_theta[1] + sum(@view trial_theta[2:(nz + 1)]) -
            sum(@view trial_theta[(nz + 2):end])
        if n > 0
            trial_theta[nz + 2] = min(trial_theta[nz + 2], buffers.upper - (n - 1) * buffers.gap)
            for j in 2:n
                index = nz + 1 + j
                trial_theta[index] = clamp(trial_theta[index], trial_theta[index - 1] + buffers.gap,
                    buffers.upper - (n - j) * buffers.gap)
            end
        end
        if plan.proper
            for j in 1:nz
                trial_theta[1 + j] = clamp(trial_theta[1 + j], trial_theta[nz + 1 + j], buffers.upper)
            end
            # The low-frequency propagation gain cannot exceed one.
            gain = min(gain, zero(eltype(theta)))
        end
        trial_theta[1] = gain + sum(@view trial_theta[(nz + 2):end]) -
            sum(@view trial_theta[2:(nz + 1)])
        bode_magnitude!(formula, trial_theta, input, trial, buffers.weights, buffers.target)
        trial_error = norm(trial.residual)
        if isfinite(trial_error) && trial_error < error &&
                all(p -> isfinite(p) && p > 0, trial.poles) &&
                all(z -> isfinite(z) && z > 0, trial.zeros) && all(isfinite, trial.jacobian)
            improvement = (error - trial_error) / error
            theta .= trial_theta
            current, trial = trial, current
            error = trial_error
            damping = max(damping / 3, eps(eltype(theta)))
            improvement < 1e-12 && break
        else
            damping *= 10
        end
    end
    # Evaluate the factors directly in high precision.
    # Poles remain separate, with a ratio of at least 1.2 between neighbours.
    residues = setprecision(BigFloat, 512) do
        poles, zeros = BigFloat.(current.poles), BigFloat.(current.zeros)
        gain = exp(BigFloat(first(theta)))
        [eltype(theta)(gain * prod(z - poles[j] for z in zeros; init=one(BigFloat)) /
            prod(poles[k] - poles[j] for k in eachindex(poles) if k != j;
                init=one(BigFloat))) for j in eachindex(poles)]
    end
    return (constant=plan.proper ? zero(eltype(theta)) : exp(first(theta)),
        poles=-current.poles, residues=residues, delay=plan.delay)
end

function initialize_buffers(::Formula{:marti}, ::Val{:delay}, ::Type{T}, input, plan,
        buffers) where {T}
    return merge(buffers, (phase=zeros(T, length(input.omega)),
        log_frequency=log.(input.omega), log_gain=log.(abs.(input.propagation))))
end

function refine_delay(formula::Formula{:marti}, omega, propagation, root, length_m;
        section=nothing)
    estimate = length_m * imag(root) / last(omega)
    buffers = initialize_buffers(formula, Val(:delay), eltype(omega),
        (; omega, propagation), (;), (;))
    if section === nothing
        # Bode gain-phase relation on the log-frequency grid, with constant tails.
        x, gain, phase = buffers.log_frequency, buffers.log_gain, buffers.phase
        for k in eachindex(x)
            integral = zero(eltype(phase))
            for j in 1:(length(x) - 1)
                left = j == k ? (gain[j + 1] - gain[j]) / (x[j + 1] - x[j]) :
                    (gain[j] - gain[k]) / sinh(x[j] - x[k])
                right = j + 1 == k ? (gain[j + 1] - gain[j]) / (x[j + 1] - x[j]) :
                    (gain[j + 1] - gain[k]) / sinh(x[j + 1] - x[k])
                integral += (left + right) * (x[j + 1] - x[j]) / 2
            end
            k > 1 && (integral += (first(gain) - gain[k]) *
                log(tanh((x[k] - first(x)) / 2)))
            k < length(x) && (integral -= (last(gain) - gain[k]) *
                log(tanh((last(x) - x[k]) / 2)))
            phase[k] = integral / π
        end
    else
        previous = zero(estimate)
        for k in eachindex(omega)
            value = section.constant + sum(r / (im * omega[k] - p)
                for (p, r) in zip(section.poles, section.residues); init=zero(complex(estimate)))
            phase = angle(value)
            buffers.phase[k] = k == 1 ? phase : buffers.phase[k - 1] +
                rem2pi(phase - previous, RoundNearest)
            previous = phase
        end
    end
    best, error = estimate, eltype(omega)(Inf)
    for j in 1:300
        # Include the initial delay among the 300 samples, also for a lossless line.
        delay = j == 150 ? estimate : estimate * (0.9 + 0.2 * (j - 1) / 299)
        previous, unwrapped, squared = zero(estimate), zero(estimate), zero(estimate)
        for k in eachindex(omega)
            current = angle(propagation[k] * cis(omega[k] * delay))
            unwrapped = k == 1 ? current : unwrapped + rem2pi(current - previous, RoundNearest)
            squared += (unwrapped - buffers.phase[k])^2
            previous = current
        end
        if squared < error
            best, error = delay, squared
        end
    end
    return best
end

function modal_samples!(::Formula{:marti}, workspace, transformation, voltage)
    buffers = workspace.buffers
    parameters = transform(PhaseDomain, workspace.input.segment.parameters)
    impedance_factor = lu(voltage)
    admittance_factor = lu(transformation)
    zm, ym = buffers.modal_impedance, buffers.modal_admittance
    for k in eachindex(workspace.input.frequencies)
        mul!(zm, @view(parameters.Z.values[:, :, k]), transformation)
        ldiv!(impedance_factor, zm)
        mul!(ym, @view(parameters.Y.values[:, :, k]), voltage)
        ldiv!(admittance_factor, ym)
        for (quantity, coefficients) in enumerate((zm, ym))
            diagonal, offdiagonal = zero(real(zm[1, 1])), zero(real(zm[1, 1]))
            for column in axes(coefficients, 2), row in axes(coefficients, 1)
                value = abs(coefficients[row, column])
                if row == column
                    diagonal = max(diagonal, value)
                else
                    offdiagonal = max(offdiagonal, value)
                end
            end
            ratio = iszero(offdiagonal) ? zero(diagonal) : offdiagonal / diagonal
            buffers.offdiagonal_ratio[quantity] = max(
                buffers.offdiagonal_ratio[quantity], ratio)
        end
        for mode in axes(zm, 1)
            root = sqrt(zm[mode, mode] * ym[mode, mode])
            real(root) < 0 && (root = -root)
            buffers.roots[mode, k] = root
            buffers.impedance[mode, k] = sqrt(zm[mode, mode] / ym[mode, mode])
            buffers.propagation[mode, k] = exp(-root * workspace.input.length)
        end
    end
    return workspace
end

function initialize_buffers(::Formula{:marti}, ::Val{:errors}, ::Type{T}, input, plan,
        buffers) where {T}
    n, count = plan.modes, length(input.frequencies)
    return merge(buffers, (fitted=Vector{Complex{T}}(undef, count),
        data=Vector{Complex{T}}(undef, count), relative=Vector{T}(undef, count),
        rms=zeros(T, n, n), max=zeros(T, n, n), f_max=zeros(T, n, n)))
end

function modal_fitting_error!(formula::Formula{:marti}, workspace)
    buffers = workspace.buffers
    n = size(buffers.impedance, 1)
    records = map((buffers.impedance_sections, buffers.propagation_sections),
            (buffers.impedance, buffers.propagation), (true, false)) do sections, samples, impedance
        errors = initialize_buffers(formula, Val(:errors), eltype(buffers.omega),
            workspace.input, (modes=n,), (;))
        fitted, data, relative = errors.fitted, errors.data, errors.relative
        for mode in 1:n
            section = sections[mode]
            for k in eachindex(buffers.omega)
                w = buffers.omega[k]
                fitted[k] = (section.constant + sum(r / (im * w - p)
                    for (p, r) in zip(section.poles, section.residues); init=zero(w))) *
                    cis(-w * section.delay)
                data[k] = samples[mode, k]
            end
            if impedance
                fitted .= inv.(fitted)
                data .= inv.(data)
            else
                count = buffers.propagation_count[mode]
                buffers.phase_error_degrees[mode] = maximum(
                    abs.(rad2deg.(angle.(fitted[1:count] ./ data[1:count]))))
            end
            relative .= abs.(abs.(fitted) .- abs.(data)) ./ maximum(abs, data)
            errors.rms[mode, mode] = sqrt(sum(abs2, relative) / length(relative))
            errors.max[mode, mode], index = findmax(relative)
            errors.f_max[mode, mode] = workspace.input.frequencies[index]
        end
        return (rms=errors.rms, max=errors.max, f_max=errors.f_max)
    end
    return (admittance=records[1], propagation=records[2])
end

function fit(formula::Formula{:marti}, workspace::LineModelWorkspace,
        parameters::NamedTuple, options::FormulationOptions)
    validate(options, formula)
    frequency = options.data.transformation.frequency
    frequencies = workspace.input.frequencies
    k = frequency === nothing ? length(frequencies) : argmin(abs.(frequencies .- frequency))
    current = copy(Ti(workspace.input.segment)[:, :, k:k])
    paired_voltage = copy(Tv(workspace.input.segment)[:, :, k:k])
    orient_modes!(paired_voltage, current, true)
    buffers = workspace.buffers
    for mode in axes(current, 2)
        column = @view current[:, mode, 1]
        buffers.transformation_imaginary_norm[mode] = norm(imag.(column)) / norm(column)
    end
    transformation = real.(current[:, :, 1])
    for column in eachcol(transformation)
        column ./= norm(column)
    end
    voltage = inv(transpose(transformation))
    modal_samples!(formula, workspace, transformation, voltage)
    for mode in axes(transformation, 2)
        impedance = @view buffers.impedance[mode, :]
        first_below = findfirst(a -> abs(a) < options.data.bode.threshold,
            @view buffers.propagation[mode, :])
        count = first_below === nothing ? length(frequencies) : first_below - 1
        buffers.propagation_count[mode] = count
        omega = @view buffers.omega[1:count]
        propagation = @view buffers.propagation[mode, 1:count]
        validate((omega, abs.(propagation), options.data.bode.tolerance_db), formula)
        delay = refine_delay(formula, omega, propagation,
            buffers.roots[mode, count], workspace.input.length)
        buffers.impedance_sections[mode] =
            bode_fit(formula, buffers.omega, impedance, options.data.bode.tolerance_db)
        lossless = all(iszero, real.(@view buffers.roots[mode, :]))
        section = bode_fit(formula, buffers.omega,
            buffers.propagation[mode, :] .* cis.(buffers.omega .* delay),
            options.data.bode.tolerance_db; delay,
            threshold=options.data.bode.threshold, proper=!lossless)
        delay = refine_delay(formula, omega, propagation,
            buffers.roots[mode, count], workspace.input.length; section)
        buffers.propagation_sections[mode] = merge(section, (; delay))
    end
    model = ModalModel(transformation, frequencies[k], buffers.impedance_sections,
        buffers.propagation_sections)
    errors = modal_fitting_error!(formula, workspace)
    passive = (checked=false, passive=nothing, bands=Tuple{eltype(frequencies), eltype(frequencies)}[],
        worst=nothing, f_worst=nothing, perturbation=nothing)
    diagnostics = (marti=(offdiagonal_ratio=buffers.offdiagonal_ratio,
        phase_error_degrees=buffers.phase_error_degrees,
        transformation_imaginary_norm=buffers.transformation_imaginary_norm,
        propagation_f_max=frequencies[buffers.propagation_count]),)
    return model, errors, passive, (targets=nothing, achieved=nothing), diagnostics
end

:marti
