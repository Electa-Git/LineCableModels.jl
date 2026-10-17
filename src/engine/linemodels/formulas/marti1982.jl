"""
$(TYPEDSIGNATURES)

**Identification.** J. Marti constant real modal transformation (1982).

**Expression.** Modal characteristic impedance and propagation magnitudes have
Bode asymptotes with slopes in multiples of 20 dB per decade. Their corners
supply starting real poles and zeros. Variable projection refines negative real
poles by least squares against complex data. Each projected solve identifies
the constant and real residues. Propagation retains a separate travel delay.

The options select the transformation frequency [Hz] and magnitude tolerance [dB].
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
        offdiagonal_ratio=zeros(R, 2), phase_error_degrees=zeros(R, n)))
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

function bode_fit(formula::Formula{:marti}, omega, data, tolerance; delay=zero(eltype(omega)))
    validate((omega, abs.(data), tolerance), formula)
    input = (; omega, data)
    asymptotes = initialize_buffers(formula, Val(:bode), eltype(omega), input, (;), (;))
    bode_asymptotes!(formula, asymptotes, omega, tolerance)
    plan = (poles=-asymptotes.poles, delay=delay)
    buffers = initialize_buffers(formula, Val(:residues), eltype(omega), input, plan, (;))
    return refine_poles!(formula, input, plan, buffers)
end

function initialize_buffers(::Formula{:marti}, ::Val{:residues}, ::Type{T}, input, plan,
        buffers) where {T}
    count, columns = length(input.omega), length(plan.poles) + 1
    return merge(buffers, (system=Matrix{T}(undef, 2count, columns),
        target=Vector{T}(undef, 2count), scales=Vector{T}(undef, columns),
        residual=Vector{T}(undef, 2count), projected=Vector{T}(undef, 2count),
        jacobian=Matrix{T}(undef, 2count, columns - 1),
        normal=Matrix{T}(undef, columns - 1, columns - 1),
        gradient=Vector{T}(undef, columns - 1),
        theta=log.(-plan.poles), trial_theta=Vector{T}(undef, columns - 1),
        trial_poles=Vector{T}(undef, columns - 1)))
end

function bode_residues(::Formula{:marti}, input, plan, buffers)
    omega, data = input.omega, input.data
    poles = plan.poles
    count = length(omega)
    # Stack the complex data with relative weights and solve for a real constant
    # and real residues together. Scaling balances the constant and pole columns.
    for k in eachindex(omega)
        weight = inv(abs(data[k]))
        target = data[k] * weight
        buffers.target[k] = real(target)
        buffers.target[count + k] = imag(target)
        buffers.system[k, 1] = weight
        buffers.system[count + k, 1] = zero(weight)
        for j in eachindex(poles)
            value = weight / (im * omega[k] - poles[j])
            buffers.system[k, j + 1] = real(value)
            buffers.system[count + k, j + 1] = imag(value)
        end
    end
    for j in axes(buffers.system, 2)
        buffers.scales[j] = norm(@view buffers.system[:, j])
        @views buffers.system[:, j] ./= buffers.scales[j]
    end
    scaled = buffers.system \ buffers.target
    mul!(buffers.projected, buffers.system, scaled)
    buffers.projected .-= buffers.target
    coefficients = scaled ./ buffers.scales
    return (constant=first(coefficients), poles=copy(poles), residues=coefficients[2:end],
        delay=plan.delay)
end

function refine_poles!(formula::Formula{:marti}, input, plan, buffers)
    section = bode_residues(formula, input, plan, buffers)
    isempty(plan.poles) && return section
    theta, trial_theta = buffers.theta, buffers.trial_theta
    normal, gradient = buffers.normal, buffers.gradient
    residual, projected, jacobian = buffers.residual, buffers.projected, buffers.jacobian
    residual .= projected
    error = norm(residual)
    damping = eltype(theta)(1e-3)
    trial = (poles=buffers.trial_poles, delay=plan.delay)
    for iteration in 1:50
        for j in eachindex(theta)
            increment = cbrt(eps(eltype(theta))) * max(one(eltype(theta)), abs(theta[j]))
            trial_theta .= theta
            trial_theta[j] += increment
            trial.poles .= -exp.(trial_theta)
            bode_residues(formula, input, trial, buffers)
            @views jacobian[:, j] .= (projected .- residual) ./ increment
        end
        mul!(normal, transpose(jacobian), jacobian)
        mul!(gradient, transpose(jacobian), residual)
        for j in eachindex(theta)
            normal[j, j] += damping * max(normal[j, j], eps(eltype(theta)))
        end
        trial_theta .= theta .- (normal \ gradient)
        trial.poles .= -exp.(trial_theta)
        if !all(p -> isfinite(p) && p < 0, trial.poles)
            damping *= 10
            continue
        end
        candidate = bode_residues(formula, input, trial, buffers)
        candidate_error = norm(projected)
        if candidate_error < error
            improvement = (error - candidate_error) / error
            theta .= trial_theta
            residual .= projected
            section = candidate
            error = candidate_error
            damping = max(damping / 3, eps(eltype(theta)))
            improvement < 1e-8 && break
        else
            damping *= 10
        end
    end
    return section
end

function minimum_phase!(::Formula{:marti}, buffers)
    # Bode gain-phase relation after x = log(omega):
    # phase(x) = integral (log|F(exp(u))| - log|F(exp(x))|) / sinh(u-x) du / pi.
    # Subtraction removes the principal-value singularity. Endpoint gains extend
    # constantly outside the measured band, whose two tails integrate analytically.
    x = buffers.log_frequency
    gain = buffers.log_gain
    phase = buffers.phase
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
    return phase
end

function initialize_buffers(::Formula{:marti}, ::Val{:delay}, ::Type{T}, input, plan,
        buffers) where {T}
    return merge(buffers, (phase=zeros(T, length(input.omega)),
        log_frequency=log.(input.omega), log_gain=log.(abs.(input.propagation))))
end

function refine_delay(formula::Formula{:marti}, omega, propagation, root, length_m)
    estimate = length_m * imag(root) / last(omega)
    buffers = initialize_buffers(formula, Val(:delay), eltype(omega),
        (; omega, propagation), (;), (;))
    minimum_phase!(formula, buffers)
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
    transformation = real.(Ti(workspace.input.segment)[:, :, k])
    for column in eachcol(transformation)
        column ./= norm(column)
    end
    voltage = inv(transpose(transformation))
    modal_samples!(formula, workspace, transformation, voltage)
    buffers = workspace.buffers
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
        buffers.impedance_sections[mode] = bode_fit(formula, buffers.omega,
            impedance, options.data.bode.tolerance_db)
        buffers.propagation_sections[mode] = bode_fit(formula, omega,
            propagation .* cis.(omega .* delay), options.data.bode.tolerance_db; delay)
    end
    model = ModalModel(transformation, frequencies[k], buffers.impedance_sections,
        buffers.propagation_sections)
    errors = modal_fitting_error!(formula, workspace)
    passive = (checked=false, passive=nothing, bands=Tuple{eltype(frequencies), eltype(frequencies)}[],
        worst=nothing, f_worst=nothing, perturbation=nothing)
    diagnostics = (marti=(offdiagonal_ratio=buffers.offdiagonal_ratio,
        phase_error_degrees=buffers.phase_error_degrees,
        propagation_f_max=frequencies[buffers.propagation_count]),)
    return model, errors, passive, (targets=nothing, achieved=nothing), diagnostics
end

:marti
