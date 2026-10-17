"""
$(TYPEDSIGNATURES)

Select time-domain simulation of a fitted line model.
"""
description(::Type{<:Formula{:td}}; compact::Bool=false) =
    compact ? "Time domain" : "Time-domain simulation of a fitted line model"

function formulation_options(::Expression{<:Formula{:td}, typeof(respond!)})
    return FormulationOptions((interp=:single, enforce_sign=false, initial_solve=false))
end

function formulation_options(expression::Expression{<:Formula{:td}, typeof(respond!)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :td formulation options"))
    controls = merge(defaults, supplied.data)
    controls.interp::Symbol
    controls.enforce_sign::Bool
    controls.initial_solve::Bool
    controls.interp in (:single, :double) ||
        throw(ArgumentError("interp must be :single or :double"))
    return FormulationOptions(controls)
end

validate(problem::TransientAnalysisProblem{<:LineModel}, ::Formula{:td}) = problem

function validate(problem::TransientAnalysisProblem{<:LineModel{<:ModalModel}},
        selected::Formula{:td})
    controls = selected.options.data
    controls.interp === :single ||
        throw(ArgumentError("interp must be :single for a ModalModel"))
    !controls.enforce_sign ||
        throw(ArgumentError("enforce_sign must be false for a ModalModel"))
    model = problem.line.model
    n = length(problem.line.coordinates)
    size(model.transformation) == (n, n) &&
        length(model.impedance) == length(model.propagation) == n && n > 0 ||
        throw(ArgumentError("modal sections and transformation must match the terminals"))
    all(isfinite, model.transformation) ||
        throw(ArgumentError("modal transformation must be finite"))
    for sections in (model.impedance, model.propagation), section in sections
        length(section.poles) == length(section.residues) &&
            all(isfinite, section.poles) && all(isfinite, section.residues) &&
            isfinite(section.constant) && isfinite(section.delay) ||
            throw(ArgumentError("modal sections must have finite coefficients and paired poles and residues"))
    end
    all(section -> iszero(section.delay), model.impedance) ||
        throw(ArgumentError("modal characteristic impedance must have zero delay"))
    delay = minimum(section -> section.delay, model.propagation)
    problem.window.step <= delay || throw(ArgumentError(
        "time step ($(problem.window.step) s) larger than the smallest delay ($delay s)"))
    return problem
end

function TransientWorkspace(problem::TransientAnalysisProblem{<:LineModel{<:ModalModel}},
        selected::Formula{:td})
    validate(problem, selected)
    n = length(problem.line.coordinates)
    input = (model=problem.line.model, sending=problem.sending, receiving=problem.receiving,
        window=problem.window)
    plan = (n=n, nodes=2n, coordinates=problem.line.coordinates,
        nt=round(Int, problem.window.tmax / problem.window.step) + 1,
        unknowns=2n + unknown_count(problem.sending) + unknown_count(problem.receiving))
    buffers = initialize_buffers(selected, eltype(input.model.transformation), input, plan, (;))
    return TransientWorkspace(input, plan, buffers)
end

function initialize_buffers(::Formula{:td}, ::Val{:impedance}, ::Type{T}, section,
        plan, buffers) where {T}
    h = plan.step
    # For x' = p*x + r*i, trapezoidal integration gives
    # (1-p*h/2)*x_n = (1+p*h/2)*x_prev + r*h/2*(i_n+i_prev).
    # a=(2+p*h)/(2-p*h), b=r*h/(2-p*h), x_n=a*x_prev+b*(i_n+i_prev).
    # Z*i = R_eq*i_n + H_n, R_eq=k0+sum(b), H_n=sum(a*x_prev+b*i_prev).
    # Store H per section: H_next=a*H_n+b*(1+a)*i_n.
    a = @. (2 + section.poles * h) / (2 - section.poles * h)
    b = @. section.residues * h / (2 - section.poles * h)
    return merge(buffers, (a=a, b=b, c=b .* (1 .+ a),
        resistance=section.constant + sum(b), history=zeros(T, 2, length(a)),
        load_history=zeros(T, 2, length(a))))
end

function convolution_coefficients(::Formula{:td}, p, h)
    z = p * h
    a = exp(z)
    # Exact integration of exp(p*(h-u)) against the linear input on 0 <= u <= h.
    if abs(z) <= 0.1
        cn = h * evalpoly(z, (1/2, 1/6, 1/24, 1/120, 1/720,
            1/5040, 1/40320, 1/362880, 1/3628800))
        cp = h * evalpoly(z, (1/2, 1/3, 1/8, 1/30, 1/144,
            1/840, 1/5760, 1/45360, 1/403200))
    else
        difference = expm1(z)
        cn = h * (difference / z - 1) / z
        cp = h * (a - difference / z) / z
    end
    return a, cn, cp
end

function initialize_buffers(selected::Formula{:td}, ::Val{:propagation}, ::Type{T}, section,
        plan, buffers) where {T}
    coefficients = [convolution_coefficients(selected, p, plan.step) for p in section.poles]
    return merge(buffers, (a=T[first(value) for value in coefficients],
        cn=T[section.residues[j] * coefficients[j][2] for j in eachindex(coefficients)],
        cp=T[section.residues[j] * coefficients[j][3] for j in eachindex(coefficients)],
        history=zeros(T, 2, length(coefficients))))
end

function initialize_buffers(selected::Formula{:td}, ::Type{T}, input::NamedTuple,
        plan::NamedTuple, buffers) where {T}
    n, nt = plan.n, plan.nt
    h = input.window.step
    model = input.model
    impedance = [initialize_buffers(selected, Val(:impedance), T, section, (step=h,), (;))
        for section in model.impedance]
    propagation = [initialize_buffers(selected, Val(:propagation), T, section, (step=h,), (;))
        for section in model.propagation]
    delay_steps = T[section.delay / h for section in model.propagation]
    for mode in eachindex(delay_steps)
        nearest = round(delay_steps[mode])
        abs(delay_steps[mode] - nearest) <= 8eps(delay_steps[mode]) &&
            (delay_steps[mode] = nearest)
    end
    conductance = model.transformation *
        Diagonal(inv.([section.resistance for section in impedance])) * transpose(model.transformation)
    arrays = merge(buffers, (M=zeros(T, plan.unknowns, plan.unknowns),
        rhs=zeros(T, plan.unknowns), excitation=zeros(T, plan.unknowns, nt),
        conductance=conductance, impedance=impedance, propagation=propagation,
        lags=floor.(Int, delay_steps), fractions=delay_steps .- floor.(delay_steps),
        matched=falses(2), terminal_a=zeros(T, plan.nodes), terminal_b=zeros(T, plan.nodes),
        terminal_history=zeros(T, plan.nodes), incident=zeros(T, 2, n),
        injection=zeros(T, 2, n), delayed=zeros(T, 2, n),
        outgoing=zeros(T, 2, n, nt),
        t=T[(k - 1) * h for k in 1:nt],
        voltages=zeros(T, 2, n, nt), currents=zeros(T, 2, n, nt)))
    workspace = TransientWorkspace(input, plan, arrays)
    for side in 0:1, j in 1:n, i in 1:n
        arrays.M[side*n+i, side*n+j] = conductance[i, j]
    end
    unknown = stamp!(selected, input.sending, workspace, 0, plan.nodes)
    stamp!(selected, input.receiving, workspace, n, unknown)
    return merge(arrays, (factor=lu(arrays.M),))
end

function stamp!(selected::Formula{:td}, assignments::Tuple, workspace, offset, unknown)
    for assignment in assignments
        index = something(findfirst(==(assignment.terminal), workspace.plan.coordinates))
        unknown = stamp!(selected, assignment, workspace, offset + index, unknown)
    end
    return unknown
end

function stamp!(::Formula{:td}, ::Matched, workspace, offset, unknown)
    arrays, n = workspace.buffers, workspace.plan.n
    arrays.matched[offset ÷ n + 1] = true
    for column in 1:n, row in 1:n
        arrays.M[offset+row, offset+column] += arrays.conductance[row, column]
    end
    return unknown
end

function stamp!(::Formula{:td}, assignment::TerminalImpedance, workspace, node, unknown)
    arrays = workspace.buffers
    if iszero(assignment.R) && iszero(assignment.L)
        unknown += 1
        arrays.M[unknown, node] = arrays.M[node, unknown] = 1
    else
        inductive = 2assignment.L / workspace.input.window.step
        g = inv(assignment.R + inductive)
        arrays.M[node, node] += g
        arrays.terminal_a[node] = (inductive - assignment.R) * g
        arrays.terminal_b[node] = (1 + arrays.terminal_a[node]) * g
    end
    return unknown
end

function stamp!(::Formula{:td}, source::Waveform, workspace, node, unknown)
    arrays = workspace.buffers
    unknown += 1
    arrays.M[node, unknown] = arrays.M[unknown, node] = 1
    for k in eachindex(arrays.t)
        arrays.excitation[unknown, k] = source(arrays.t[k])
    end
    return unknown
end

function delayed_wave(::Formula{:td}, arrays, far, mode, k)
    right = k - arrays.lags[mode]
    fraction = arrays.fractions[mode]
    (right < 1 || (right == 1 && fraction > 0)) && return zero(eltype(arrays.outgoing))
    previous = right > 1 ? arrays.outgoing[far, mode, right-1] : zero(eltype(arrays.outgoing))
    return (1 - fraction) * arrays.outgoing[far, mode, right] + fraction * previous
end

function modal_history!(selected::Formula{:td}, workspace, k)
    arrays, n = workspace.buffers, workspace.plan.n
    transformation = workspace.input.model.transformation
    for mode in 1:n, side in 1:2
        section = arrays.propagation[mode]
        u = delayed_wave(selected, arrays, 3 - side, mode, k)
        incident = workspace.input.model.propagation[mode].constant * u
        for j in eachindex(section.a)
            section.history[side, j] = section.a[j] * section.history[side, j] +
                section.cn[j] * u + section.cp[j] * arrays.delayed[side, mode]
            incident += section.history[side, j]
        end
        arrays.delayed[side, mode] = u
        arrays.incident[side, mode] = incident
        impedance = arrays.impedance[mode]
        line_history, load_history = zero(incident), zero(incident)
        for j in eachindex(impedance.a)
            line_history += impedance.history[side, j]
            load_history += impedance.load_history[side, j]
        end
        arrays.injection[side, mode] = (incident + line_history) / impedance.resistance
        total = arrays.injection[side, mode] + load_history / impedance.resistance
        for phase in 1:n
            arrays.rhs[(side-1)*n+phase] += transformation[phase, mode] * total
        end
    end
    return nothing
end

function advance_modal!(selected::Formula{:td}, workspace::TransientWorkspace, indices)
    arrays, plan = workspace.buffers, workspace.plan
    transformation = workspace.input.model.transformation
    # T maps currents as i=T*i_m and voltages as v_m=transpose(T)*v.
    # The phase conductance is G=T*inv(R_eq)*transpose(T).
    for k in indices
        for row in eachindex(arrays.rhs)
            arrays.rhs[row] = arrays.excitation[row, k]
        end
        for node in 1:plan.nodes
            arrays.rhs[node] -= arrays.terminal_history[node]
        end
        modal_history!(selected, workspace, k)
        ldiv!(arrays.factor, arrays.rhs)
        for phase in 1:plan.n, side in 1:2
            node = (side - 1) * plan.n + phase
            v = arrays.rhs[node]
            arrays.voltages[side, phase, k] = v
            arrays.currents[side, phase, k] = 0
            arrays.terminal_history[node] = arrays.terminal_a[node] *
                arrays.terminal_history[node] + arrays.terminal_b[node] * v
        end
        for mode in 1:plan.n, side in 1:2
            vm = zero(eltype(arrays.rhs))
            for phase in 1:plan.n
                vm += transformation[phase, mode] * arrays.voltages[side, phase, k]
            end
            impedance = arrays.impedance[mode]
            modal_current = vm / impedance.resistance - arrays.injection[side, mode]
            load = zero(vm)
            for j in eachindex(impedance.a)
                load += impedance.load_history[side, j]
            end
            load = (vm - load) / impedance.resistance
            for j in eachindex(impedance.a)
                impedance.history[side, j] = impedance.a[j] * impedance.history[side, j] +
                    impedance.c[j] * modal_current
                if arrays.matched[side]
                    impedance.load_history[side, j] = impedance.a[j] *
                        impedance.load_history[side, j] + impedance.c[j] * load
                end
            end
            arrays.outgoing[side, mode, k] = 2vm - arrays.incident[side, mode]
            for phase in 1:plan.n
                arrays.currents[side, phase, k] += transformation[phase, mode] * modal_current
            end
        end
    end
    return nothing
end

function respond!(selected::Formula{:td}, workspace::TransientWorkspace, parameters, options)
    start = options.data.initial_solve ? 1 : 2
    advance_modal!(selected, workspace, start:workspace.plan.nt)
    return workspace.buffers.t, workspace.buffers.voltages, workspace.buffers.currents
end

:td
