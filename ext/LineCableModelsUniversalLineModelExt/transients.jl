solid_ground(assignment::TerminalImpedance) = iszero(assignment.R) && iszero(assignment.L)
solid_ground(::Waveform) = false

series_nodes(assignment::TerminalImpedance) =
    !iszero(assignment.R) && !iszero(assignment.L) ? 1 : 0
series_nodes(::Waveform) = 0
series_nodes(::Matched) = 0
series_nodes(assignments::Tuple) = sum(series_nodes, assignments; init=0)

function terminal_nodes(::TransientAnalysis.Formula{:td}, assignments::Tuple, coordinates, offset)
    nodes = ones(Int, length(coordinates))
    for assignment in assignments
        if solid_ground(assignment)
            index = something(findfirst(==(assignment.terminal), coordinates))
            nodes[index] = 0
        end
    end
    for i in eachindex(nodes)
        if !iszero(nodes[i])
            offset += 1
            nodes[i] = offset
        end
    end
    return nodes, offset
end

function terminal_nodes(::TransientAnalysis.Formula{:td}, ::Matched, coordinates, offset)
    count = offset + length(coordinates)
    return collect((offset + 1):count), count
end

function connect!(net::Network, assignments::Tuple, nodes, coordinates, model, extra)
    for assignment in assignments
        index = something(findfirst(==(assignment.terminal), coordinates))
        extra = connect!(net, assignment, nodes[index], extra)
    end
    return extra, nothing
end

function connect!(net::Network, ::Matched, nodes, coordinates, model, extra)
    add_rational_admittance!(net, nodes, model.pY, model.RY, model.DY)
    return extra, (admittance=:fitted_characteristic, poles=length(model.pY))
end

function connect!(net::Network, assignment::TerminalImpedance, node, extra)
    if iszero(assignment.R)
        iszero(assignment.L) || add_inductor!(net, node, 0, assignment.L)
    elseif iszero(assignment.L)
        add_resistor!(net, node, 0, assignment.R)
    else
        extra += 1
        add_resistor!(net, node, extra, assignment.R)
        add_inductor!(net, extra, 0, assignment.L)
    end
    return extra
end

function connect!(net::Network, source::Step, node, extra)
    add_vsource!(net, node, StepSource(source.amplitude, source.t0))
    return extra
end

function connect!(net::Network, source::Cosine, node, extra)
    add_vsource!(net, node, CosineSource(source.amplitude, 2π * source.frequency, source.phase))
    return extra
end

function connect!(net::Network, source::DoubleExponential, node, extra)
    add_vsource!(net, node, DoubleExpSource(source.amplitude, source.alpha, source.beta))
    return extra
end

function compute(selected::TransientAnalysis.Formula{:td},
        problem::TransientAnalysisProblem{<:LineModel{<:WidebandModel}}, execution)
    model = ULMModel(problem.line)
    n = model.n
    coordinates = problem.line.coordinates
    nodes0, count0 = terminal_nodes(selected, problem.sending, coordinates, 0)
    nodesL, countL = terminal_nodes(selected, problem.receiving, coordinates, count0)
    extra0 = series_nodes(problem.sending)
    extraL = series_nodes(problem.receiving)
    net = Network(countL + extra0 + extraL)
    (; interp, enforce_sign, initial_solve) = selected.options.data
    add_line!(net, model, nodes0, nodesL; interp, enforce_sign)
    extra, sending = connect!(net, problem.sending, nodes0, coordinates, model, countL)
    _, receiving = connect!(net, problem.receiving, nodesL, coordinates, model, extra)
    res = simulate(net, problem.window.step, problem.window.tmax; initial_solve)
    nt = length(res.t)
    voltages = zeros(eltype(res.v), 2, n, nt)
    currents = similar(voltages)
    for (side, volts, amps) in ((1, res.lines[1].v0, res.lines[1].i0),
            (2, res.lines[1].vL, res.lines[1].iL))
        for i in 1:n, k in 1:nt
            voltages[side, i, k] = volts[k, i]
            currents[side, i, k] = amps[k, i]
        end
    end
    retained = (td=(; interp, enforce_sign, initial_solve, steps=nt,
        terminations=(; sending, receiving)), line=details(problem.line))
    return res.t, voltages, currents, coordinates, retained
end
