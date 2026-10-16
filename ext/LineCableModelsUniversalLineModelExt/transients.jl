solid_ground(assignment::TerminalImpedance) = iszero(assignment.R) && iszero(assignment.L)
solid_ground(::Waveform) = false

series_nodes(assignment::TerminalImpedance) =
    !iszero(assignment.R) && !iszero(assignment.L) ? 1 : 0
series_nodes(::Waveform) = 0

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
    throw(ArgumentError("Matched() needs a fitted characteristic-admittance branch, which the UniversalLineModel network does not provide"))
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
    extra0 = sum(series_nodes, problem.sending; init=0)
    extraL = sum(series_nodes, problem.receiving; init=0)
    net = Network(countL + extra0 + extraL)
    (; interp, enforce_sign, initial_solve) = selected.options.data
    add_line!(net, model, nodes0, nodesL; interp, enforce_sign)
    extra = countL
    for (assignments, nodes) in ((problem.sending, nodes0), (problem.receiving, nodesL))
        for assignment in assignments
            index = something(findfirst(==(assignment.terminal), coordinates))
            extra = connect!(net, assignment, nodes[index], extra)
        end
    end
    res = simulate(net, problem.window.step, problem.window.tmax; initial_solve)
    nt = length(res.t)
    voltages = zeros(eltype(res.v), 2, n, nt)
    currents = similar(voltages)
    for (side, nodes, values) in ((1, nodes0, res.lines[1].i0), (2, nodesL, res.lines[1].iL))
        for i in 1:n, k in 1:nt
            node = nodes[i]
            voltages[side, i, k] = iszero(node) ? 0 : res.v[k, node]
            currents[side, i, k] = values[k, i]
        end
    end
    retained = (td=(; interp, enforce_sign, initial_solve, steps=nt), line=details(problem.line))
    return res.t, voltages, currents, coordinates, retained
end
