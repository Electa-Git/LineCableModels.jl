"""
$(TYPEDSIGNATURES)

Select the numerical Laplace transform of the line two-port.
"""
description(::Type{<:Formula{:nlt}}; compact::Bool=false) =
    compact ? "NLT" : "Numerical Laplace transform of the line two-port"

function formulation_options(::Expression{<:Formula{:nlt}, typeof(respond!)})
    return FormulationOptions((line=Formulation(), N=4096, Tfactor=1.25, taper=:hanning))
end

function formulation_options(expression::Expression{<:Formula{:nlt}, typeof(respond!)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :nlt formulation options"))
    controls = merge(defaults, supplied.data)
    controls.line::LineParametersFormulation
    controls.N::Int
    controls.N >= 2 || throw(DomainError(controls.N, "N must be at least 2"))
    isfinite(controls.Tfactor) && controls.Tfactor > 1 ||
        throw(DomainError(controls.Tfactor, "Tfactor must be finite and greater than 1"))
    controls.taper in (:hanning, :blackman, :lanczos) ||
        throw(ArgumentError("taper must be :hanning, :blackman or :lanczos"))
    return FormulationOptions(controls)
end

validate(problem::TransientAnalysisProblem{<:LineParametersProblem}, ::Formula{:nlt}) = problem

function compute(selected::Formula{:nlt}, problem::TransientAnalysisProblem,
        execution::ComputationOptions)
    options = selected.options.data
    T = options.Tfactor * problem.window.tmax
    c = log(options.N^2) / T
    f = [(2k - 1) / (2T) for k in 1:options.N]
    line = problem.line
    scan = compute(LineParametersProblem(line.system, line.temperature, line.earth_props, f, c),
        options.line)
    workspace = TransientWorkspace(scan, problem, selected)
    t, voltages, currents = respond!(selected, workspace, selected.parameters, selected.options)
    alias = exp(-2c * T)
    retained = (nlt=(N=options.N, T=T, c=c, Tfactor=options.Tfactor,
        taper=options.taper, aliasing=alias / (1 - alias)),
        line=computation_details(LineParametersFormulation, scan))
    return t, voltages, currents, workspace.plan.coordinates, retained
end

unknown_count(assignment::TerminalImpedance) =
    iszero(assignment.R) && iszero(assignment.L) ? 1 : 0
unknown_count(::Waveform) = 1
unknown_count(::Matched) = 0
unknown_count(assignments::Tuple) = sum(unknown_count, assignments; init=0)

function validate(assignments::Tuple, ::Formula{:nlt}, coordinates)
    for assignment in assignments
        assignment.terminal in coordinates || throw(ArgumentError(
            "terminal $(assignment.terminal): the line formulation removed it"))
    end
    return assignments
end
validate(assignment::Matched, ::Formula{:nlt}, coordinates) = assignment

"""
$(TYPEDSIGNATURES)

Build the transient workspace from a computed line scan and terminal assignments.
"""
function TransientWorkspace(scan::LineParameters, problem::TransientAnalysisProblem,
        selected::Formula{:nlt})
    coordinates = details(scan).data.coordinates
    validate(problem.sending, selected, coordinates)
    validate(problem.receiving, selected, coordinates)
    options = selected.options.data
    T = options.Tfactor * problem.window.tmax
    input = (scan=scan, sending=problem.sending, receiving=problem.receiving,
        window=problem.window, T=T, c=log(options.N^2) / T, ℓ=line_length(problem.line))
    n = size(scan.Z, 1)
    plan = (n=n, N=options.N, nt=round(Int, problem.window.tmax / problem.window.step) + 1,
        coordinates=coordinates, nodes=2n,
        unknowns=2n + unknown_count(problem.sending) + unknown_count(problem.receiving))
    S = float(promote_type(eltype(scan.Z), eltype(scan.Y)))
    buffers = initialize_buffers(selected, S, input, plan, (;))
    return TransientWorkspace(input, plan, buffers)
end

function initialize_buffers(selected::Formula{:nlt}, ::Type{T}, input::NamedTuple,
        plan::NamedTuple, buffers) where {T}
    n, N, nt = plan.n, plan.N, plan.nt
    R = typeof(real(zero(T)))
    window = Val(selected.options.data.taper)
    return merge(buffers, (M=zeros(T, plan.unknowns, plan.unknowns),
        rhs=zeros(T, plan.unknowns), ZY=Matrix{T}(undef, n, n),
        Tv=Matrix{T}(undef, n, n), Tinv=Matrix{T}(undef, n, n),
        factor=Matrix{T}(undef, n, n), scaled=Matrix{T}(undef, n, n),
        Zslice=Matrix{T}(undef, n, n), Ys=Matrix{T}(undef, n, n),
        Ym=Matrix{T}(undef, n, n), Yc=Matrix{T}(undef, n, n),
        roots=Vector{T}(undef, n),
        taper=R[taper(selected, window, (2k - 1) / (2N)) for k in 1:N],
        voltage_sum=Matrix{T}(undef, 2, n), current_sum=Matrix{T}(undef, 2, n),
        voltage_spectrum=Array{T}(undef, N, 2, n),
        current_spectrum=Array{T}(undef, N, 2, n),
        t=R[(m - 1) * input.window.step for m in 1:nt],
        voltages=Array{R}(undef, 2, n, nt), currents=Array{R}(undef, 2, n, nt)))
end

taper(::Formula{:nlt}, ::Val{:hanning}, ratio) = (1 + cos(π * ratio)) / 2
taper(::Formula{:nlt}, ::Val{:blackman}, ratio) =
    0.42 + 0.5cos(π * ratio) + 0.08cos(2π * ratio)
taper(::Formula{:nlt}, ::Val{:lanczos}, ratio) = sinc(ratio)

laplace_transform(source::Step, s) = source.amplitude * exp(-s * source.t0) / s
function laplace_transform(source::Cosine, s)
    ω = 2π * source.frequency
    return source.amplitude * (s * cos(source.phase) - ω * sin(source.phase)) / (s^2 + ω^2)
end
laplace_transform(source::DoubleExponential, s) =
    source.amplitude * (1 / (s + source.alpha) - 1 / (s + source.beta))

function stamp!(selected::Formula{:nlt}, assignments::Tuple, workspace, s, offset, unknown)
    for assignment in assignments
        node = offset + something(findfirst(==(assignment.terminal), workspace.plan.coordinates))
        unknown = stamp!(selected, assignment, workspace, s, node, unknown)
    end
    return unknown
end

function stamp!(::Formula{:nlt}, ::Matched, workspace, s, offset, unknown)
    n = workspace.plan.n
    for j in 1:n, i in 1:n
        workspace.buffers.M[offset+i, offset+j] += workspace.buffers.Yc[i, j]
    end
    return unknown
end

function stamp!(::Formula{:nlt}, assignment::TerminalImpedance, workspace, s, node, unknown)
    if iszero(assignment.R) && iszero(assignment.L)
        unknown += 1
        workspace.buffers.M[unknown, node] = 1
        workspace.buffers.M[node, unknown] = 1
    else
        workspace.buffers.M[node, node] += inv(assignment.R + s * assignment.L)
    end
    return unknown
end

function stamp!(::Formula{:nlt}, source::Waveform, workspace, s, node, unknown)
    unknown += 1
    workspace.buffers.M[unknown, node] = 1
    workspace.buffers.M[node, unknown] = 1
    workspace.buffers.rhs[unknown] = laplace_transform(source, s)
    return unknown
end

function sample_network!(selected::Formula{:nlt}, workspace::TransientWorkspace, k, s)
    input, plan, arrays = workspace.input, workspace.plan, workspace.buffers
    n = plan.n
    Z = @view input.scan.Z.values[:, :, k]
    Y = @view input.scan.Y.values[:, :, k]
    mul!(arrays.ZY, Z, Y)
    F = eigen!(arrays.ZY)
    copyto!(arrays.Tv, F.vectors)
    copyto!(arrays.factor, arrays.Tv)
    fill!(arrays.Tinv, 0)
    for j in 1:n
        arrays.Tinv[j, j] = 1
        d = sqrt(complex(F.values[j]))
        arrays.roots[j] = real(d) < 0 ? -d : d
    end
    ldiv!(lu!(arrays.factor), arrays.Tinv)
    for j in 1:n
        d = arrays.roots[j]
        hm = exp(-d * input.ℓ)
        Am = d * (1 + hm^2) / (1 - hm^2)
        for i in 1:n
            arrays.scaled[i, j] = arrays.Tv[i, j] * Am
        end
    end
    mul!(arrays.Ys, arrays.scaled, arrays.Tinv)
    for j in 1:n
        d = arrays.roots[j]
        hm = exp(-d * input.ℓ)
        Bm = -2d * hm / (1 - hm^2)
        for i in 1:n
            arrays.scaled[i, j] = arrays.Tv[i, j] * Bm
        end
    end
    mul!(arrays.Ym, arrays.scaled, arrays.Tinv)
    for j in 1:n, i in 1:n
        arrays.scaled[i, j] = arrays.Tv[i, j] * arrays.roots[j]
    end
    mul!(arrays.Yc, arrays.scaled, arrays.Tinv)
    copyto!(arrays.Zslice, Z)
    Zfactor = lu!(arrays.Zslice)
    ldiv!(Zfactor, arrays.Ys)
    ldiv!(Zfactor, arrays.Ym)
    ldiv!(Zfactor, arrays.Yc)
    fill!(arrays.M, 0)
    fill!(arrays.rhs, 0)
    for j in 1:n, i in 1:n
        arrays.M[i, j] = arrays.Ys[i, j]
        arrays.M[n+i, n+j] = arrays.Ys[i, j]
        arrays.M[i, n+j] = arrays.Ym[i, j]
        arrays.M[n+i, j] = arrays.Ym[i, j]
    end
    unknown = stamp!(selected, input.sending, workspace, s, 0, plan.nodes)
    stamp!(selected, input.receiving, workspace, s, n, unknown)
    ldiv!(lu!(arrays.M), arrays.rhs)
    σ = arrays.taper[k]
    for i in 1:n
        arrays.voltage_spectrum[k, 1, i] = σ * arrays.rhs[i]
        arrays.voltage_spectrum[k, 2, i] = σ * arrays.rhs[n+i]
        I0 = zero(s)
        IL = zero(s)
        for j in 1:n
            I0 += arrays.Ys[i, j] * arrays.rhs[j] + arrays.Ym[i, j] * arrays.rhs[n+j]
            IL += arrays.Ym[i, j] * arrays.rhs[j] + arrays.Ys[i, j] * arrays.rhs[n+j]
        end
        arrays.current_spectrum[k, 1, i] = σ * I0
        arrays.current_spectrum[k, 2, i] = σ * IL
    end
    return workspace
end

function invert_spectra!(::Formula{:nlt}, workspace)
    input, plan, arrays = workspace.input, workspace.plan, workspace.buffers
    Δω = π / input.T
    for m in 1:plan.nt
        tt = arrays.t[m]
        factor = 2 / input.T * exp(input.c * tt)
        fill!(arrays.voltage_sum, 0)
        fill!(arrays.current_sum, 0)
        @inbounds for k in 1:plan.N
            e = cis((2k - 1) * Δω * tt)
            for i in 1:plan.n, side in 1:2
                arrays.voltage_sum[side, i] += arrays.voltage_spectrum[k, side, i] * e
                arrays.current_sum[side, i] += arrays.current_spectrum[k, side, i] * e
            end
        end
        for i in 1:plan.n, side in 1:2
            arrays.voltages[side, i, m] = factor * real(arrays.voltage_sum[side, i])
            arrays.currents[side, i, m] = factor * real(arrays.current_sum[side, i])
        end
    end
    return arrays.t, arrays.voltages, arrays.currents
end

function respond!(selected::Formula{:nlt}, workspace::TransientWorkspace, parameters, options)
    Δω = π / workspace.input.T
    for k in 1:workspace.plan.N
        s = complex(workspace.input.c, (2k - 1) * Δω)
        sample_network!(selected, workspace, k, s)
    end
    return invert_spectra!(selected, workspace)
end

:nlt
