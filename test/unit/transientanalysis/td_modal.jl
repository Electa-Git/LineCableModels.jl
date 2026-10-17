@testmodule TDModalFixtures begin
    using LineCableModels
    using LinearAlgebra: lu
    import LineCableModels.TransientAnalysis as TA

    function modal_line(; transformation=ones(1, 1), impedances=[
                (constant=100.0, poles=Float64[], residues=Float64[], delay=0.0)],
            propagation=[(constant=1.0, poles=Float64[], residues=Float64[], delay=1.035e-5)])
        n = size(transformation, 1)
        payload = ModalModel(transformation, 1e6, impedances, propagation)
        statistic = (rms=zeros(n, n), max=zeros(n, n), f_max=zeros(n, n))
        return LineModel(payload, ["wire:$i" for i in 1:n], 1000.0, [1.0, 1e6],
            (admittance=statistic, propagation=statistic),
            (checked=false, passive=nothing, bands=[], worst=nothing,
                f_worst=nothing, perturbation=nothing),
            (targets=nothing, achieved=nothing), ComputationDetails((;)))
    end

    function driven_line(line; receiving=(), source=Step("wire:1", 1.0),
            window=(tmax=1e-4, step=1e-7), initial_solve=true, resistance=0.0)
        problem = TransientAnalysisProblem(line; sending=(source,), receiving, window)
        selected = TransientFormulation(:td; options=(; initial_solve)).formula
        workspace = TA.TransientWorkspace(problem, selected)
        # The public waveform is ideal. Stamp its series source resistance in this fixture.
        workspace.buffers.M[end, end] = -resistance
        buffers = merge(workspace.buffers, (factor=lu(workspace.buffers.M),))
        return selected, TA.TransientWorkspace(workspace.input, workspace.plan, buffers)
    end
end

@testitem "TransientAnalysis / modal TD three bare wires (b) and allocations (c)" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.TransientAnalysis as TA
    import LineCableModels.Engine: compare

    f = 10.0 .^ range(0, 6; length=31)
    line = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(line, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    fitted = compute(LineModelProblem(segment), LineModelFormulation(:marti))
    window = (tmax=2e-4, step=1e-7)
    sending = (Step("cable:1:core", 1.0),
        TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    receiving = (TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    problem = TransientAnalysisProblem(fitted; sending, receiving, window)
    reference = compute(TransientAnalysisProblem(line; sending, receiving, window),
        TransientFormulation(:nlt; options=(N=1024,
            line=Formulation(earth_impedance=:unified, earth_admittance=:unified)));
        options=(timing=true,))
    response = compute(problem, TransientFormulation(:td); options=(timing=true,))

    @testset "Comparison (b)" begin
        @test response.t == reference.t
        @test response.coordinates == reference.coordinates
        @test all(isfinite, response.voltage)
        @test all(isfinite, response.current)
        diagnostics = details(fitted).data.marti
        errors = fitting_error(fitted)
        @info "E3 E2 fit evidence" admittance=errors.admittance propagation=errors.propagation offdiagonal_ratio=diagnostics.offdiagonal_ratio phase_error_degrees=diagnostics.phase_error_degrees
        for (quantity, actual, exact) in ((voltage, response.voltage, reference.voltage),
                (current, response.current, reference.current))
            scale = maximum(abs, exact)
            rms = compare(reference, response, quantity).absolute ./ scale
            peak_deviation = maximum(rms)
            pointwise = dropdims(maximum(abs, actual .- exact; dims=3); dims=3) ./ scale
            @info "E3 test (b), rows sending and receiving, columns terminals" quantity scale scaled_rms=repr(rms) peak_deviation scaled_pointwise=repr(pointwise) maximum_scaled_front=maximum(pointwise)
            # The D05 measure scales each terminal RMS by the largest reference peak.
            # Peak deviation is the maximum of that measure across ends and terminals.
            # Pointwise front maxima include NLT inversion and modal fit errors.
            if peak_deviation > 0.01
                @info "E3 deviations above one percent include E2 fit and transformation errors, trapezoidal discretization, and NLT inversion" quantity admittance_rms=errors.admittance.rms propagation_rms=errors.propagation.rms offdiagonal_ratio=diagnostics.offdiagonal_ratio phase_error_degrees=diagnostics.phase_error_degrees
            end
            for side in 1:2, terminal in 1:3
                @test rms[side, terminal] <= 0.10
            end
        end
        for side in 1:2, terminal in 2:3
            @test response.current[side, terminal, :] ≈ -response.voltage[side, terminal, :] atol=1e-12
        end
        @info "E3 computation wall seconds" nlt=details(reference).data.timing.wall_seconds td=details(response).data.timing.wall_seconds
    end

    @testset "One hundred warmed steps (c)" begin
        selected = TransientFormulation(:td).formula
        workspace = TA.TransientWorkspace(problem, selected)
        @test (@inferred TA.advance_modal!(selected, workspace, 2:101)) === nothing
        bytes = @allocated TA.advance_modal!(selected, workspace, 102:201)
        @info "E3 test (c)" steps=100 bytes allocations=(iszero(bytes) ? 0 : missing)
        @test bytes == 0
    end
end

@testitem "TransientAnalysis / modal TD lossless step (a)" tags=[:unit, :engine] setup=[TDModalFixtures] begin
    import LineCableModels.TransientAnalysis as TA
    line = TDModalFixtures.modal_line()
    delay = only(line.model.propagation).delay
    zc = only(line.model.impedance).constant
    for (receiving, expected) in (((), 1.0), (Matched(), 0.5))
        selected, workspace = TDModalFixtures.driven_line(line; receiving, resistance=zc)
        t, volts, amps = TA.respond!(selected, workspace, selected.parameters, selected.options)
        arrival = findfirst(>(1e-3), volts[2, 1, :])
        @test arrival !== nothing
        arrival_s = t[arrival]
        tail = t .>= 2delay
        error = maximum(abs, volts[2, 1, tail] .- expected)
        @info "E3 test (a)" receiving expected delay arrival_s error
        @test 0 <= arrival_s - delay <= workspace.input.window.step
        @test error <= 1e-3
        @test all(iszero, volts[2, 1, t .< delay])
        @test volts[1, 1, 1] ≈ 0.5 atol=1e-14
        @test amps[1, 1, 1] ≈ 0.5 / zc atol=1e-14
        if expected == 0.5
            reflected = maximum(abs, workspace.buffers.outgoing[2, 1, :])
            @info "E3 test (a) matched reflected wave" reflected
            @test reflected <= 1e-3
            @test maximum(abs, volts[1, 1, :] .- 0.5) <= 1e-3
        else
            @test maximum(abs, amps[2, 1, :]) <= 1e-14
        end
    end

    matched = compute(TransientAnalysisProblem(line; sending=Matched(), receiving=Matched(),
        window=(tmax=1e-4, step=1e-7)), TransientFormulation(:td))
    @test all(iszero, matched.voltage)
    @test all(iszero, matched.current)
    @test details(matched).data.td.terminations ==
        (sending=(impedance=:fitted_characteristic, poles=0),
            receiving=(impedance=:fitted_characteristic, poles=0))
end

@testitem "TransientAnalysis / modal TD sections and terminals" tags=[:unit, :engine] setup=[TDModalFixtures] begin
    import LineCableModels.TransientAnalysis as TA
    using LinearAlgebra: Diagonal

    @testset "One impedance section" begin
        h, p, r, k0 = 1e-6, -2000.0, 2e5, 100.0
        line = TDModalFixtures.modal_line(;
            impedances=[(constant=k0, poles=[p], residues=[r], delay=0.0)],
            propagation=[(constant=1.0, poles=Float64[], residues=Float64[], delay=1e-3)])
        selected, workspace = TDModalFixtures.driven_line(line;
            window=(tmax=1e-4, step=h))
        t, volts, amps = TA.respond!(selected, workspace, selected.parameters, selected.options)
        rate = p - r / k0
        ratio = (2 + rate * h) / (2 - rate * h)
        steady = inv(k0 - r / p)
        first_current = inv(k0 + r * h / (2 - p * h))
        expected = [steady + (first_current - steady) * ratio^(k-1) for k in eachindex(t)]
        @test amps[1, 1, :] ≈ expected rtol=1e-12
        @test all(==(1.0), volts[1, 1, :])
        @test all(iszero, volts[2, 1, :])

        dynamic = TDModalFixtures.modal_line(;
            impedances=line.model.impedance,
            propagation=[(constant=0.2, poles=[-1e4], residues=[8e3], delay=1.035e-5)])
        selected, workspace = TDModalFixtures.driven_line(dynamic; receiving=Matched())
        TA.respond!(selected, workspace, selected.parameters, selected.options)
        @test maximum(abs, workspace.buffers.outgoing[2, 1, :]) <= 1e-12
        @test all(isfinite, workspace.buffers.currents)
        # Constant outgoing wave 2, convolved independently with the sampled delayed ramp.
        tau = only(dynamic.model.propagation).delay
        lower = ceil(Int, tau / workspace.input.window.step) + 1
        a, cn, cp = TA.convolution_coefficients(selected, -1e4, workspace.input.window.step)
        state = 2cn
        @test workspace.buffers.voltages[2, 1, lower] ≈ 0.2 + 8e3state / 2 rtol=1e-12
        for k in (lower + 1):(lower + 4)
            state = a * state + 2(cn + cp)
            @test workspace.buffers.voltages[2, 1, k] ≈ 0.2 + 8e3state / 2 rtol=1e-12
        end
        @test TA.convolution_coefficients(selected, 0.0, h) == (1.0, h/2, h/2)
        @test TA.convolution_coefficients(selected, -1e-10, h)[2] ≈ h/2 rtol=1e-14
    end

    @testset "Fractional delayed history" begin
        line = TDModalFixtures.modal_line()
        source = Cosine("wire:1", 1.0, 2e4; phase=0.2)
        selected, workspace = TDModalFixtures.driven_line(line; receiving=Matched(), source)
        t, volts, amps = TA.respond!(selected, workspace, selected.parameters, selected.options)
        tau, h = only(line.model.propagation).delay, workspace.input.window.step
        for k in 120:140
            q = (t[k] - tau) / h
            lower = floor(Int, q)
            weight = q - lower
            expected = (1 - weight) * cos(2π * 2e4 * lower * h + 0.2) +
                weight * cos(2π * 2e4 * (lower + 1) * h + 0.2)
            @test volts[2, 1, k] ≈ expected atol=1e-12
        end
    end

    @testset "Nonorthogonal transformation" begin
        transformation = [1.0 0.25; 0.1 1.0]
        line = TDModalFixtures.modal_line(; transformation,
            impedances=[(constant=z, poles=Float64[], residues=Float64[], delay=0.0)
                for z in (100.0, 160.0)],
            propagation=[(constant=1.0, poles=Float64[], residues=Float64[], delay=tau)
                for tau in (1.03e-5, 1.27e-5)])
        problem = TransientAnalysisProblem(line;
            sending=(Step("wire:1", 1.0), Step("wire:2", 0.25)), receiving=Matched(),
            window=(tmax=8e-5, step=1e-7))
        response = compute(problem, TransientFormulation(:td; options=(initial_solve=true,)))
        g = transformation * Diagonal([1/100, 1/160]) * transpose(transformation)
        @test response.voltage[2, :, end] ≈ [1.0, 0.25] atol=1e-12
        @test response.current[1, :, end] ≈ g * [1.0, 0.25] atol=1e-12
        @test response.current[2, :, end] ≈ -g * [1.0, 0.25] atol=1e-12
        @test response.coordinates == line.coordinates
        @test details(response).data.line === details(line)
        reverse_problem = TransientAnalysisProblem(line; sending=Matched(),
            receiving=problem.sending, window=problem.window)
        reverse = compute(reverse_problem, TransientFormulation(:td; options=(initial_solve=true,)))
        @test reverse.voltage[1, :, :] ≈ response.voltage[2, :, :] atol=1e-12
        @test reverse.current[1, :, :] ≈ response.current[2, :, :] atol=1e-12
    end

    @testset "R, L, series R and L, and solid ground" begin
        h, tau, zc = 1e-7, 1e-5, 100.0
        line = TDModalFixtures.modal_line(;
            propagation=[(constant=1.0, poles=Float64[], residues=Float64[], delay=tau)])
        for (r, l) in ((40.0, 0.0), (0.0, 1e-3), (40.0, 1e-3), (0.0, 0.0))
            selected, workspace = TDModalFixtures.driven_line(line;
                source=Step("wire:1", 0.5), receiving=(TerminalImpedance("wire:1"; R=r, L=l),),
                window=(tmax=3e-5, step=h))
            t, volts, amps = TA.respond!(selected, workspace, selected.parameters, selected.options)
            arrival = round(Int, tau / h) + 1
            ratio = (2l / h - r - zc) / (2l / h + r + zc)
            steady = inv(r + zc)
            first_current = inv(r + zc + 2l / h)
            for k in arrival:(arrival + 20)
                load_current = steady + (first_current - steady) * ratio^(k-arrival)
                @test amps[2, 1, k] ≈ -load_current atol=1e-13
                @test volts[2, 1, k] ≈ 1 - zc * load_current atol=1e-12
            end
        end
    end

    @testset "Waveforms, initial solve and option validation" begin
        line = TDModalFixtures.modal_line()
        window = (tmax=1e-4, step=1e-7)
        t = [(k-1) * window.step for k in 1:1001]
        for (source, expected) in ((Step("wire:1", 2.0; t0=1e-6),
                    [tt >= 1e-6 ? 2.0 : 0.0 for tt in t]),
                (Cosine("wire:1", 2.0, 50.0; phase=0.2),
                    2 .* cos.(2π * 50 .* t .+ 0.2)),
                (DoubleExponential("wire:1", 2.0, 1e4, 1e6),
                    2 .* (exp.(-1e4 .* t) .- exp.(-1e6 .* t))))
            problem = TransientAnalysisProblem(line; sending=(source,), receiving=(), window)
            response = compute(problem, TransientFormulation(:td; options=(initial_solve=true,)))
            @test response.voltage[1, 1, :] ≈ expected atol=1e-12
        end
        problem = TransientAnalysisProblem(line; sending=(Step("wire:1", 1.0),), receiving=(), window)
        response = compute(problem, TransientFormulation(:td))
        @test all(iszero, response.voltage[:, :, 1])
        @test all(iszero, response.current[:, :, 1])
        @test response.voltage[1, 1, 2] == 1.0
        @test_throws ArgumentError("interp must be :single for a ModalModel") compute(problem,
            TransientFormulation(:td; options=(interp=:double,)))
        @test_throws ArgumentError("enforce_sign must be false for a ModalModel") compute(problem,
            TransientFormulation(:td; options=(enforce_sign=true,)))
        delay = only(line.model.propagation).delay
        invalid = TransientAnalysisProblem(line; sending=problem.sending, receiving=(),
            window=(tmax=1e-4, step=2delay))
        @test_throws ArgumentError("time step ($(2delay) s) larger than the smallest delay ($delay s)") compute(
            invalid, TransientFormulation(:td))
        boundary = TransientAnalysisProblem(line; sending=problem.sending, receiving=Matched(),
            window=(tmax=10delay, step=delay))
        @test all(isfinite, compute(boundary, TransientFormulation(:td)).voltage)
    end
end
