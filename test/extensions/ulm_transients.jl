@testitem "UniversalLineModel / fitted line transients" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel
    import LineCableModels.Engine: compare

    f = 10.0 .^ range(0, 6; length=31)
    line = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(line, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    fitting = LineModelFormulation(:ulm; options=(yc=(npoles=8,), h=(npoles=8,)))
    fitted = compute(LineModelProblem(segment), fitting)
    model = ULMModel(fitted)
    window = (tmax=2e-4, step=min(1e-7, minimum(model.tau) / 10))
    sending = (Step("cable:1:core", 1.0),
        TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    receiving = (TerminalImpedance("cable:2:core"; R=1.0),
        TerminalImpedance("cable:3:core"; R=1.0))
    reference_problem = TransientAnalysisProblem(line; sending, receiving, window)
    problem = TransientAnalysisProblem(fitted; sending, receiving, window)
    reference = compute(reference_problem, TransientFormulation(:nlt; options=(N=1024,
        line=Formulation(earth_impedance=:unified, earth_admittance=:unified)));
        options=(timing=true,))
    response = compute(problem, TransientFormulation(:td); options=(timing=true,))
    @test response.t == reference.t
    @test response.coordinates == reference.coordinates
    @test size(response.voltage) == size(reference.voltage)
    @test details(response).data.td ==
        (interp=:single, enforce_sign=false, initial_solve=false, steps=length(response.t))
    @test details(response).data.line === details(fitted)
    @info "D05 wall times" nlt=details(reference).data.timing.wall_seconds td=details(response).data.timing.wall_seconds
    for (selector, values) in ((voltage, reference.voltage), (current, reference.current))
        error = compare(reference, response, selector).absolute
        peaks = dropdims(maximum(abs, values; dims=3); dims=3)
        scale = maximum(peaks)
        @info "D05 comparison" quantity=selector absolute_rms=repr(error) reference_peak=repr(peaks) fraction=repr(error ./ peaks) scaled=repr(error ./ scale)
        for side in 1:2, terminal in 1:3
            @test error[side, terminal] <= 0.05scale
        end
    end
    voltage_error = compare(reference, response, voltage).absolute
    voltage_peaks = dropdims(maximum(abs, reference.voltage; dims=3); dims=3)
    for side in 1:2, terminal in 1:3
        @test voltage_error[side, terminal] <= 0.05voltage_peaks[side, terminal]
    end
    # Both line ends use current entering the line, as does the NLT two-port.
    for side in 1:2, terminal in 2:3
        @test response.current[side, terminal, :] ≈ -response.voltage[side, terminal, :] rtol=1e-9
    end
    message = "Matched() needs a fitted characteristic-admittance branch, which the UniversalLineModel network does not provide"
    for (sending, receiving) in ((Matched(), receiving), (sending, Matched()))
        matched = TransientAnalysisProblem(fitted; sending, receiving, window)
        @test_throws ArgumentError(message) compute(matched, TransientFormulation(:td))
    end
    large_step = 2minimum(model.tau)
    invalid = TransientAnalysisProblem(fitted; sending, receiving,
        window=(tmax=window.tmax, step=large_step))
    @test_throws ArgumentError("time step ($(large_step) s) larger than the smallest delay ($(minimum(model.tau)) s)") compute(
        invalid, TransientFormulation(:td))

    @testset "Source and series impedance mapping" begin
        for (source, native) in ((Step("cable:1:core", 2.0; t0=1e-6), StepSource(2.0, 1e-6)),
                (Cosine("cable:1:core", 2.0, 50.0; phase=0.2), CosineSource(2.0, 2π * 50.0, 0.2)),
                (DoubleExponential("cable:1:core", 2.0, 1e4, 1e6), DoubleExpSource(2.0, 1e4, 1e6)))
            source_end = (source, TerminalImpedance("cable:2:core"; L=1e-3),
                TerminalImpedance("cable:3:core"; R=1.0, L=1e-3))
            mapped = TransientAnalysisProblem(fitted; sending=source_end, receiving, window)
            result = compute(mapped, TransientFormulation(:td;
                options=(interp=:double, enforce_sign=true, initial_solve=true)))
            net = Network(7)
            add_line!(net, model, [1, 2, 3], [4, 5, 6]; interp=:double, enforce_sign=true)
            add_vsource!(net, 1, native)
            add_inductor!(net, 2, 0, 1e-3)
            add_resistor!(net, 3, 7, 1.0)
            add_inductor!(net, 7, 0, 1e-3)
            add_resistor!(net, 5, 0, 1.0)
            add_resistor!(net, 6, 0, 1.0)
            direct = simulate(net, window.step, window.tmax; initial_solve=true)
            for side in 1:2, terminal in 1:3
                @test result.voltage[side, terminal, :] == direct.v[:, (side - 1) * 3 + terminal]
            end
            @test result.current[1, :, :] == permutedims(direct.lines[1].i0)
            @test result.current[2, :, :] == permutedims(direct.lines[1].iL)
        end
    end

    @testset "Solid ground in both formulations" begin
        ground = TerminalImpedance("cable:2:core")
        for (side, sending_end, receiving_end) in (
                (1, (sending[1], ground, sending[3]), receiving),
                (2, sending, (ground, receiving[2])))
            exact = compute(TransientAnalysisProblem(line; sending=sending_end,
                    receiving=receiving_end, window),
                TransientFormulation(:nlt; options=(N=1024,
                    line=Formulation(earth_impedance=:unified, earth_admittance=:unified))))
            @test all(iszero, exact.voltage[side, 2, :])
            result = compute(TransientAnalysisProblem(fitted; sending=sending_end,
                receiving=receiving_end, window), TransientFormulation(:td))
            @test all(iszero, result.voltage[side, 2, :])
        end
    end
end
