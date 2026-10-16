@testitem "TransientAnalysis / terminal assignments" tags=[:unit, :engine] begin
    import LineCableModels: validate

    grounded = TerminalImpedance("cable:1:core")
    @test (grounded.terminal, grounded.R, grounded.L) == ("cable:1:core", 0, 0)
    @test (grounded.R, grounded.L) === (0.0, 0.0)
    @test TerminalImpedance("cable:1:core"; R=10).R === 10.0
    impedance = TerminalImpedance("cable:2:core"; R=2, L=0.5)
    @test (impedance.R, impedance.L) == (2.0, 0.5)
    @test validate(impedance) === impedance
    for value in (-1, Inf, NaN)
        @test_throws DomainError TerminalImpedance("cable:1:core"; R=value)
        @test_throws DomainError TerminalImpedance("cable:1:core"; L=value)
    end
    resistances = TerminalImpedance("cable:1:core"; R=Grid((1, 2)))
    @test resistances isa Gridspace{TerminalImpedance}
    @test [assignment.R for assignment in resistances] == [1, 2]
    inductances = TerminalImpedance("cable:1:core"; L=Grid((0.1, 0.2)))
    @test [assignment.L for assignment in inductances] == [0.1, 0.2]

    step = Step("cable:1:core", 10; t0=0.1)
    cosine = Cosine("cable:1:core", 20, 50; phase=0.2)
    pulse = DoubleExponential("cable:1:core", 30, 1, 2)
    @test all(source -> source isa Waveform, (step, cosine, pulse))
    @test (step.terminal, step.amplitude, step.t0) == ("cable:1:core", 10, 0.1)
    @test (cosine.amplitude, cosine.frequency, cosine.phase) == (20, 50, 0.2)
    @test (pulse.amplitude, pulse.alpha, pulse.beta) == (30, 1, 2)
    @test Step("cable:1:core", 1).t0 == 0
    @test Cosine("cable:1:core", 1, 50).phase == 0
    @test (Step("cable:1:core", 1).amplitude, Step("cable:1:core", 1).t0) === (1.0, 0.0)
    @test (Cosine("cable:1:core", 1, 50).amplitude,
        Cosine("cable:1:core", 1, 50).frequency,
        Cosine("cable:1:core", 1, 50).phase) === (1.0, 50.0, 0.0)
    @test (pulse.amplitude, pulse.alpha, pulse.beta) === (30.0, 1.0, 2.0)
    @test validate(step) === step
    @test validate(cosine) === cosine
    @test validate(pulse) === pulse
    for value in (Inf, NaN)
        @test_throws DomainError Step("cable:1:core", value)
        @test_throws DomainError Step("cable:1:core", 1; t0=value)
        @test_throws DomainError Cosine("cable:1:core", value, 50)
        @test_throws DomainError Cosine("cable:1:core", 1, 50; phase=value)
        @test_throws DomainError DoubleExponential("cable:1:core", value, 1, 2)
    end
    for value in (-1, 0, Inf, NaN)
        @test_throws DomainError Cosine("cable:1:core", 1, value)
        @test_throws DomainError DoubleExponential("cable:1:core", 1, value, 2)
        @test_throws DomainError DoubleExponential("cable:1:core", 1, 1, value)
    end
    @test_throws DomainError DoubleExponential("cable:1:core", 1, 2, 2)
    @test_throws DomainError DoubleExponential("cable:1:core", 1, 3, 2)
    for assignments in (Step("cable:1:core", Grid((1, 2))),
            Cosine("cable:1:core", 1, Grid((50, 60))),
            DoubleExponential("cable:1:core", 1, 1, Grid((2, 3))))
        @test assignments isa Gridspace
        @test length(collect(assignments)) == 2
    end
    @test Matched() === Matched()
end
