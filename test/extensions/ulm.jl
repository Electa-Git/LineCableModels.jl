@testitem "UniversalLineModel / fits, diagnostics and conversion" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel
    using Logging: ConsoleLogger, with_logger

    f = 10.0 .^ range(0, 6; length=31)
    problem = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(problem, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    data = LineFrequencyData(segment)
    @test data.f == f
    @test data.s == 2π * im .* f
    @test data.len == line_length(segment)
    @test data.Hm == permutedims(exp.(-gamma(segment) .* line_length(segment)))
    @test data.gamma == permutedims(gamma(segment))
    @test all(data.Yc[k] == Yc(segment, PhaseDomain)[:, :, k] for k in eachindex(f))
    @test all(data.H[k] == H(segment, PhaseDomain; field=:current)[:, :, k] for k in eachindex(f))
    @test all(data.T[k] == Ti(segment)[:, :, k] for k in eachindex(f))

    fitting = LineModelFormulation(:ulm; options=(yc=(npoles=8,), h=(npoles=8,)))
    model = with_logger(ConsoleLogger(stdout)) do
        compute(LineModelProblem(segment), fitting;
            options=(timing=true, verbosity=(default=1,)))
    end
    errors = fitting_error(model)
    @test maximum(errors.admittance.rms) < 1e-2
    @test maximum(errors.propagation.rms) < 1e-2
    @test all(0 .< delays(model) .< 1.5 * line_length(segment) / 299792458.0)
    @test passivity(model).checked
    @test passivity(model).perturbation !== nothing
    @test all(band -> length(band) == 2, passivity(model).bands)
    @test model.dc == (targets=nothing, achieved=nothing)
    @test !isempty(sprint(show, MIME"text/plain"(), report(ObservedResult(model, (fitting_error,)))))

    rebuilt = ULMModel(model)
    @test Set(keys(rebuilt.info)) == Set([:f])
    @test rebuilt.info[:f] == f
    payload = WidebandModel(rebuilt)
    restored = LineModel(payload, model.coordinates, model.line_length, model.frequencies,
        model.fitting_error, model.passivity, model.dc, model.details)
    roundtrip = ULMModel(restored)
    for s in data.s
        @test eval_yc(roundtrip, s) ≈ eval_yc(rebuilt, s) rtol=1e-12 atol=1e-12
        @test eval_h(roundtrip, s) ≈ eval_h(rebuilt, s) rtol=1e-12 atol=1e-12
    end
    for identifier in (:ulm, :fdcm)
        invalid = LineModelFormulation(identifier; options=(dc=(admittance=zeros(3, 3),),))
        @test_throws ArgumentError("dc targets are not implemented yet") compute(LineModelProblem(segment), invalid)
    end
    fdcm = compute(LineModelProblem(segment), LineModelFormulation(:fdcm;
        options=(yc=(npoles=8,), h=(npoles=8,), passivity=false)))
    @test fdcm isa LineModel{<:WidebandModel}
    @test passivity(fdcm).checked
    @test passivity(fdcm).perturbation === nothing
    @test all(isfinite, fitting_error(fdcm).propagation.rms)
end
