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
    reference = WidebandModel(fit_ulm(data; geom=nothing, yc=(npoles=8,),
        h=(npoles=8,), passivity=(;)))
    @test isequal(model.model.admittance, reference.admittance)
    @test isequal(model.model.propagation, reference.propagation)
    @test !details(model).data.step_check
    @test isnan(details(model).data.step_err)
    @test isnan(details(model).data.step_err_late)
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
    @test_throws ArgumentError("unknown :ulm dc options") LineModelFormulation(:ulm;
        options=(dc=(admittance=zeros(3, 3),),))
    invalid = LineModelFormulation(:fdcm; options=(dc=(method=:add_pole,),))
    @test_throws ArgumentError("dc targets are not implemented yet") compute(LineModelProblem(segment), invalid)
    fdcm = compute(LineModelProblem(segment), LineModelFormulation(:fdcm;
        options=(yc=(npoles=8,), h=(npoles=8,), passivity=false)))
    @test fdcm isa LineModel{<:WidebandModel}
    @test passivity(fdcm).checked
    @test passivity(fdcm).perturbation === nothing
    @test all(isfinite, fitting_error(fdcm).propagation.rms)
    @test !details(fdcm).data.step_check
    @test isnan(details(fdcm).data.step_err)
    @test isnan(details(fdcm).data.step_err_late)
    unused_line = LineParametersProblem(problem.system, problem.temperature + 10,
        problem.earth_props, f)
    unused_problem = LineModelProblem(segment; line=unused_line)
    unused_fdcm = @test_logs (:warn, ":fdcm ignores the supplied line problem") compute(
        unused_problem, LineModelFormulation(:fdcm;
            options=(yc=(npoles=8,), h=(npoles=8,), passivity=false)))
    @test isequal(unused_fdcm.model.admittance, fdcm.model.admittance)
    @test isequal(unused_fdcm.model.propagation, fdcm.model.propagation)
    @test isequal((details(unused_fdcm).data.step_check,
        details(unused_fdcm).data.step_err, details(unused_fdcm).data.step_err_late),
        (false, NaN, NaN))
    unused_marti = @test_logs (:warn, ":marti ignores the supplied line problem") compute(
        unused_problem, LineModelFormulation(:marti))
    @test typeof(unused_marti.model) <: ModalModel
end

@testitem "UniversalLineModel / DC correction from lossy cable parameters" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel

    # three_phase_system uses coaxial_design: insulation rho=1e8 Ohm m,
    # core radius 5 mm, insulation outer radius 10 mm, grounded sheath.
    f = 10.0 .^ range(0, 6; length=31)
    line = TestFixtures.line_parameters_problem(; frequencies=f)
    line_formulation = LineParametersFormulation(insulation_admittance=:lossy)
    phase = compute(line, line_formulation)
    @test real(phase.Y.values[1, 1, 1]) ≈ 2π / (1e8 * log(2)) rtol=1e-12
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=line_length(line))
    problem = LineModelProblem(segment; line)
    data = LineFrequencyData(segment)
    extension = Base.get_extension(LineCableModels, :LineCableModelsUniversalLineModelExt)
    # Nonorthogonal current modes have decay rates of 1 and 2 per metre.
    # This independent modal fixture distinguishes current and voltage propagation.
    coupled = extension.wave_matrices((ComplexF64[1 -1; -1 2], ComplexF64[5 4; 4 4]), 1.0)
    @test coupled.admittance ≈ [3 2; 2 2] rtol=1e-12
    @test coupled.propagation ≈ [exp(-1) exp(-2)-exp(-1); 0 exp(-2)] rtol=1e-12
    options = (yc=(npoles=8,), h=(npoles=8, step_check=false),
        line=line_formulation, passivity=false)
    baseline = compute(problem, LineModelFormulation(:ulm; options))
    native = WidebandModel(fit_ulm(data;
        geom=extension.LineCableModelsLine(line, line_formulation),
        yc=options.yc, h=options.h, passivity=false))
    explicit_empty = compute(problem, LineModelFormulation(:ulm;
        options=merge(options, (dc=(;),))))
    @test isequal(baseline.model.admittance, native.admittance)
    @test isequal(baseline.model.propagation, native.propagation)
    @test isequal(explicit_empty.model.admittance, baseline.model.admittance)
    @test isequal(explicit_empty.model.propagation, baseline.model.propagation)
    @test baseline.dc == explicit_empty.dc == (targets=nothing, achieved=nothing)
    @test !haskey(details(baseline).data, :dc)

    for dc in ((method=:add_pole,), (method=:functional_form,), (method=:add_pole, k=0.1))
        selected = LineModelFormulation(:ulm; options=merge(options, (; dc)))
        @test_throws ArgumentError("dc requires LineModelProblem(segment; line=...) and its line formulation option") compute(
            LineModelProblem(segment), selected)
        model = compute(problem, selected)
        rebuilt = ULMModel(model)
        evidence = details(model).data.dc
        convergence = evidence.convergence
        @test convergence.frequencies == [1e-6, 1e-9, 1e-12, 1e-15, 1e-18]
        @test convergence.pairs == ((1, 2), (2, 3), (3, 4), (4, 5))
        @test convergence.accepted_frequency == 1e-9
        for quantity in (:admittance, :propagation)
            target = model.dc.targets[quantity]
            achieved = model.dc.achieved[quantity]
            relative_error = maximum(abs, achieved - target) / maximum(abs, target)
            @test relative_error <= 1e-10
            @test all(isreal, target)
            @test convergence.imaginary_relative[quantity] < 1e-6
            @test evidence[quantity].method == dc.method
            if dc.method == :add_pole
                @test 0 < evidence[quantity].k < 1
                @test evidence[quantity].a0 ≈ 2π * first(f) * evidence[quantity].k rtol=1e-12
                if haskey(dc, :k)
                    @test evidence[quantity].k ≈ dc.k rtol=1e-12
                    @test isnothing(evidence[quantity].reference_range)
                else
                    @test evidence[quantity].reference_range == (first(f) / 1000, first(f))
                end
            end
            @info "ULM DC correction" method=dc.method quantity relative_error correction=evidence[quantity]
        end
        @test model.dc.achieved.admittance == eval_yc(rebuilt, 0)
        @test model.dc.achieved.propagation == eval_h(rebuilt, 0)
        @test maximum(convergence.relative_differences.impedance) <= 1e-6
        @test maximum(convergence.relative_differences.admittance) <= 1e-6
        @test all(isfinite, fitting_error(model).admittance.rms)
        @test all(isfinite, fitting_error(model).propagation.rms)
        @info "ULM DC limit" method=dc.method convergence
    end

    geom = extension.LineCableModelsLine(line, line_formulation)
    points = 2π * im .* [1e-6, 1e-9, 1e-12]
    samples = UniversalLineModel.zy(geom, points)
    @test length(geom.scans) == 3
    @test all(samples[k] === geom.scans[points[k]] for k in eachindex(points))
    @test UniversalLineModel.zy(geom, reverse(points)) == reverse(samples)
    reference = extension.dc_reference(geom, data)
    @test length(reference.admittance.s) == length(reference.admittance.F) == 31
    @test reference.propagation.s == reference.admittance.s
    @test length(reference.propagation.F) == 31
    @test length(geom.scans) == 34
    @test extrema(abs.(reference.admittance.s) ./ (2π)) == (first(f) / 1000, first(f))
    @test_throws ArgumentError UniversalLineModel.zy(geom, ComplexF64[])
    @test_throws ArgumentError UniversalLineModel.zy(geom, [1 + im, 2 + im])
    for dc in ((method=:unknown,), (k=0.1,), (method=:functional_form, k=0.1),
            (method=:add_pole, k=0.0), (method=:add_pole, k=1.0),
            (method=:add_pole, k=NaN), (method=:add_pole, reference=(;)))
        @test_throws ArgumentError LineModelFormulation(:ulm; options=(; dc))
    end
    @test_throws ArgumentError LineModelFormulation(:marti; options=(dc=(method=:add_pole,),))
end

@testitem "UniversalLineModel / DC limit with small insulation conductance" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel

    # Same coaxial geometry and grounded sheaths as the DC correction fixture.
    f = 10.0 .^ range(0, 6; length=31)
    dielectric = Material(kind=:insulator, rho=1e14, eps_r=3.0, mu_r=1.0)
    designs = [build(CableDesign, "phase-$phase", Stack(
        terminal(:core, Region(:metal, Disk(0.005), TestFixtures.conductor_material())),
        Region(:dielectric, Shell(0.005), dielectric),
        terminal(:sheath, Region(:return, Shell(0.001), TestFixtures.conductor_material())),
        Region(:cover, Shell(0.001), dielectric))) for phase in 1:3]
    system = build(LineCableSystem, designs, [(-0.08, -1.0), (0.0, -1.08), (0.08, -1.0)];
        system_id="small-conductance-three-phase", line_length=600.0,
        connections=[Dict(:core=>phase, :sheath=>0) for phase in 1:3])
    line = TestFixtures.line_parameters_problem(system; frequencies=f)
    line_formulation = LineParametersFormulation(insulation_admittance=:lossy)
    phase = compute(line, line_formulation)
    @test real(phase.Y.values[1, 1, 1]) ≈ 2π / (1e14 * log(2)) rtol=1e-12
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=line_length(line))
    problem = LineModelProblem(segment; line)
    for method in (:add_pole, :functional_form)
        model = compute(problem, LineModelFormulation(:ulm;
            options=(yc=(npoles=8,), h=(npoles=8, step_check=false),
                line=line_formulation, passivity=false, dc=(; method))))
        # At zero frequency every delay factor is one. Sum all pole terms by entry.
        roundoff_bound = 16 * eps(Float64) * maximum(sum(
            abs.(group.residues[:, :, k] ./ group.poles[k])
            for group in model.model.propagation for k in eachindex(group.poles)))
        @info "ULM small-conductance DC rounding bound" method roundoff_bound
        convergence = details(model).data.dc.convergence
        accepted = findfirst(==(convergence.accepted_frequency), convergence.frequencies)
        @test convergence.accepted_frequency < 1e-9
        for quantity in (:impedance, :admittance)
            differences = convergence.relative_differences[quantity]
            @test differences[accepted - 1] <= 1e-6
        end
        for quantity in (:admittance, :propagation)
            target = model.dc.targets[quantity]
            achieved = model.dc.achieved[quantity]
            relative_error = maximum(abs, achieved - target) / maximum(abs, target)
            if quantity == :propagation
                @test maximum(abs, achieved - target) <= roundoff_bound
            else
                @test relative_error <= 1e-10
            end
            @test all(isreal, target)
            @test convergence.imaginary_relative[quantity] < 1e-6
            @info "ULM small-conductance DC target" method quantity relative_error
        end
        @info "ULM small-conductance DC limit" method convergence
    end
end

@testitem "UniversalLineModel / time rows and exact step check" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel
    import LineCableModels: validate

    extension = Base.get_extension(LineCableModels, :LineCableModelsUniversalLineModelExt)
    f = 10.0 .^ range(0, 6; length=31)
    line = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    line_formulation = Formulation(earth_impedance=:unified, earth_admittance=:unified)
    phase = compute(line, line_formulation)
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    problem = LineModelProblem(segment; line)
    formulation = LineModelFormulation(:ulm;
        options=(yc=(npoles=8,), h=(td=true, td_N=1024),
            line=line_formulation, passivity=false))
    @test validate(problem, formulation.formula) === problem
    @info "ULM time-row fit scans"
    model = compute(problem, formulation; options=(verbosity=(default=1,),))
    diagnostics = details(model).data
    @test diagnostics.step_check
    @test isfinite(diagnostics.step_err)
    @test isfinite(diagnostics.step_err_late)
    @info "ULM time rows" step_err=diagnostics.step_err step_err_late=diagnostics.step_err_late

    geom = extension.LineCableModelsLine(line, line_formulation)
    fitted = ULMModel(model)
    t = collect(range(0.98minimum(fitted.tau),
        maximum(fitted.tau) + 0.3minimum(fitted.tau); length=5))
    exact = h_step_response_exact(geom, line_length(segment), t; N=1024)
    @info "ULM five-time fitted error" relative_error=maximum(abs, h_step_response(fitted, t) - exact)/maximum(abs, exact)
    @test UniversalLineModel.nphases(geom) == 3
    T = 1.25maximum(t)
    c = log(1024^2) / T
    points = [complex(c, (2k - 1) * (π / T)) for k in 1:1024]
    @test Set(keys(geom.scans)) == Set(points)
    @test h_step_response_exact(geom, line_length(segment), t; N=1024) == exact
    @test length(geom.scans) == 1024
    slow_point = complex(1.0, 2π * 100.0)
    Z, Y = UniversalLineModel.zy(geom, slow_point)
    direct = compute(LineParametersProblem(line.system, line.temperature, line.earth_props,
        [imag(slow_point) / (2π)], real(slow_point)), line_formulation)
    @test Z == direct.Z.values[:, :, 1]
    @test Y == direct.Y.values[:, :, 1]
    @test UniversalLineModel.zy(geom, slow_point) === geom.scans[slow_point]
    k = 16
    Zreal, Yreal = UniversalLineModel.zy(geom, 2π * im * f[k])
    phase_segment = transform(PhaseDomain, segment.parameters)
    @test size(Zreal) == size(Yreal) == size(phase_segment.Z.values)[1:2]
    @test Zreal ≈ phase_segment.Z.values[:, :, k] rtol=1e-12
    @test Yreal ≈ phase_segment.Y.values[:, :, k] rtol=1e-12

    @info "ULM default step-check scans"
    default_model = compute(problem, LineModelFormulation(:ulm;
        options=(yc=(npoles=8,), line=line_formulation, passivity=false));
        options=(verbosity=(default=1,),))
    @test details(default_model).data.step_check
    @test isfinite(details(default_model).data.step_err)
    @test isfinite(details(default_model).data.step_err_late)
    @info "ULM default step check" step_err=details(default_model).data.step_err step_err_late=details(default_model).data.step_err_late
    disabled = compute(problem, LineModelFormulation(:ulm;
        options=(yc=(npoles=8,), h=(step_check=false,), line=line_formulation,
            passivity=false)))
    @test !details(disabled).data.step_check
    @test isnan(details(disabled).data.step_err)
    @test isnan(details(disabled).data.step_err_late)

    # Both inversions use the same odd harmonics, abscissa and Hanning window.
    # At N=4096, N*eps(Float64)*exp(c*tmax) is below 6e-7.
    # A tolerance of 1e-5 relative to the peak allows for solves and summation order.
    inversion_tolerance = 1e-5
    late_end = 10 / (2π * 1e-3)
    windows = ((label=:arrival, start=first(t), stop=last(t), intervals=2048),
        (label=:late_1, start=20maximum(fitted.tau), stop=1e-3, intervals=64),
        (label=:late_2, start=1e-3, stop=0.1, intervals=64),
        (label=:late_3, start=0.1, stop=10.0, intervals=64),
        (label=:late_4, start=10.0, stop=late_end, intervals=64))
    references = map(windows) do window
        transient = TransientAnalysisProblem(line;
            sending=(Step("cable:1:core", 1.0), TerminalImpedance("cable:2:core"),
                TerminalImpedance("cable:3:core")), receiving=Matched(),
            window=(tmax=window.stop, step=window.stop / window.intervals))
        response = compute(transient, TransientFormulation(:nlt;
            options=(line=line_formulation, N=4096, Tfactor=1.25)))
        hook = h_step_response_exact(geom, line_length(segment), response.t; N=4096)
        received = permutedims(response.voltage[2, :, :])
        indices = findall(>=(window.start), response.t)
        peak = maximum(abs, received[indices, :])
        column_error = maximum(abs, hook[indices, :, 1] - received[indices, :]) / peak
        # NLT receiving voltages select row 1 of ULM's current-wave response.
        relative_error = maximum(abs, hook[indices, 1, :] - received[indices, :]) / peak
        @info "ULM independent NLT reference" window=window.label N=4096 T=1.25window.stop relative_error column_error inversion_tolerance
        @test isfinite(relative_error)
        @test isfinite(column_error)
        if window.label !== :arrival
            @test relative_error <= inversion_tolerance
        end
        return (t=response.t, exact=hook, received=received, relative_error=relative_error)
    end

    convergence = map(Iterators.product((1024, 2048, 4096), (1.25, 2.0))) do (N, Tfactor)
        values = h_step_response_exact(geom, line_length(segment), t; N, Tfactor)
        return (N=N, Tfactor=Tfactor, values=values)
    end
    convergence_peak = maximum(entry -> maximum(abs, entry.values), convergence)
    refined_values = last(convergence).values
    for entry in convergence
        relative_error = maximum(abs, entry.values - refined_values) / convergence_peak
        @info "ULM exact five-time convergence" N=entry.N Tfactor=entry.Tfactor T=entry.Tfactor*maximum(t) relative_error times=t column_1=entry.values[:, :, 1]
    end
    envelope_max = reduce((left, right) -> max.(left, right),
        map(entry -> entry.values, convergence))
    envelope_min = reduce((left, right) -> min.(left, right),
        map(entry -> entry.values, convergence))
    spread_per_time = [maximum(envelope_max[it, :, :] - envelope_min[it, :, :]) /
        convergence_peak for it in eachindex(t)]
    @info "ULM exact five-time spread" relative_spread=maximum(spread_per_time) spread_per_time

    arrival = first(references)
    default_fitted = ULMModel(default_model)
    # This line has no DC shunt conductance, so its exact final H is the identity.
    for (label, values, final_values) in (
            (:exact, arrival.exact[:, :, 1], [1.0, 0.0, 0.0]),
            (:td, h_step_response(fitted, arrival.t)[:, :, 1],
                real.(eval_h(fitted, 0im)[:, 1])),
            (:default, h_step_response(default_fitted, arrival.t)[:, :, 1],
                real.(eval_h(default_fitted, 0im)[:, 1])))
        for conductor in 1:3
            signal = abs.(values[:, conductor])
            onset = findfirst(>(1e-3), signal)
            half = iszero(final_values[conductor]) ? nothing :
                findfirst(>(0.5abs(final_values[conductor])), signal)
            onset_time = isnothing(onset) ? missing : arrival.t[onset]
            half_time = isnothing(half) ? missing : arrival.t[half]
            @info "ULM front timing" label conductor onset_time half_time final_value=final_values[conductor]
        end
    end
    @info "ULM front scan" td_delays=Tuple(fitted.tau) default_delays=Tuple(default_fitted.tau) time_step=arrival.t[2] td_highest_frequency=(2*1024-1)/(2*1.25maximum(t)) default_highest_frequency=(2*4096-1)/(2*1.25maximum(t)) fitted_band=(first(f),last(f)) sample_count=length(f)

    native = fit_ulm(hvdc_cable(), line_length(segment); f)
    @info "ULM native cable baseline" step_err=native.info[:step_err] step_err_late=native.info[:step_err_late]
    @test isfinite(native.info[:step_err])
    @test isfinite(native.info[:step_err_late])

    other_line = LineParametersProblem(line.system, line.temperature + 10,
        line.earth_props, f)
    @test_throws ArgumentError("line consistency mismatch at inputs.temperature") compute(
        LineModelProblem(segment; line=other_line), formulation)
    other_formulation = LineModelFormulation(:ulm;
        options=(line=Formulation(earth_impedance=:carson1926),))
    @test_throws ArgumentError("line consistency mismatch at formulations.formulas.earth_impedance.identifier") compute(
        problem, other_formulation)
    @test_throws ArgumentError("td=true requires LineModelProblem(segment; line=...)") compute(
        LineModelProblem(segment), formulation)
    unrecorded = PropagationParameters(segment.parameters, segment.line_length,
        ComputationDetails(Base.structdiff(details(segment).data,
            NamedTuple{(:inputs, :formulations)})))
    @test_logs (:warn, "line consistency cannot be checked: the segment has no input or formulation record") validate(
        LineModelProblem(unrecorded; line), formulation.formula)
end

@testitem "UniversalLineModel / fitting selection sweep" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel

    f = 10.0 .^ range(0, 6; length=31)
    problem = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(problem, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    segments = ParametricResult(Combinatorial(ModalAnalysisFormulation()), [segment])
    problems = Gridspace{LineModelProblem}(segments)
    formulations = LineModelFormulation(Grid((:ulm, :fdcm));
        options=(yc=(npoles=8,), h=(npoles=8,), passivity=false))
    started = time_ns()
    run = try
        compute(problems, formulations)
    finally
        @info "Line model selection sweep" requested_fits=2 wall_seconds=(time_ns() - started) / 1e9
    end
    @test run isa ParametricResult
    @test length(run) == 2
    @test all(value -> value isa LineModel, run)
    @test run.axes.problems === problems
    @test run.axes.formulations == collect(formulations)
    for j in 1:length(formulations)
        model = run[1, j]
        @test model === run[j]
        @test details(model).data.formula == formula_id(run.axes.formulations[j].formula)
        @test details(model).data.options == formulation_options(run.axes.formulations[j].formula).data
        @test line_length(model) == line_length(segment)
        for quantity in values(fitting_error(model)), statistic in values(quantity)
            @test all(isfinite, statistic)
        end
    end
end

@testitem "UniversalLineModel / nested fitting option sweep" tags=[:extension, :ulm] setup=[TestFixtures] begin
    using UniversalLineModel

    f = 10.0 .^ range(0, 6; length=31)
    line = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(line, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    segments = ParametricResult(Combinatorial(ModalAnalysisFormulation()), [segment])
    problems = Gridspace{LineModelProblem}(segments)
    formulations = LineModelFormulation(:ulm;
        options=(yc=(npoles=Grid((6, 10)),), h=(npoles=Grid((8, 12)),),
            passivity=Grid((false, true))))
    run = compute(problems, formulations)
    @test length(formulations) == length(run) == 8
    @test eltype(run) <: LineModel
    @test run.axes.problems === problems
    @test run.axes.formulations == collect(formulations)
    for (j, (yc, h, passivity)) in enumerate(Iterators.product((6, 10), (8, 12), (false, true)))
        model = run[1, j]
        @test model === run[j]
        @test formulation_options(run.axes.formulations[j].formula).data ==
            (yc=(npoles=yc,), h=(npoles=h,), passivity=passivity)
        @test details(model).data.options ==
            formulation_options(run.axes.formulations[j].formula).data
        @test model.passivity.checked
        @test (model.passivity.perturbation !== nothing) == passivity
        for quantity in values(fitting_error(model)), statistic in values(quantity)
            @test all(isfinite, statistic)
        end
    end
end
