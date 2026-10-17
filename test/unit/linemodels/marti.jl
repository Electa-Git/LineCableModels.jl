@testitem "LineModels / Marti Bode corners and magnitude" tags=[:unit, :engine] begin
    import LineCableModels.LineModels as LM

    omega = 10.0 .^ range(0, 7; length=200)
    data = @. 100 * (1 + im * omega / 1e3) / (1 + im * omega / 1e5)
    formula = LM.Formula(:marti; options=(bode=(tolerance_db=3.0,),))
    tolerance = formula.options.data.bode.tolerance_db
    section = LM.bode_fit(formula, omega, data, tolerance)
    fitted = [section.constant + sum(r / (im * w - p)
        for (p, r) in zip(section.poles, section.residues); init=0.0) for w in omega]
    magnitude_error_db = maximum(abs.(20 .* log10.(abs.(fitted ./ data))))
    @test length(section.poles) == 1
    recovered_pole = -only(section.poles)
    recovered_zero = -(only(section.poles) - only(section.residues) / section.constant)
    @info "Marti test (a)" recovered_zero recovered_pole magnitude_error_db tolerance
    @test abs(recovered_pole / 1e5 - 1) <= 0.2
    @test abs(recovered_zero / 1e3 - 1) <= 0.2
    @test magnitude_error_db <= tolerance
    @test section.delay == 0
end

@testitem "LineModels / Marti lossless segment" tags=[:unit, :engine] begin
    import LineCableModels.Engine as E

    f = 10.0 .^ range(0, 6; length=31)
    omega = 2π .* f
    inductance, capacitance, length_m = 1e-6, 1e-11, 1000.0
    exact_delay = length_m * sqrt(inductance * capacitance)
    maps = ones(ComplexF64, 1, 1, length(f))
    domain = E.ModalDomain(ModalOperators(maps, maps),
        reshape(im .* omega .* sqrt(inductance * capacitance), 1, :))
    parameters = LineParameters(domain,
        SeriesImpedance(reshape(im .* omega .* inductance, 1, 1, :)),
        ShuntAdmittance(reshape(im .* omega .* capacitance, 1, 1, :)), f,
        ComputationDetails((phase_coordinates=["wire"],)))
    segment = PropagationParameters(parameters; line_length=length_m)
    result = compute(LineModelProblem(segment), LineModelFormulation(:marti))
    model = result.model
    impedance, propagation = only(model.impedance), only(model.propagation)
    relative_delay_error = abs(propagation.delay / exact_delay - 1)
    @info "Marti test (b)" impedance.constant propagation.constant exact_delay propagation.delay relative_delay_error
    @test model.transformation == ones(1, 1)
    @test model.transformation_frequency == last(f)
    @test impedance.constant ≈ sqrt(inductance / capacitance) rtol=1e-14
    @test isempty(impedance.poles) && isempty(impedance.residues)
    @test impedance.delay == 0
    @test propagation.constant ≈ 1.0 rtol=1e-14
    @test isempty(propagation.poles) && isempty(propagation.residues)
    @test relative_delay_error <= 1e-9
    @test maximum(fitting_error(result).admittance.max) < 1e-14
    @test maximum(fitting_error(result).propagation.max) < 1e-14
    @test !passivity(result).checked && passivity(result).passive === nothing
    @test result.dc == (targets=nothing, achieved=nothing)
    selected = compute(LineModelProblem(segment), LineModelFormulation(:marti;
        options=(transformation=(frequency=17.0,),)))
    @test selected.model.transformation_frequency == f[argmin(abs.(f .- 17.0))]
end

@testitem "LineModels / Marti three bare wires diagnostics" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.LineModels as LM

    f = 10.0 .^ range(0, 6; length=31)
    problem = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=1000.0)
    phase = compute(problem, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=1000.0)
    tolerances = (3.0, 1.0, 0.5, 0.1)
    results = map(tolerances) do tolerance
        compute(LineModelProblem(segment), LineModelFormulation(:marti;
            options=(bode=(tolerance_db=tolerance,),)))
    end
    # Preserve the original full-band fitting sequence bit for bit.
    input = (segment=segment, coordinates=details(segment).data.phase_coordinates,
        frequencies=f, length=1000.0)
    reference = LM.LineModelWorkspace(input, (;),
        LM.initialize_buffers(LM.Formula(:marti), Float64, input, (;), (;)))
    transformation = first(results).model.transformation
    LM.modal_samples!(LM.Formula(:marti), reference, transformation,
        inv(transpose(transformation)))
    omega = reference.buffers.omega
    @test all(a -> abs(a) >= 1e-2, reference.buffers.propagation)
    for (tolerance, result) in zip(tolerances, results)
        formula = LM.Formula(:marti; options=(bode=(tolerance_db=tolerance,),))
        for mode in 1:3
            impedance = @view reference.buffers.impedance[mode, :]
            propagation = @view reference.buffers.propagation[mode, :]
            delay = LM.refine_delay(formula, omega, propagation,
                reference.buffers.roots[mode, end], input.length)
            expected = (LM.bode_fit(formula, omega, impedance, tolerance),
                LM.bode_fit(formula, omega, propagation .* cis.(omega .* delay),
                    tolerance; delay))
            for (actual, original) in zip(
                    (result.model.impedance[mode], result.model.propagation[mode]), expected)
                @test reinterpret(UInt64, [actual.constant; actual.poles; actual.residues; actual.delay]) ==
                    reinterpret(UInt64, [original.constant; original.poles; original.residues; original.delay])
            end
        end
        @test details(result).data.marti.propagation_f_max == fill(last(f), 3)
    end
    admittance_rms = [maximum(fitting_error(result).admittance.rms) for result in results]
    propagation_rms = [maximum(fitting_error(result).propagation.rms) for result in results]
    for (index, tolerance) in enumerate(tolerances)
        result = results[index]
        errors = fitting_error(result)
        @info "Marti tolerance table" tolerance admittance_rms=admittance_rms[index] propagation_rms=propagation_rms[index] admittance_max=maximum(errors.admittance.max) propagation_max=maximum(errors.propagation.max) impedance_poles=[length(s.poles) for s in result.model.impedance] propagation_poles=[length(s.poles) for s in result.model.propagation]
    end
    @test all(<=(0), diff(admittance_rms))
    @test all(<=(0), diff(propagation_rms))
    @test admittance_rms[3] <= 0.01
    @test propagation_rms[3] <= 0.01
    eligible = findfirst(index -> admittance_rms[index] <= 0.01 &&
        propagation_rms[index] <= 0.01, 1:3)
    @test formulation_options(LineModelFormulation(:marti).formula).data.bode.tolerance_db ==
        tolerances[eligible]
    result = first(results)
    errors = fitting_error(result)
    transformation = result.model.transformation
    voltage = inv(transpose(transformation))
    expected_ratios = zeros(2)
    for (quantity, samples, left, right) in
            ((1, Z(phase), voltage, transformation), (2, Y(phase), transformation, voltage))
        for k in eachindex(f)
            coefficients = left \ (samples[:, :, k] * right)
            offdiagonal = maximum(abs(coefficients[i, j]) for i in 1:3, j in 1:3 if i != j)
            diagonal = maximum(abs(coefficients[i, i]) for i in 1:3)
            expected_ratios[quantity] = max(expected_ratios[quantity], offdiagonal / diagonal)
        end
    end
    @info "Marti test (c) reference ratios and fitting errors" expected_ratios admittance=errors.admittance propagation=errors.propagation
    diagnostic = details(result).data.marti
    @test diagnostic.offdiagonal_ratio ≈ expected_ratios rtol=1e-12
    @test diagnostic.offdiagonal_ratio ≈ [0.036, 4.6e-5] rtol=0.003
    @info "Marti test (c)" offdiagonal_ratio=diagnostic.offdiagonal_ratio phase_error_degrees=diagnostic.phase_error_degrees admittance=errors.admittance propagation=errors.propagation
    @test result.model isa ModalModel
    @test size(result.model.transformation) == (3, 3)
    @test result.model.transformation_frequency == last(f)
    @test all(isfinite, diagnostic.offdiagonal_ratio)
    @test all(isfinite, diagnostic.phase_error_degrees)
    @test all(>=(0), diagnostic.offdiagonal_ratio)
    @test all(>=(0), diagnostic.phase_error_degrees)
    for quantity in values(errors), statistic in values(quantity)
        @test all(isfinite, statistic)
        @test all(i == j || iszero(statistic[i, j]) for i in 1:3, j in 1:3)
    end
    @test all(section -> all(<(0), section.poles), result.model.impedance)
    @test all(section -> all(<(0), section.poles), result.model.propagation)
    @test passivity(result) == (checked=false, passive=nothing, bands=[],
        worst=nothing, f_worst=nothing, perturbation=nothing)
    @test result.dc == (targets=nothing, achieved=nothing)
end

@testitem "LineModels / Marti attenuated propagation" tags=[:unit, :engine] setup=[TestFixtures] begin
    import LineCableModels.LineModels as LM

    # A 10000 km three-wire segment underflows at the high end of this band.
    f = 10.0 .^ range(-3, 6; length=227)
    length_m = 1e7
    problem = TestFixtures.three_bare_wires_problem(; frequencies=f, line_length=length_m)
    phase = compute(problem, LineParametersFormulation())
    modal = compute(ModalAnalysisProblem(phase), ModalAnalysisFormulation(:default))
    segment = PropagationParameters(modal; line_length=length_m)
    for threshold in (1e-2, 1e-1)
        formulation = LineModelFormulation(:marti; options=(bode=(; threshold),))
        result = compute(LineModelProblem(segment), formulation)
        diagnostic = details(result).data.marti
        input = (segment=segment, coordinates=details(segment).data.phase_coordinates,
            frequencies=f, length=length_m)
        reference = LM.LineModelWorkspace(input, (;),
            LM.initialize_buffers(formulation.formula, Float64, input, (;), (;)))
        transformation = result.model.transformation
        LM.modal_samples!(formulation.formula, reference, transformation,
            inv(transpose(transformation)))
        @test any(iszero, reference.buffers.propagation)
        @test all(isfinite, diagnostic.phase_error_degrees)
        for quantity in values(fitting_error(result)), statistic in values(quantity)
            @test all(isfinite, statistic)
        end
        for mode in 1:3
            samples = @view reference.buffers.propagation[mode, :]
            first_below = findfirst(a -> abs(a) < threshold, samples)
            @test first_below !== nothing
            last_sample = first_below - 1
            @test all(a -> abs(a) >= threshold, samples[1:last_sample])
            @test abs(samples[last_sample + 1]) < threshold
            @test diagnostic.propagation_f_max[mode] == f[last_sample]
            omega = reference.buffers.omega[1:last_sample]
            propagation = samples[1:last_sample]
            delay = LM.refine_delay(formulation.formula, omega, propagation,
                reference.buffers.roots[mode, last_sample], length_m)
            expected = LM.bode_fit(formulation.formula, omega,
                propagation .* cis.(omega .* delay), 0.5; delay)
            @test isequal(result.model.propagation[mode], expected)
            impedance = LM.bode_fit(formulation.formula, reference.buffers.omega,
                reference.buffers.impedance[mode, :], 0.5)
            @test isequal(result.model.impedance[mode], impedance)
        end
        @info "Marti attenuated propagation bands" length_m threshold f_min=first(f) f_max=diagnostic.propagation_f_max underflows=count(iszero, reference.buffers.propagation)
    end
    for threshold in (0.0, 1.0, -0.1, 1.1, NaN, Inf, -Inf)
        formula = LM.Formula(:marti; options=(bode=(; threshold),))
        @test_throws ArgumentError LM.validate(formula.options, formula)
    end
    attenuated = PropagationParameters(modal; line_length=1e14)
    @test_throws ArgumentError compute(LineModelProblem(attenuated), LineModelFormulation(:marti))
end
