@testitem "Monte Carlo / modal and transient compositions" tags=[:unit, :measurements] setup=[TestFixtures] begin
    using Measurements
    using Statistics: mean, std
    import LineCableModels as LCM
    base = TestFixtures.known_case_problem(Val(:two_bare_wires); frequencies=[50.0,500.0,5000.0])
    space = Gridspace{LineParametersProblem}(temperature -> LineParametersProblem(base.system;
        temperature, earth_props=base.earth_props, frequencies=base.frequencies),
        (Grid(20.0,0.1),))
    problem = ParametricProblem(space)
    line = Formulation(earth_impedance=:unified,earth_admittance=:unified)
    mc = MonteCarlo(line; trials=8,seed=42,distribution=:uniform,
        return_samples=true,return_histograms=true)
    keyword = compute(problem,mc;modal=:default)
    positional = compute(problem,mc,ModalAnalysisFormulation();modal_options=(rotate=true,))
    @test keyword isa MonteCarloResult
    @test isconcretetype(typeof(keyword))
    # Trial actions have concrete result types. Population return types depend on
    # runtime retention options; Symbol selections also resolve at runtime.
    nominal_modal = @inferred compute(base,keyword.formulation.inner)
    explicit_modal = @inferred compute(base,line,ModalAnalysisFormulation())
    @test gamma(nominal_modal) == gamma(explicit_modal)
    @test all(isequal(samples(keyword)[1][key],samples(positional)[1][key]) for key in keys(samples(keyword)[1]))
    @test gamma(keyword[1]) == gamma(positional[1])
    assignment = details(keyword).data.mode_assignment[1]
    @test size(assignment.overlap) == (2,3,8)
    @test assignment.cluster_sizes == ones(Int,2,3)
    @test assignment.frequencies == base.frequencies
    @test minimum(assignment.overlap) > 0.99
    @test all(assignment.permutation .== reshape([1,2],:,1,1))
    @test details(keyword).data.failure_summary[1].failed == 0
    @test trial_count(keyword,1) == 8
    @test isfinite(confidence(keyword,1).cdf_bound)
    @test confidence(keyword,1).conditioning === :successful_realizations
    order = [2,1]
    sign = reshape([-1.0,1.0],1,:,1)
    maps = ModalOperators(Tv(nominal_modal)[:,order,:].*sign,Ti(nominal_modal)[:,order,:].*sign)
    permuted = LineParameters(LCM.Engine.ModalDomain(maps,gamma(nominal_modal)[order,:]),
        nominal_modal.Z.values[order,order,:].*sign./permutedims(sign,(2,1,3)),
        nominal_modal.Y.values[order,order,:].*sign./permutedims(sign,(2,1,3)),
        frequencies(nominal_modal);details=details(nominal_modal))
    axis = (reference=nominal_modal,frequencies=frequencies(nominal_modal))
    buffers = LCM.Commons.initialize_buffers(mc,Float64,nominal_modal,(trials=1,axis=axis),(;))
    aligned = @inferred LCM.Engine.ModalAnalysis.align_modes(permuted,nominal_modal)
    @test aligned.permutation == repeat(order,1,3)
    @test minimum(aligned.overlap) ≈ 1.0
    @test Ti(aligned.parameters) ≈ Ti(nominal_modal)
    @test Tv(aligned.parameters) ≈ Tv(nominal_modal)
    @test aligned.parameters.Z ≈ nominal_modal.Z
    @test aligned.parameters.Y ≈ nominal_modal.Y
    @test gamma(aligned.parameters) == gamma(nominal_modal)
    LCM.UQ._record_sample!(buffers,permuted,1,axis)
    @test buffers.permutations[:,:,1] == repeat(order,1,3)
    @test buffers.products.Ti[:,:,:,1] ≈ Ti(nominal_modal)
    @test buffers.products.Tv[:,:,:,1] ≈ Tv(nominal_modal)
    @test buffers.products.gamma[:,:,1] == gamma(nominal_modal)
    for key in (:gamma,:Zc,:Yc,:Ti,:Tv)
        draws = samples(keyword)[1][key]
        summary = statistics(keyword)[1][key]
        for index in CartesianIndices(size(draws)[1:end-1])
            values = vec(draws[index.I...,:])
            @test summary.real[index].mean ≈ mean(real.(values))
            @test summary.imag[index].std ≈ std(imag.(values))
        end
    end
    for key in (:alpha,:beta,:velocity)
        @test all(isfinite, getproperty.(statistics(keyword)[1][key],:mean))
    end
    @test_throws ArgumentError compute(problem,mc;modal_options=(rotate=true,))
    @test_throws ArgumentError compute(problem,LinearError(line);modal=:default)
    @test_throws ArgumentError compute(problem,LinearError(line),ModalAnalysisFormulation())
    @test_throws ArgumentError("LinearError does not support modal composition; first-order mode tracking is missing") compute(problem,LinearError(ModalAnalysisFormulation()))

    sending = (Step("cable:1:core",1.0),TerminalImpedance("cable:2:core";R=1.0))
    receiving = (TerminalImpedance("cable:1:core";R=50.0),TerminalImpedance("cable:2:core";R=1.0))
    transients = TransientAnalysisProblem(space;sending,receiving,window=(tmax=2e-6,step=1e-8))
    transient = TransientFormulation(:nlt;options=(N=32,line))
    response = compute(ParametricProblem(transients),MonteCarlo(transient;trials=5,seed=13,
        distribution=:uniform,return_samples=true,return_histograms=true,retain_details=true))
    nominal_transient = @inferred compute(TransientAnalysisProblem(base;sending,receiving,
        window=(tmax=2e-6,step=1e-8)),transient)
    @test response[1] isa TransientResponse
    @test size(response[1].voltage) == size(nominal_transient.voltage)
    for key in (:voltage,:current)
        summaries = statistics(response)[1][key]
        @test size(summaries) == size(getproperty(response[1],key))
        @test all(summary -> isfinite(summary.mean) && isfinite(summary.std),summaries)
        for endpoint in 1:2, terminal in 1:2
            @test all(isfinite,getproperty.(summaries[endpoint,terminal,:],:mean))
        end
        for kpi in (:peak,:peak_time,:arrival_time,:final_value)
            product = statistics(response)[1].kpis[key][kpi]
            @test size(product) == (2,2)
            @test all(summary -> ismissing(summary) || isfinite(summary.mean),product)
        end
    end
    @test details(response).data.failure_summary[1].accepted == 5
    @test isfinite(confidence(response,1).cdf_bound)
    @test_throws ArgumentError("LinearError does not support transient responses; first-order transient propagation is missing") compute(ParametricProblem(transients),LinearError(transient))
end

@testitem "Monte Carlo / failed mode identity is counted without replacement" tags=[:unit, :measurements] setup=[TestFixtures] begin
    using Measurements
    using LinearAlgebra: Diagonal
    import LineCableModels as LCM
    struct RotatingLine <: AbstractFormulation end
    function LCM.compute(problem::LineParametersProblem, ::RotatingLine;
            options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
        angle = problem.temperature > 20 ? π/4 : 0.0
        basis = [cos(angle) -sin(angle); sin(angle) cos(angle)]
        diagonal = Diagonal(ComplexF64[0.1+0.2im,0.3+0.4im])
        matrix = basis*diagonal*transpose(basis)
        coefficients = repeat(reshape(matrix,2,2,1),1,1,length(problem.frequencies))
        return LineParameters(coefficients,copy(coefficients),problem.frequencies;
            details=ComputationDetails(;inputs=LCM.Engine.completed_inputs(problem)))
    end
    base = TestFixtures.known_case_problem(Val(:two_bare_wires);frequencies=[50.0,500.0,5000.0])
    space = Gridspace{LineParametersProblem}(temperature -> LineParametersProblem(base.system;
        temperature,earth_props=base.earth_props,frequencies=base.frequencies),
        (Grid(20.0,0.1),))
    mc = MonteCarlo(RotatingLine();trials=12,seed=42,distribution=:uniform,
        on_error=:resample,retain_details=true,return_samples=true)
    result = compute(ParametricProblem(space),mc;modal=:default)
    summary = details(result).data.failure_summary[1]
    @test summary.attempts == 12
    @test 0 < summary.failed < 12
    @test summary.accepted + summary.failed == 12
    @test trial_count(result,1) == summary.accepted
    @test size(samples(result)[1].gamma,3) == summary.accepted
    @test all(failure -> occursin("overlap",failure.error.message) &&
        occursin("without replacement",failure.error.message),details(result).data.failures[1])
    @test all(summary -> summary.n == trial_count(result,1),statistics(result)[1].gamma.real)
    rejected = MonteCarlo(RotatingLine();trials=5,seed=42,
        distribution=(_rng,mean,sigma) -> mean+sigma)
    @test_throws ArgumentError compute(ParametricProblem(space),rejected;modal=:default)
end

@testitem "Monte Carlo / transient responses from fitted modal lines" tags=[:unit, :measurements] setup=[TDModalFixtures] begin
    using Measurements
    line = TDModalFixtures.modal_line()
    sending = (Step("wire:1",Grid(1.0,0.01)),)
    problems = TransientAnalysisProblem(line;sending,receiving=(),window=(tmax=3e-5,step=1e-6))
    result = compute(ParametricProblem(problems),MonteCarlo(TransientFormulation(:td);
        trials=5,seed=42,distribution=:uniform))
    @test size(result[1].voltage)[1:2] == (2,1)
    @test all(summary -> isfinite(summary.mean) && isfinite(summary.std),statistics(result)[1].voltage)
    @test confidence(result,1).conditioning === :none
    zero_problem = TransientAnalysisProblem(line;sending=(Step("wire:1",0.0),),
        receiving=(),window=(tmax=3e-5,step=1e-6))
    zero_response = compute(ParametricProblem(zero_problem),MonteCarlo(TransientFormulation(:td);
        trials=5,seed=42,return_samples=true,return_histograms=true))
    @test all(ismissing,statistics(zero_response)[1].kpis.voltage.arrival_time)
    @test all(ismissing,histograms(zero_response)[1].kpis.current.arrival_time)
    @test all(isnan,samples(zero_response)[1].kpis.voltage.arrival_time)
end
