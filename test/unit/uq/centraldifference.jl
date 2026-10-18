@testmodule CentralDifferenceFixtures begin
    using LineCableModels
    function line(scale, spacing)
        scale > 0 || throw(DomainError(scale, "radial dimensions must be positive"))
        metal = Material(kind=:conductor,rho=2e-8,eps_r=1.0,mu_r=1.0,T0=20.0,alpha=.004)
        dielectric = Material(kind=:insulator,rho=1e8,eps_r=3.0)
        designs = [build(CableDesign,"wire-$i",Stack(
            terminal(:core,Region(:metal,Disk(.005scale),metal)),
            Region(:insulation,Shell(.005),dielectric))) for i in 1:2]
        system = build(LineCableSystem,designs,[(0.0,-1.0),(spacing,-1.3)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        return LineParametersProblem(system;frequencies=[50.0,500.0,5000.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0),temperature=20.0)
    end
    formulation() = Formulation(earth_impedance=:unified,earth_admittance=:unified)
end

@testitem "Linear error / central differences / line Jacobian and refinement" tags=[:unit,:measurements] setup=[CentralDifferenceFixtures] begin
    using Measurements
    calls = Tuple{Float64,Float64}[]
    builder = (scale,spacing) -> begin
        scale isa Measurement || push!(calls,(scale,spacing))
        CentralDifferenceFixtures.line(scale,spacing)
    end
    space = Gridspace{LineParametersProblem}(builder,
        (Grid(1.0,(1.0,)),Grid(.2,AbsoluteError(.002))))
    problem = ParametricProblem(space)
    inner = CentralDifferenceFixtures.formulation()
    direct = compute(problem,LinearError(inner))
    @test isempty(calls)
    central = compute(problem,LinearError(inner;derivative=:central))
    @test calls == [(1.0,.2),(1.0001,.2),(.9999,.2),(1.0,.20002),(1.0,.19998000000000002)]
    empty!(calls)
    refined = compute(problem,LinearError(inner;options=(derivative=:central,relative_step=.005)))
    @test length(calls) == 5
    # h/x = 1e-4 for both inputs. Smooth central truncation is O((h/x)^2),
    # about 1e-8. A 1e-5 relative tolerance also covers solver and subtraction error.
    tolerance = 1e-5
    agreement_errors = Float64[0.0]
    refinement_errors = Float64[0.0]
    for quantity in (Z,Y), component in (real,imag)
        exact = uncertainty.(component.(observe(direct[1],quantity)))
        approximate = uncertainty.(component.(observe(central[1],quantity)))
        finer = uncertainty.(component.(observe(refined[1],quantity)))
        @test all(isfinite, approximate)
        @test all(isapprox.(approximate,exact;rtol=tolerance,atol=0))
        @test all(isapprox.(finer,approximate;rtol=tolerance,atol=0))
        for index in eachindex(exact)
            iszero(exact[index]) || push!(agreement_errors,abs(approximate[index]/exact[index]-1))
            iszero(approximate[index]) || push!(refinement_errors,abs(finer[index]/approximate[index]-1))
        end
    end
    @info "Central difference uncertainty comparison" worst_agreement=maximum(agreement_errors) worst_refinement=maximum(refinement_errors) tolerance
    @test Measurements.cov(R(central[1])[1,1,1],R(central[1])[2,2,1]) ≈
        Measurements.cov(R(direct[1])[1,1,1],R(direct[1])[2,2,1]) rtol=tolerance
    @test isconcretetype(eltype(central.values))
    record = details(central).data.points[1].data
    @test record.evaluations == 5
    @test length(record.perturbations) == 4
    @test getproperty.(record.inputs,:argument_path) == [(1,),(2,)]
    @test nominal.(central[1].Z) ≈ compute(CentralDifferenceFixtures.line(1.0,.2),inner).Z
    for invalid in (:unknown, nothing, "central")
        @test_throws ArgumentError LinearError(inner;derivative=invalid)
    end
    for invalid in (0,-1,Inf,NaN,true,"0.01")
        @test_throws ArgumentError LinearError(inner;relative_step=invalid)
    end
    @test_throws ArgumentError LinearError(inner;derivative=:central,options=(derivative=:central,))
    @test LinearError(inner).derivative isa LineCableModels.UQ.MeasurementsPropagation
    @test LinearError(inner;derivative=:central).derivative.relative_step === 0.01
    @test_throws ArgumentError compute(problem,LinearError(inner;derivative=:central);modal_options=(rotate=true,))
end

@testitem "Linear error / central differences / modal and NLT compositions" tags=[:unit,:measurements] setup=[CentralDifferenceFixtures] begin
    using Measurements
    import LineCableModels as LCM
    space = Gridspace{LineParametersProblem}(scale -> CentralDifferenceFixtures.line(scale,.2),
        (Grid(1.0,AbsoluteError(.01)),))
    problem = ParametricProblem(space)
    line = CentralDifferenceFixtures.formulation()
    selected = LinearError(line;derivative=:central,retain_details=true)
    keyword = compute(problem,selected;modal=:default,modal_options=(rotate=true,))
    positional = compute(problem,selected,ModalAnalysisFormulation();modal_options=(rotate=true,))
    @test keyword isa LinearErrorResult
    @test nominal.(gamma(keyword[1])) == nominal.(gamma(positional[1]))
    @test uncertainty.(gamma(keyword[1])) == uncertainty.(gamma(positional[1]))
    record = details(keyword).data.points[1].data
    @test record.evaluations == 3
    @test length(record.computations) == 3
    @test size(record.mode_assignment.overlap) == (2,3)
    @test minimum(record.mode_assignment.overlap) > .99
    @test all(record.mode_assignment.permutation .== reshape([1,2],:,1))
    for quantity in (Z,Y,gamma,Zc,Yc,Ti,Tv)
        values = quantity(keyword[1])
        @test all(isfinite,nominal.(values))
        @test all(isfinite,uncertainty.(values))
    end
    for quantity in (alpha,beta,velocity)
        @test all(isfinite, uncertainty.(quantity(PropagationParameters(keyword[1]))))
    end
    # Direct modal problems use the same nominal assignment as composed problems.
    modal_space = Gridspace{ModalAnalysisProblem}(scale ->
        ModalAnalysisProblem(compute(CentralDifferenceFixtures.line(scale,.2),line)),
        (Grid(1.0,AbsoluteError(.01)),))
    direct = compute(ParametricProblem(modal_space),LinearError(ModalAnalysisFormulation();derivative=:central))
    @test nominal.(gamma(direct[1])) == nominal.(gamma(keyword[1]))
    @test uncertainty.(gamma(direct[1])) ≈ uncertainty.(gamma(keyword[1]))
    @test_throws ArgumentError("LinearError does not support modal composition; first-order mode tracking is missing") compute(problem,LinearError(line);modal=:default)

    sending = (Step("cable:1:core",1.0),TerminalImpedance("cable:2:core";R=1.0))
    receiving = (TerminalImpedance("cable:1:core";R=50.0),TerminalImpedance("cable:2:core";R=1.0))
    transient_space = TransientAnalysisProblem(space;sending,receiving,window=(tmax=2e-6,step=1e-8))
    transient = TransientFormulation(:nlt;options=(N=32,line))
    response = compute(ParametricProblem(transient_space),LinearError(transient;derivative=:central))
    @test response isa LinearErrorResult
    @test response[1] isa TransientResponse
    for quantity in (:voltage,:current)
        values = getproperty(response[1],quantity)
        @test all(isfinite,nominal.(values))
        @test all(isfinite,uncertainty.(values))
        @test any(>(0),uncertainty.(values))
    end
    @test details(response).data.points[1].data.evaluations == 3
    @test_throws ArgumentError("LinearError does not support transient responses; first-order transient propagation is missing") compute(ParametricProblem(transient_space),LinearError(transient))
end

@testitem "Linear error / central differences / input diagnostics and independence" tags=[:unit,:measurements] setup=[CentralDifferenceFixtures,TestFixtures] begin
    using Measurements
    import LineCableModels as LCM
    inner = CentralDifferenceFixtures.formulation()
    selected = LinearError(inner;derivative=:central)
    bounded = Gridspace{LineParametersProblem}(scale -> begin
        scale <= 1.0 || throw(DomainError(scale,"scale must not exceed one"))
        CentralDifferenceFixtures.line(scale,.2)
    end,(Grid(1.0,AbsoluteError(.01)),))
    failure = try
        compute(ParametricProblem(bounded),selected)
    catch exception
        exception
    end
    @test failure isa LCM.UQ.PerturbationError
    diagnostic = details(failure).data
    @test diagnostic.input == "input[1]"
    @test diagnostic.argument_path == (1,)
    @test diagnostic.step == .0001
    @test diagnostic.stage === :build
    @test failure.cause isa DomainError
    @test occursin("input[1]",sprint(showerror,failure))
    @test occursin("scale must not exceed one",sprint(showerror,failure))

    shared = Measurements.correlated_values([1.0,.2],[.0001 .00001; .00001 .000004])
    correlated = Gridspace{LineParametersProblem}(CentralDifferenceFixtures.line,
        (Grid((shared[1],)),Grid((shared[2],))))
    @test_throws ArgumentError compute(ParametricProblem(correlated),selected)
    @test_throws r"shared or correlated" compute(ParametricProblem(correlated),selected)
    # A supplied uncertain core also contains an undeclared source graph.
    measured = compute(ParametricProblem(Gridspace{LineParametersProblem}(
        scale -> CentralDifferenceFixtures.line(scale,.2),(Grid(1.0,AbsoluteError(.01)),))),LinearError(inner))[1]
    transported = Gridspace{ModalAnalysisProblem}(ModalAnalysisProblem,(Grid((measured,)),))
    @test_throws r"shared or correlated" compute(ParametricProblem(transported),
        LinearError(ModalAnalysisFormulation();derivative=:central))
    calls = Ref(0)
    exact = Gridspace{LineParametersProblem}(scale -> begin
        calls[] += 1
        CentralDifferenceFixtures.line(scale,.2)
    end,(Grid(1.0,AbsoluteError(0.0)),))
    result = compute(ParametricProblem(exact),selected)
    @test calls[] == 1
    @test details(result).data.points[1].data.evaluations == 1
    @test all(iszero,uncertainty.(result[1].Z))

    cable_space = Gridspace{CableConstantsProblem}(temperature ->
        CableConstantsProblem(TestFixtures.coaxial_design();temperature),
        (Grid((20.0,40.0),AbsoluteError(1.0)),))
    cable_problem = ParametricProblem(cable_space)
    direct_cables = compute(cable_problem,LinearError(CableConstantsFormulation()))
    central_cables = compute(cable_problem,LinearError(CableConstantsFormulation();derivative=:central))
    @test length(central_cables) == 2
    @test all(record -> record.data.evaluations == 3,details(central_cables).data.points)
    for index in eachindex(central_cables), quantity in (R,L,C,G)
        @test uncertainty.(quantity(central_cables[index])) ≈
            uncertainty.(quantity(direct_cables[index])) rtol=1e-5 atol=0
    end
end
