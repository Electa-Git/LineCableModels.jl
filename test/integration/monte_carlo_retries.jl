@testitem "UQ / retries preserve accepted trials, source identities and seeded replay" tags=[:integration] setup=[TestFixtures] begin
    using Measurements
    using Statistics
    design = TestFixtures.coaxial_design()
    attempted_temperatures = Float64[]
    space = Gridspace{CableConstantsProblem}(
        temperature -> begin
            push!(attempted_temperatures, temperature)
            length(attempted_temperatures) <= 2 && throw(DomainError(
                temperature, "injected rejection before resolving the cable"))
            CableConstantsProblem(design; temperature)
        end,
        (Grid((20.0, 40.0), AbsoluteError((0.5,))),),
    )
    inner = CableConstantsFormulation()
    formulation = MonteCarlo(inner; trials=4, seed=71, distribution=:uniform,
        return_samples=true, return_histograms=true, bins=2,
        options=(on_error=:retry, retain_details=true, max_failures=3))
    sampled = compute(ParametricProblem(space), formulation)
    @test length(sampled) == 2
    @test sampled.trial_counts == [4, 4]
    @test length(attempted_temperatures) == 12
    @test attempted_temperatures[[7,12]] == [20.,40.]
    @test length(unique(sampled.point_seeds)) == 2
    @test all(isconcretetype, (eltype(sampled.values), eltype(sampled.stats),
        eltype(sampled.sample_values), eltype(sampled.histogram_values)))
    @test length.(sampled.details.data.trials) == [4, 4]
    @test length.(sampled.details.data.failures) == [2, 0]
    first_summary, second_summary = sampled.details.data.failure_summary
    @test first_summary.attempts == 6
    @test first_summary.accepted == 4
    @test first_summary.failed == 2
    @test first_summary.acceptance_rate == 4 / 6
    @test first_summary.by_type == [(type="DomainError", count=2)]
    @test first_summary.by_stage == [(stage=:build, count=2)]
    @test second_summary.attempts == second_summary.accepted == 4
    @test second_summary.failed == 0
    for (index, failure) in enumerate(first(sampled.details.data.failures))
        @test failure.attempt == index
        @test failure.target_trial == 1
        @test failure.stage === :build
        @test failure.sample !== nothing
        @test occursin("injected rejection", failure.error.message)
        @test !isempty(failure.error.stack)
        @test length(failure.error.stack) <= 8
        @test all(frame -> frame.line > 0 && !isempty(frame.function_name), failure.error.stack)
    end
    for point in 1:2
        accepted_temperatures = attempted_temperatures[point==1 ? (3:6) : (8:11)]
        expected = [compute(CableConstantsProblem(design; temperature), inner)
                    for temperature in accepted_temperatures]
        for quantity in (R, L, C, G)
            retained = observe(sampled, samples, quantity, point)
            @test retained == hcat((observe(value, quantity) for value in expected)...)
            @test nominal.(observe(sampled.values[point], quantity)) ≈ vec(mean(retained; dims=2))
            @test uncertainty.(observe(sampled.values[point], quantity)) ≈ vec(std(retained; dims=2))
        end
    end
    before_observation=length(attempted_temperatures)
    observed=observables(sampled)
    @test length(attempted_temperatures)==before_observation
    @test observed[1].gridpoint.sampling.diagnostics.failure_summary==first_summary
    @test observed[1].gridpoint.sampling.diagnostics.failures[1].sample==sampled.details.data.failures[1][1].sample
    recorded = copy(attempted_temperatures)
    empty!(attempted_temperatures)
    replay = compute(ParametricProblem(space), formulation)
    @test attempted_temperatures == recorded
    @test samples(replay) == samples(sampled)
    @test replay.point_seeds == sampled.point_seeds
    @test replay.details.data.failure_summary == sampled.details.data.failure_summary

end

@testitem "UQ / retries do not hide non-domain errors or run past their limit" tags=[:integration] begin
    using Measurements
    inner = CableConstantsFormulation()
    for (injected, policy, expected_attempts) in (
        (DomainError(-1, "invalid input marker"), :fail, 1),
        (ArgumentError("programming error marker"), :retry, 1),
        (DomainError(-1, "unrecoverable input marker"), :retry, 3),
    )
        attempts = Ref(0)
        space = Gridspace{CableConstantsProblem}(
            _ -> begin
                attempts[] += 1
                throw(injected)
            end,
            (Grid((1.0,)),),
        )
        formulation = MonteCarlo(inner; trials=2, seed=71,
            options=(on_error=policy, retain_details=true, max_failures=3))
        failure = try
            compute(ParametricProblem(space), formulation)
        catch exception
            exception
        end
        @test attempts[] == expected_attempts
        if expected_attempts == 1
            @test failure === injected
        else
            @test failure isa ErrorException
            message = sprint(showerror, failure)
            @test occursin("retry limit of 3", message)
            @test occursin("3 attempts (0 accepted)", message)
            @test occursin("during build", message)
            @test occursin("unrecoverable input marker", message)
        end
    end
    @test_throws ArgumentError MonteCarlo(inner; options=(on_error=:retry,))
end

@testitem "UQ / negative physical radius / failure and retry records" tags=[:integration] begin
    using Measurements
    attempts=Float64[]
    function physical_problem(radius)
        push!(attempts,radius)
        radius>0 || throw(DomainError(radius,"physical radius must be positive"))
        metal=Material(kind=:conductor,rho=2e-8)
        dielectric=Material(kind=:insulator,rho=1e8,eps_r=3.)
        design=build(CableDesign,"positive-radius",Stack(
            terminal(:core,Region(:metal,Disk(radius),metal)),
            Region(:insulation,Shell(.005),dielectric)))
        CableConstantsProblem(design)
    end
    problem=ParametricProblem(Gridspace{CableConstantsProblem}(physical_problem,
        (Grid(.005,AbsoluteError(.001)),)))
    inner=CableConstantsFormulation()
    negative=(_rng,mean,_sigma)->-mean
    @test_throws DomainError compute(problem,MonteCarlo(inner;trials=2,seed=103,
        distribution=negative,on_error=:fail))
    @test attempts==[-.005]
    empty!(attempts)
    failure=try
        compute(problem,MonteCarlo(inner;trials=2,seed=103,distribution=negative,
            on_error=:retry,max_failures=3,retain_details=true))
    catch error
        error
    end
    @test attempts==fill(-.005,3)
    @test failure isa ErrorException
    @test occursin("3 attempts (0 accepted)",sprint(showerror,failure))
    @test occursin("physical radius must be positive",sprint(showerror,failure))
    # This deterministic control tests conditioning and source records only; no claim
    # is made that accepted retries follow an unconditioned input law.
    calls=Ref(0)
    alternating=(_rng,mean,_sigma)->begin
        calls[]+=1
        isodd(calls[]) ? -mean : mean
    end
    empty!(attempts)
    sampled=compute(problem,MonteCarlo(inner;trials=2,seed=103,distribution=alternating,
        on_error=:retry,max_failures=3,retain_details=true,return_samples=true))
    @test sampled.trial_counts==[2]
    @test attempts==[-.005,.005,-.005,.005,.005]
    @test only(sampled.values).details.data.inputs.interpretation===:nominal_declaration
    @test only(only(sampled.values).details.data.inputs.uncertain_arguments).standard_deviation==.001
    @test length(only(sampled.details.data.failures))==2
    @test only(sampled.details.data.failure_summary).accepted==2
    @test only(sampled.details.data.failure_summary).failed==2
end
