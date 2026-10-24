@testitem "UQ / redrawing infeasible draws preserves accepted trials, source identities and seeded replay" tags=[:integration, :measurements] setup=[TestFixtures] begin
    using Measurements
    using Statistics
    design = TestFixtures.coaxial_design()
    attempted_temperatures = Float64[]
    space = Gridspace{CableConstantsProblem}(
        temperature -> begin
            push!(attempted_temperatures, temperature)
            # The nominal problems of both points build first. The first two draws fail.
            length(attempted_temperatures) in 3:4 && throw(DomainError(
                temperature, "injected rejection before resolving the cable"))
            CableConstantsProblem(design; temperature)
        end,
        (Grid((20.0, 40.0), AbsoluteError((0.5,))),),
    )
    inner = CableConstantsFormulation()
    formulation = MonteCarlo(inner; trials=4, seed=71, distribution=:uniform,
        return_samples=true, return_histograms=true, bins=2, options=(max_failures=3,))
    sampled = compute(ParametricProblem(space), formulation; options=(retain_details=true,))
    @test length(sampled) == 2
    @test sampled.trial_counts == [4, 4]
    @test length(attempted_temperatures) == 12
    @test attempted_temperatures[[1,2]] == [20.,40.]
    @test length(unique(sampled.point_seeds)) == 2
    @test all(isconcretetype, (eltype(sampled.values), eltype(sampled.stats),
        eltype(sampled.sample_values), eltype(sampled.histogram_values)))
    @test length.(sampled.details.data.trials) == [4, 4]
    @test keys(sampled.details.data) == (:trials,)
    for point in 1:2
        accepted_temperatures = attempted_temperatures[point==1 ? (5:8) : (9:12)]
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
    recorded = copy(attempted_temperatures)
    empty!(attempted_temperatures)
    replay = compute(ParametricProblem(space), formulation; options=(retain_details=true,))
    @test attempted_temperatures == recorded
    @test samples(replay) == samples(sampled)
    @test replay.point_seeds == sampled.point_seeds

end

@testitem "UQ / redrawing does not hide other errors or run past its limit" tags=[:integration, :uq] setup=[TestFixtures] begin
    using Measurements
    inner = CableConstantsFormulation()
    design = TestFixtures.coaxial_design()
    for (injected, expected_attempts) in (
        (ArgumentError("programming error marker"), 2),
        (DomainError(-1, "unrecoverable input marker"), 4),
    )
        attempts = Ref(0)
        # The nominal problem builds, and every draw throws.
        space = Gridspace{CableConstantsProblem}(
            _ -> begin
                attempts[] += 1
                attempts[] > 1 && throw(injected)
                CableConstantsProblem(design)
            end,
            (Grid((1.0,)),),
        )
        formulation = MonteCarlo(inner; trials=2, seed=71, options=(max_failures=3,))
        failure = try
            compute(ParametricProblem(space), formulation; options=(retain_details=true,))
        catch exception
            exception
        end
        @test attempts[] == expected_attempts
        if expected_attempts == 2
            @test failure === injected
        else
            # The point fails, and a study with no computed point throws.
            @test failure isa ErrorException
            message = sprint(showerror, failure)
            @test occursin("all 1 points of the space failed", message)
            # The point's `DomainError` holds its uncertain arguments: none.
            @test occursin("DomainError with NamedTuple[]:", message)
            @test occursin("discarded 3 infeasible draws in a row", message)
            @test occursin("max_failures", message)
            @test occursin("unrecoverable input marker", message)
        end
    end
end

@testitem "UQ / negative physical radius / redraws, their limit and source records" tags=[:integration, :measurements] begin
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
    # Point 1 has a feasible nominal radius and negative draws only, so it fails. Point 2
    # draws its nominal radius.
    negative=(_rng,mean,_sigma)->mean<.0055 ? -mean : mean
    two=ParametricProblem(Gridspace{CableConstantsProblem}(physical_problem,
        (Grid((.005,.006),AbsoluteError(.001)),)))
    sampled=compute(two,MonteCarlo(inner;trials=2,seed=103,distribution=negative,
        max_failures=3);options=(retain_details=true,))
    @test attempts==[.005,.006,-.005,-.005,-.005,.006,.006]
    @test sampled.trial_counts==[2]
    gridpoint,failure=only(sampled.failures)
    @test gridpoint.problem_index==1
    @test failure isa DomainError
    @test only(failure.val).standard_deviation==.001
    @test occursin("discarded 3 infeasible draws in a row",sprint(showerror,failure))
    @test occursin("physical radius must be positive",sprint(showerror,failure))
    # This deterministic control tests the redraw and source records only. No claim
    # is made that accepted draws follow an unconditioned input law.
    calls=Ref(0)
    alternating=(_rng,mean,_sigma)->begin
        calls[]+=1
        isodd(calls[]) ? -mean : mean
    end
    empty!(attempts)
    sampled=compute(problem,MonteCarlo(inner;trials=2,seed=103,distribution=alternating,
        max_failures=3,return_samples=true);options=(retain_details=true,))
    @test sampled.trial_counts==[2]
    @test attempts==[.005,-.005,.005,-.005,.005]
    @test only(sampled.values).details.data.inputs.interpretation===:nominal_declaration
    @test only(only(sampled.values).details.data.inputs.uncertain_arguments).standard_deviation==.001
end
