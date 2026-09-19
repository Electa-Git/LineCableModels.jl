@testitem "Engine / shared temperature law reaches scalar, Gridspace, constants and export" tags=[:integration] setup=[FormulaContractModels] begin
    copper = Material(:conductor,1.72e-8,1,1,20,0.004)
    dielectric = Material(:insulator,1e7,2.3,1,20,-0.003;tan_delta=0.025)
    function cable_with(metal, passive)
        build(CableDesign,"temperature-law",terminal(:core,
            core(metal;r=0.005),insulation(passive;t=0.005)))
    end
    function system_with(design)
        build(LineCableSystem,[design,design],[(0.0,-1.0),(0.1,-1.0)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
    end
    design = cable_with(copper,dielectric)
    system = system_with(design)
    problem = LineParametersProblem(system;temperature=80.0,frequencies=[50.0,1000.0],earth_props=homogeneous(rho=100.0))
    declaration=FormulaContractModels.ScaledResistivity()
    selected = Formulation(temperature_dependence=declaration,
        insulation_admittance=:lossy,options=(ideal_transposition=false,))
    reference_design = cable_with(
        Material(:conductor,2copper.rho,copper.eps_r,copper.mu_r),
        Material(:insulator,2dielectric.rho,dielectric.eps_r,dielectric.mu_r;
            tan_delta=dielectric.tan_delta))
    reference_problem = LineParametersProblem(system_with(reference_design);
        temperature=80.0,frequencies=problem.frequencies,earth_props=problem.earth_props)
    identity = Formulation(temperature_dependence=nothing,
        insulation_admittance=:lossy,options=selected.options)
    reference = compute(reference_problem,identity)
    actual = compute(problem,selected)
    @test !isempty(declaration.seen)
    @test actual.Z.values ≈ reference.Z.values rtol=2e-13
    @test actual.Y.values ≈ reference.Y.values rtol=2e-13
    @test details(actual).data.formulations.methods.temperature_dependence.identifier === :ScaledResistivity
    grid = Formulation(temperature_dependence=Grid((declaration,nothing)),
        insulation_admittance=:lossy,options=selected.options)
    results = compute(problem,grid)
    @test results[1].Z.values == actual.Z.values
    @test results[1].Y.values == actual.Y.values
    unchanged = compute(problem,identity)
    @test results[2].Z.values == unchanged.Z.values
    @test results[2].Y.values == unchanged.Y.values
    @test !isapprox(results[1].Z.values,results[2].Z.values;rtol=1e-5)
    constants = compute(CableConstantsProblem(design;temperature=80.0),
        CableConstantsFormulation(temperature_dependence=declaration,insulation_admittance=:lossy))
    reference_constants = compute(CableConstantsProblem(reference_design;temperature=80.0),
        CableConstantsFormulation(temperature_dependence=nothing,insulation_admittance=:lossy))
    for request in (R,L,C,G)
        @test request(constants) ≈ request(reference_constants) rtol=2e-13
    end
    const IE = LineCableModels.ImportExport
    exported = only(LineCableModels.PSCAD._pscad_components(design,50.0,selected,80.0))
    expected = only(LineCableModels.PSCAD._pscad_components(reference_design,50.0,identity,80.0))
    @test exported.conductor.material.rho == expected.conductor.material.rho
    @test exported.dielectric.shunt_conductance ≈ expected.dielectric.shunt_conductance
    @test exported.dielectric.shunt_capacitance ≈ expected.dielectric.shunt_capacitance
    hot = LineParametersProblem(system;temperature=250.0,frequencies=[50.0],earth_props=problem.earth_props)
    @test_throws DomainError compute(hot,Formulation())
    @test all(isfinite,compute(hot,identity).Z)
    @test all(isfinite,compute(hot,selected).Z)
end

@testitem "UQ / current AC cable / actual sampling and independently recomputed statistics" tags=[:integration] setup=[TestFixtures] begin
    using Random
    design=TestFixtures.coaxial_design()
    temperatures=(20.0,60.0)
    space=Gridspace{CableConstantsProblem}(t->CableConstantsProblem(design;temperature=t),
        (Grid(temperatures,AbsoluteError(1.0)),))
    problem=ParametricProblem(space)
    N=8
    seed=2027
    sampled=compute(problem,MonteCarlo(CableConstantsFormulation();trials=N,seed,
        distribution=:uniform,return_samples=true,return_histograms=true,retain_details=true))
    @test sampled.trial_counts==[N,N]
    @test length(unique(sampled.point_seeds))==2
    @test all(isempty,sampled.details.data.failures)
    # Rebuild each retained draw through scalar compute; statistics below are
    # arithmetic checks on these samples, not distribution-accuracy claims.
    for (index,point) in enumerate(LineCableModels.points(space))
        rng=Xoshiro(sampled.point_seeds[index])
        expected=[compute(CableConstantsProblem(design;
            temperature=only(LineCableModels.realize_arguments(rng,point,:uniform))),
            CableConstantsFormulation()) for _ in 1:N]
        for quantity in (:R,:L,:C,:G)
            @test vec(getproperty(sampled.sample_values[index],quantity))==
                [only(getproperty(value,quantity)) for value in expected]
        end
    end
    quantile7(sorted,p)=begin
        position=1+(length(sorted)-1)*p
        i=floor(Int,position);fraction=position-i
        i==length(sorted) ? last(sorted) : (1-fraction)*sorted[i]+fraction*sorted[i+1]
    end
    for point in eachindex(temperatures),quantity in (:R,:L,:C,:G)
        draws=vec(getproperty(sampled.sample_values[point],quantity))
        summary=only(getproperty(sampled.stats[point],quantity))
        density=only(getproperty(sampled.histogram_values[point],quantity))
        @test length(draws)==N
        expected_mean=sum(draws)/N
        expected_std=sqrt(sum(x->(x-expected_mean)^2,draws)/(N-1))
        sorted=sort(draws)
        @test summary.mean ≈ expected_mean rtol=1e-10 atol=0
        @test summary.std ≈ expected_std rtol=1e-10 atol=eps(maximum(abs,draws))
        @test summary.min==first(sorted) && summary.max==last(sorted)
        @test summary.q05 ≈ quantile7(sorted,.05) rtol=1e-10 atol=0
        @test summary.median ≈ quantile7(sorted,.5) rtol=1e-10 atol=0
        @test summary.q95 ≈ quantile7(sorted,.95) rtol=1e-10 atol=0
        counts=zeros(Int,length(density.density))
        for value in draws
            index=min(searchsortedlast(density.edges,value),length(counts))
            1<=index<=length(counts) || error("sample outside declared histogram support")
            counts[index]+=1
        end
        @test density.density ≈ counts./(N.*diff(density.edges)) rtol=1e-10 atol=0
        @test sum(density.density.*diff(density.edges)) ≈ 1.0 rtol=1e-10
        @test getproperty(sampled.values[point],quantity)[1] ≈ summary.mean rtol=1e-10 atol=0
    end
    for point in eachindex(temperatures)
        arrays=values(sampled.sample_values[point]);models=values(sampled.histogram_values[point])
        for i in eachindex(arrays),j in eachindex(arrays)
            i==j && continue
            @test arrays[i] !== arrays[j]
            @test only(models[i]).density !== only(models[j]).density
        end
        @test all(iszero,sampled.sample_values[point].G)
    end
end
