@testitem "UQ / current joint radial scale and spacing / actual line sampling and aggregation" tags=[:integration,:extension] begin
    using Measurements,Random,Statistics
    function problem(scale,spacing)
        scale>0 || throw(DomainError(scale,"radial dimensions must be positive"))
        conductor=Material(kind=:conductor,rho=2e-8,eps_r=1.0,mu_r=1.0,T0=20.0,alpha=.004)
        dielectric=Material(kind=:insulator,rho=1e8,eps_r=3.0)
        designs=[build(CableDesign,"joint-wire-$i",Stack(
            terminal(:core,Region(:metal,Disk(.005scale),conductor)),
            Region(:insulation,Shell(.005scale),dielectric))) for i in 1:2]
        system=build(LineCableSystem,designs,[(0.0,-1.0),(spacing,-1.3)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        return LineParametersProblem(system;frequencies=[50.0,1000.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0),temperature=20.0)
    end
    # Uniform inputs have half-width sqrt(3)*sigma. Even at both worst-case
    # bounds the circular surfaces are separated before clearance adjustment.
    @test hypot(.2-sqrt(3)*.002,.3)-2*.01*(1.1+sqrt(3)*.01)>0
    calls=Tuple{Float64,Float64}[]
    log_inputs=Ref(false)
    builder=(scale,spacing)->begin
        log_inputs[] && !(scale isa Measurement) && !(spacing isa Measurement) && push!(calls,(scale,spacing))
        problem(scale,spacing)
    end
    space=Gridspace{LineParametersProblem}(builder,
        (Grid((1.0,1.1),AbsoluteError(.01)),Grid(.2,AbsoluteError(.002))))
    inner=Formulation(options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    parametric=ParametricProblem(space)
    linear=compute(parametric,LinearError(inner))
    @test length(linear.values)==2
    N=4
    seed=2029
    method=MonteCarlo(inner;trials=N,seed,distribution=:uniform,
        return_samples=true,retain_details=true)
    log_inputs[]=true
    sampled=compute(parametric,method)
    log_inputs[]=false
    @test sampled.trial_counts==[N,N]
    @test all(isempty,sampled.details.data.failures)
    @test length(unique(sampled.point_seeds))==2
    # The builder log contains the inputs that reached actual scalar compute.
    # Reconstruct the exact retained draws through the current sampler, then
    # compare the resulting channel arrays and independently aggregate them.
    recorded=copy(calls)
    for (index,point) in enumerate(LineCableModels.points(space))
        rng=Xoshiro(sampled.point_seeds[index])
        expected=NamedTuple[]
        realized=Tuple{Float64,Float64}[]
        for trial in 1:N
            arguments=LineCableModels.realize_arguments(rng,point,:uniform)
            push!(realized,Tuple(arguments))
            value=compute(problem(arguments...),inner)
            push!(expected,(R=R(value),L=L(value),C=C(value),G=G(value)))
        end
        @test all(pair->count(==(pair),recorded)==1,realized)
        for channel in (:R,:L,:C,:G)
            values=cat((getproperty(value,channel) for value in expected)...;dims=4)
            @test getproperty(sampled.sample_values[index],channel)==values
            average=dropdims(sum(values;dims=4)./N;dims=4)
            observed = getproperty(LineCableModels,channel)(sampled.values[index])
            @test eltype(observed) <: Measurement
            @test nominal.(observed) ≈ average rtol=1e-10 atol=0
            @test uncertainty.(observed) ≈ dropdims(std(values;dims=4);dims=4)
        end
    end
    @test_throws DomainError problem(-.1,.2)
end

@testitem "UQ / Monte Carlo / computed marginal ownership and scientific transport" tags=[:integration,:extension] setup=[TestFixtures] begin
    using Measurements, Statistics, JSON3
    IE=LineCableModels.ImportExport
    design=TestFixtures.coaxial_design()
    space=Gridspace{CableConstantsProblem}(
        temperature -> CableConstantsProblem(design;temperature),
        (Grid((20.0,40.0),AbsoluteError(1.0)),))
    problem=ParametricProblem(space)
    inner=CableConstantsFormulation()
    baseline=compute(problem,MonteCarlo(inner;trials=4,seed=2029,distribution=:uniform,return_samples=true))
    for samples_retained in (false,true), histograms_retained in (false,true), bins in (1,3)
        result=compute(problem,MonteCarlo(inner;trials=4,seed=2029,distribution=:uniform,
            return_samples=samples_retained,return_histograms=histograms_retained,bins))
        @test uncertain(result) === result.values
        @test isempty(details(result).data)
        @test (samples(result) !== nothing) == samples_retained
        @test (histograms(result) !== nothing) == histograms_retained
        @test statistics(result) == statistics(baseline)
        for point in eachindex(result), quantity in (R,L,C,G)
            output=observe(uncertain(result,point),quantity)
            draws=observe(baseline,samples,quantity,point)
            @test eltype(output) <: Measurement
            @test nominal.(output) ≈ vec(mean(draws;dims=2))
            @test uncertainty.(output) ≈ vec(std(draws;dims=2))
            @test all(iszero,uncertainty.(output .- observe(result[point],quantity)))
        end
        record=IE.serialize_value(result)
        @test record["version"] == 2
        restored=IE.deserialize_value(JSON3.read(JSON3.write(record),Dict{String,Any}))
        @test statistics(restored) == statistics(result)
        @test restored.point_seeds == result.point_seeds
        @test restored.trial_counts == result.trial_counts
        for point in eachindex(result), quantity in (R,L,C,G)
            @test nominal.(observe(restored[point],quantity)) == nominal.(observe(result[point],quantity))
            @test uncertainty.(observe(restored[point],quantity)) ≈ uncertainty.(observe(result[point],quantity))
        end
    end
    # V1 is a supported full scientific record, unlike a moments-only publication.
    old=IE.serialize_value(baseline)
    old["version"]=1
    delete!(old,"sources")
    old["points"]=[IE.serialize_value(CableConstants(value.cores,
        nominal.(value.R),nominal.(value.L),nominal.(value.C),nominal.(value.G),
        nominal(value.frequency),value.details)) for value in baseline]
    restored=IE.deserialize_value(JSON3.read(JSON3.write(old),Dict{String,Any}))
    @test statistics(restored) == statistics(baseline)
    @test uncertainty.(first(restored).R) ≈ uncertainty.(first(baseline).R)
    @test uncertain(restored,1) === first(restored)
    rounded=deepcopy(old)
    rounded["point_seeds"]=[1.0,1.0]
    @test_throws r"exact integers" IE.deserialize_value(rounded)

    singleton=compute(problem,MonteCarlo(inner;trials=1,seed=2029,return_samples=true))
    for core in singleton, quantity in (R,L,C,G)
        @test eltype(observe(core,quantity)) <: Measurement
        @test all(iszero,uncertainty.(observe(core,quantity)))
    end
end

@testitem "UQ / linear propagation / affine means variances and covariance" tags=[:integration,:extension] begin
    using Measurements
    struct AffineProblem{T} <: AbstractProblemDefinition
        x::T
        y::T
    end
    LineCableModels.validate(problem::AffineProblem)=problem
    struct AffineFormulation <: AbstractFormulation end
    function LineCableModels.compute(problem::AffineProblem,::AffineFormulation;
            options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
        options = options isa NamedTuple ? ComputationOptions(options) : options
        isempty(options.data) || throw(ArgumentError("affine fixture accepts no execution controls"))
        x,y=problem.x,problem.y
        u,v=2x+3y,x-y
        LineParameters(reshape([complex(u,v)],1,1,1),reshape([complex(v,u)],1,1,1),[50.0])
    end
    space=Gridspace{AffineProblem}(AffineProblem,
        (Grid((2.0,4.0),AbsoluteError(.1)),Grid(3.0,AbsoluteError(.2))))
    result=compute(ParametricProblem(space),LinearError(AffineFormulation()))
    @test length(result)==2
    for (index,x) in enumerate((2.0,4.0))
        u,v=real(only(Z(result[index]))),imag(only(Z(result[index])))
        @test Measurements.value(u)==2x+9
        @test Measurements.value(v)==x-3
        @test Measurements.uncertainty(u)^2 ≈ 4*.1^2+9*.2^2
        @test Measurements.uncertainty(v)^2 ≈ .1^2+.2^2
        @test Measurements.cov(u,v) ≈ 2*.1^2-3*.2^2
        @test only(Y(result[index]))==complex(v,u)
    end
end
